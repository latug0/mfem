// Compile with: make cas_point_contact_p
//
//Version parallèle
//Cas test avec solution analytique. Plaque 2D avec point de contact
//d'intensité Force en (0,0).
//Calcul de la déformation en élasticité linéaire. 

//Calcul des erreurs en normes H1 L2 et énergie

//ref: http://e6.ijs.si/medusa/wiki/index.php/point_contact

#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

static constexpr double Force = -1.;
static constexpr double e = 1000.;	//module de young
static constexpr double nu = 0.25;	//coef de poisson
static constexpr double lambda = e*nu/((1.+nu)*(1.-2.*nu));
static constexpr double mu = e/(2.*(1.+nu));	//coef de lamé
//Solution exacte
void sol_exact(const Vector &, Vector &);
void Grad_Exact(const Vector &, DenseMatrix &);
void Stress_excat(const Vector &, Vector &);
void Strain_exacte(const Vector &, Vector &);

//Erreur en norme energie
double computeenergynorm(ParGridFunction &,
			 Coefficient &, Coefficient &);

void elasticy_mat(ElementTransformation &,const IntegrationPoint &, 
		  int, Coefficient &, Coefficient &, DenseMatrix &);
void computestrain(ElementTransformation &,const IntegrationPoint &,
		   ParGridFunction &, int,  Vector &);

//Erreur en norme H1
double computeh1norm(ParGridFunction &);

int main(int argc, char *argv[])
{

  // initialize mpi.
  int num_procs, myid;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  // parse command-line Options.
  int order=2;
  const char *mesh_file = "carre.msh";
  bool static_cond = false;
  bool amg_elast = 0;
  bool reorder_space = false;
  bool solver=true;
  int ref_levels = 8;
  OptionsParser args(argc, argv);
  args.AddOption(&mesh_file, "-m", "--mesh",
 		 "mesh file to use.");
  args.AddOption(&order, "-o", "--order",
 		 "finite element order (polynomial degree).");
  args.AddOption(&amg_elast, "-elast", "--amg-for-Elasticity", "-sys",
 		 "--amg-for-systems",
 		 "use the special amg Elasticity solver (gm/ln approaches), "
 		 "or standard amg for systems (unknown approach).");
  args.AddOption(&static_cond, "-sc", "--static-condensation", "-no-sc",
 		 "--no-static-condensation", "enable static condensation.");
  args.AddOption(&reorder_space, "-nodes", "--by-nodes", "-vdim", "--by-vdim",
 		 "use bynodes ordering of Vector space instead of byvdim");
  args.AddOption(&solver, "-sol", "--itératif", "-direct",
		 "--solver itératif", "solver direct.");
  args.AddOption(&ref_levels, "-ref", "--num_ref", "nombre de rafinement de maillage");

  args.Parse();
  if (!args.Good())
    {
      if (myid == 0)
 	{
 	  args.PrintUsage(cout);
 	}
      MPI_Finalize();
      return 1;
    }
  if (myid == 0)
    {
      args.PrintOptions(cout);
    }
  //lecture du malliage
  Mesh *mesh = new Mesh(mesh_file, 1, 1);
  int dim = mesh->Dimension();
  //  refine the mesh to increase the resolution. in this example we do
  //    'ref_levels' of uniform refinement. we choose 'ref_levels' to be the
  //    largest number that gives a final mesh with no more than 5,000
  //    elements.
  //   int ref_levels = (int)floor(log(5000./mesh->Getne())/log(2.)/dim);
  for (int l = 0; l < ref_levels; l++)
    {
      mesh->UniformRefinement();
    }
	  cout<<"Nb de dofs: "<<mesh->GetNE()<<endl;   
  //define parallel mesh by a partitioning of the serial mesh.
  ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);

  //  define a finite element space on the mesh. here we use Vector finite
  //    elements, i.e. dim copies of a scalar finite element space. the Vector
  //    dimension is specified by the last argument of the ParFiniteElementSpace
  //    constructor. for nurbs meshes, we use the (degree elevated) nurbs space
  //    associated with the mesh nodes.
  FiniteElementCollection *fec;
  ParFiniteElementSpace *fespace;
  const bool use_nodal_fespace = pmesh->NURBSext && !amg_elast;
  if (use_nodal_fespace)
    {
      fec = NULL;
      fespace = (ParFiniteElementSpace *)pmesh->GetNodes()->FESpace();
    }
  else
    {
      fec = new H1_FECollection(order, dim);
      if (reorder_space)
 	{
 	  fespace = new ParFiniteElementSpace(pmesh, fec, dim, Ordering::byVDIM);
 	}
      else
 	{
 	  fespace = new ParFiniteElementSpace(pmesh, fec, dim,  Ordering::byVDIM);
 	}
    }
  HYPRE_Int Size = fespace->GlobalTrueVSize();
  if (myid == 0)
    {
      cout << "number of finite element unknowns: " << Size << endl
	   << "assembling: " << flush;
    }
  // determine the list of true (i.e. conforming) essential boundary dofs.
  //    in this example, the boundary coGlobalTrueVSizenditions are defined by marking only
  //    boundary attribute 1 from the mesh as essential and converting it to a
  //    list of true dofs.

  // list of true dofs : define (here) dirichlet conditions
  Array<int> ess_tdof_list;
  Array<int> ess_bdr(pmesh->bdr_attributes.Max());
  ess_bdr = 1;		//on sélectionne tous les bords
  fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
  // Set up the linear form b(.) which corresponds to the right-hand side of
  //    the fem linear system. in this case, b_i equals the boundary integral
  //    of f*phi_i where f represents a "pull down" Force on the neumann part
  //    of the boundary and phi_i are the basis Functions in the finite element
  //    fespace. the Force is defined by the VectorArrayCoefficient object f,
  //    which is a Vector of Coefficient objects. the fact that f is non-zero
  //    on boundary attribute 2 is indicated by the use of piece-wise Constants
  //    Coefficient for its last component.
  VectorArrayCoefficient f(dim);
  for (int i = 0; i < dim-1; i++)
    {
      f.Set(i, new ConstantCoefficient(0.0));
    }
  ParLinearForm *b = new ParLinearForm(fespace);
  b->Assemble();
  if (myid == 0)
    {
      cout << "r.h.s. ... " << flush;
    }

  // define the solution Vector x as a finite element Grid Function
  //    corresponding to fespace. initialize x with initial guess of zero,
  //    which satisfies the boundary conditions.
  ParGridFunction x(fespace);
  VectorFunctionCoefficient boundary_dirichlet_coef(dim, sol_exact);
  x = 0.0;
  // to use if there are different dirichlet conditions.
  // beware, the values of dirichlet boundary conditions are Set here !
  //projète solution exacte sur les bords
  x.ProjectBdrCoefficient(boundary_dirichlet_coef, ess_bdr);	
  // Set up the bilinear form a(.,.) on the finite element space
  //    corresponding to the linear Elasticity integrator with piece-wise
  //    Constants Coefficient lambda and mu.
  ConstantCoefficient mu_func(mu);
  ConstantCoefficient lambda_func(lambda);
  ParBilinearForm *a = new ParBilinearForm(fespace);
  BilinearFormIntegrator *integ = new ElasticityIntegrator(lambda_func, mu_func);
  a->AddDomainIntegrator(integ);
  // assemble the parallel bilinear form and the corresponding linear
  //     system, applying any necessary transformations such as: parallel
  //     assembly, eliminating boundary conditions, applying conforming
  //     constraints for non-conforming amr, static condensation, etc.
  if (myid == 0) { cout << "matrix ... " << flush; }
  if (static_cond) { a->EnableStaticCondensation(); }
  a->Assemble();
  HypreParMatrix A;
  Vector B, X;
  a->FormLinearSystem(ess_tdof_list, x, *b, A, X, B);
  if (myid == 0)
    {
      cout << "done." << endl;
      cout << "Size of linear system: " << A.GetGlobalNumRows() << endl;
    }
  // 13. define and apply a parallel pcg solver for a x = b with the boomeramg
  //     preconditioner from hypre.
  if(solver){
    HypreBoomerAMG *amg = new HypreBoomerAMG(A);
    if (amg_elast && !a->StaticCondensationIsEnabled())
      {
 	amg->SetElasticityOptions(fespace);
      }
    else
      {
 	amg->SetSystemsOptions(dim, reorder_space);
      }
    HyprePCG *pcg = new HyprePCG(A);
    pcg->SetTol(1e-20);
    pcg->SetMaxIter(5000);
    pcg->SetPrintLevel(2);
    pcg->SetPreconditioner(*amg);
    pcg->Mult(B, X);
    delete pcg;
    delete amg;
  }
  else{
    cout<<"solver direct non implémenté"<<endl;
  }   


  // 12. recover the solution as a finite element Grid Function.
  a->RecoverFEMSolution(X, *b, x);

  // 13. for non-nurbs meshes, make the mesh curved based on the finite element
  //     space. this means that we define the mesh elements through a fespace
  //     based transformation of the reference element. this allows us to save
  //     the displaced mesh as a curved mesh when using high-order finite
  //     element displacement field. we assume that the initial mesh (read from
  //     the file) is not higher order curved mesh c<ompared to the chosen fe
  //     space.
  if (!use_nodal_fespace)
    {
      pmesh->SetNodalFESpace(fespace);
    }

  // compute errors
HYPRE_Int global_dofs_ = fespace->GlobalTrueVSize();
	  ParGridFunction zero(fespace);
	  zero = 0.0;
      double ener_reff = computeenergynorm(zero, lambda_func, mu_func);
  double h1_error = computeh1norm(x);
  double ener_error = computeenergynorm(x, lambda_func, mu_func);
  VectorFunctionCoefficient sol_exact_coef(dim, sol_exact);
  double l2_error = x.ComputeL2Error(sol_exact_coef);

  double h = pmesh->GetElementSize(1);
  if (myid == 1)
    {
      cout<<"erreur en norme l2: "<<l2_error<<endl;
      cout<<"erreur en norme énergie: "<<ener_error/ener_reff<<endl;
      cout<<"erreur en norme h1: "<<h1_error<<endl;  
      cout<<"taille de maille: "<<h<<endl; 
	  cout<<"Nb de dofs: "<<fespace->GetNDofs()<<endl;   
      cout << "numbers of elements: " << pmesh->GetNE() <<endl;
    }
  //  free the used memory.
  delete a;
  delete b;
  if (fec)
    {
      delete fec;
    }
  delete mesh;
  delete pmesh;
  MPI_Finalize();
  return 0;
}

//===================== solution exacte =====================
void sol_exact(const Vector &x, Vector &u)
{
  double r2 = x(1)*x(1)+x(0)*x(0);
  double pi = M_PI;
  u(0) = -Force/(4*pi*mu)*(2*x(0)*x(1)/r2 + 2*mu/(lambda+mu)*atan2(x(1),x(0)));
  u(1) = -Force/(4*pi*mu)*((x(1)*x(1)-x(0)*x(0))/r2 - 
 			   (lambda+2*mu)*log(r2)/(lambda+mu));
}
//===================== grad exacte =====================
void grad_exact(const Vector &x, DenseMatrix &grad)
{
  double r2 = x(1)*x(1)+x(0)*x(0);
  double pi = M_PI;
  grad(0,0) = -Force/(4*pi*mu)*(2*x(1)*(x(1)*x(1)-x(0)*x(0))/(r2*r2) -
				2*mu/(lambda+mu)*x(1)/r2);
  grad(0,1) = -Force/(4*pi*mu)*(2*x(0)*(x(0)*x(0)-x(1)*x(1))/(r2*r2) + 
				2*mu/(lambda+mu)*x(0)/r2);
  grad(1,1) = -Force/(4*pi*mu)*(4*x(1)*x(0)*x(0)/(r2*r2) - 
 				(lambda+2*mu)/(lambda+mu)*2*x(1)/r2);
  grad(1,0) = -Force/(4*pi*mu)*(-4*x(1)*x(1)*x(0)/(r2*r2) - 
 				(lambda+2*mu)/(lambda+mu)*2*x(0)/r2);
}
//===================== Stress exacte =====================
void Stress_exacte(const Vector &x, Vector &stress)
{
  double r2 = x(1)*x(1)+x(0)*x(0);
  double pi = M_PI;
  stress(0) = 2*Force/pi*x(0)*x(0)*x(1)/(r2*r2);
  stress(1) = 2*Force/pi*x(1)*x(1)*x(1)/(r2*r2);
  stress(2) = 2*Force/pi*x(0)*x(1)*x(1)/(r2*r2);
}
//===================== Strain exacte =====================
void Strain_exacte(const Vector &x, Vector &strain)
{
  int dim = x.Size();
  DenseMatrix grad(dim,dim);
  double r2 = x(1)*x(1)+x(0)*x(0);
  double pi = M_PI;
  grad(0,0) = -Force/(4*pi*mu)*(2*x(1)*(x(1)*x(1)-x(0)*x(0))/(r2*r2) -
				2*mu/(lambda+mu)*x(1)/r2);
  grad(0,1) = -Force/(4*pi*mu)*(2*x(0)*(x(0)*x(0)-x(1)*x(1))/(r2*r2) + 
				2*mu/(lambda+mu)*x(0)/r2);
  grad(1,1) = -Force/(4*pi*mu)*(4*x(1)*x(0)*x(0)/(r2*r2) - 
				(lambda+2*mu)/(lambda+mu)*2*x(1)/r2);
  grad(1,0) = -Force/(4*pi*mu)*(-4*x(1)*x(1)*x(0)/(r2*r2) - 
				(lambda+2*mu)/(lambda+mu)*2*x(0)/r2);
  strain.SetSize(dim*(dim+1)/2);
  strain(0)=grad(0,0);
  strain(1)=grad(1,1);
  if(dim==2){
    strain(2)=0.5*(grad(1,0)+grad(0,1));
  }
  else if(dim==3){
    strain(2)=grad(2,2);
    strain(3)=0.5*(grad(1,0)+grad(0,1));
    strain(4)=0.5*(grad(2,0)+grad(0,2));
    strain(5)=0.5*(grad(2,1)+grad(1,2));
  }
  else{
    cout<<"Dimention not suported"<<endl;}
}


//===================== erreur en norme h1 =====================
double computeh1norm(ParGridFunction &x){
  ParFiniteElementSpace *fes = x.ParFESpace();
  int dim = fes->GetMesh()->SpaceDimension();
  MatrixFunctionCoefficient grad_exact_coef (dim, grad_exact);
  ElementTransformation *trans;
  DenseMatrix grad, gradh;
  double error_local = 0.0, error_global = 0.0;
  Array<int> udofs;
  for (int i = 0; i < fes->GetNE() ; i++)
    {
      const FiniteElement *fe = fes->GetFE(i);
      const int order = 2*fe->GetOrder() + 3;   //<----------
      const IntegrationRule *ir = &(IntRules.Get(fe->GetGeomType(), order));
      trans = fes->GetElementTransformation(i);
      const int dof = fe->GetDof();
      const int tdim = dim*(dim+1)/2; // num. entries in a symmetric tensor

      DenseMatrix dshape(dof, dim);
      DenseMatrix gh(dim, dim),gradh (dim, dim),grad(dim,dim);
      fes->GetElementVDofs(i, udofs);
      DenseMatrix loc_data(dof, dim);

      for (int s=0 ; s<dim ; s++)
 	{
 	  Array<int> udofs_tmp(dof);
 	  udofs_tmp = 0.;
 	  for(int j=0 ; j<dof ; j++){
 	    udofs_tmp[j] = udofs[j+dof*s];}

 	  Vector loc_data_tmp;
 	  x.GetSubVector(udofs_tmp, loc_data_tmp);

 	  for (int j=0 ; j<dof ; j++)
 	    loc_data(j,s) = loc_data_tmp(j);
 	}
      for (int j = 0; j < ir->GetNPoints(); j++)
 	{
	  const IntegrationPoint &ip = ir->IntPoint(j);
 	  trans->SetIntPoint(&ip);
 	  double w = trans->Weight() * ip.weight;
 	  fe->CalcDShape(ip, dshape);
 	  MultAtB(loc_data, dshape, gh);
	  Mult(gh, trans->InverseJacobian(), gradh);
 	  grad_exact_coef.Eval(grad,*trans,ip);
 	  grad -= gradh;
 	  error_local += w * grad.FNorm2();
 	}			
    }
  MPI_Reduce(&error_local, &error_global, 1, MPI_DOUBLE, MPI_SUM, 0,
 	     MPI_COMM_WORLD);
  if(error_global>0.0){
    return sqrt(error_global);}
  else{
    //cout<<"Negative H1 error"<<endl;
    return 0;
  }
}

//==============erreur en norme energie ===================
double computeenergynorm(ParGridFunction &x,
 			 Coefficient &lambdah, Coefficient &muh)
{
  ParFiniteElementSpace *fes = x.ParFESpace();
  int dim = fes->GetMesh()->SpaceDimension();
  const int tdim = dim*(dim+1)/2; // num. entries in a symmetric tensor
  VectorFunctionCoefficient Strain_exacte_coef(tdim,Strain_exacte);
  ElementTransformation *Trans;

  double energy_local = 0.0, energy_global = 0.0;
  for (int i = 0; i < fes->GetNE() ; i++)
    {
      const FiniteElement *fe = fes->GetFE(i);
      const int order = 2*fe->GetOrder()+3; // <----------
      const IntegrationRule *ir = &IntRules.Get(fe->GetGeomType(), order);
      const int dof = fe->GetDof();
      Trans = fes->GetElementTransformation(i);
      Vector stressh(tdim), strainh(tdim);	//approché
      Vector stress(tdim), strain(tdim);	//exacte
      DenseMatrix Ch;
      for (int j = 0; j < ir->GetNPoints(); j++)
	{
	  const IntegrationPoint &ip = ir->IntPoint(j);
	  Trans->SetIntPoint(&ip);
	  double w = Trans->Weight() * ip.weight;

	  //======= strain vectors ========
	  computestrain(*Trans, ip, x, i, strainh);
	  Strain_exacte_coef.Eval(strain,*Trans,ip);
	  //======= Stress vectors ========
	  elasticy_mat(*Trans,ip,dim,lambdah,muh,Ch);
	  Ch.Mult(strainh,stressh);	//approx
	  Ch.Mult(strain,stress);	//exacte

	  strain -= strainh;
	  stress -= stressh;

	  double pdc=0.0;
	  for (int k = 0; k< dim; k++)
	    pdc += strain(k)*stress(k);

	  for (int k = dim; k < dim*(dim+1)/2; k++)
	    pdc += 2*strain(k)*stress(k);

	  energy_local += w * pdc;
 	}
    }
  MPI_Allreduce(&energy_local, &energy_global, 1, MPI_DOUBLE, MPI_SUM,
 	     MPI_COMM_WORLD);
  if(energy_global>0.0){
    return sqrt(energy_global);}
  else{
    return 0; //cout<<"Negative Energy error"<<endl;
  }
}

//===================== matrice élasticité =====================
void elasticy_mat(ElementTransformation &t,const IntegrationPoint &ip, 
 		  int dim, Coefficient &lambda, Coefficient &mu_func, DenseMatrix &c){
  double m = mu_func.Eval(t, ip);
  double l = lambda.Eval(t, ip);
  c.SetSize(dim*(dim+1)/2,dim*(dim+1)/2);
  c = 0.;
  for (int k = 0; k< dim; k++)
    {
      // extra-diagonal terms
      for (int lh = 0; lh< dim; lh++)
 	c(k,lh) = l;
      // diagonal terms
      c(k,k) = l+2*m;
    }
  // diagonal terms
  for (int k = dim; k < dim*(dim+1)/2; k++)
    c(k,k) = m;
}
//===================== déformation =====================
void computestrain(ElementTransformation &t,const IntegrationPoint &ip,
 		   ParGridFunction &x, int elem,  Vector &stress){
  ParFiniteElementSpace *fes = x.ParFESpace();
  Array<int> udofs;
  const FiniteElement *fe = fes->GetFE(elem);
  const int dof = fe->GetDof();
  const int dim = fe->GetDim();
  const int tdim = dim*(dim+1)/2; // num. entries in a symmetric tensor
  DenseMatrix dshape(dof, dim);
  DenseMatrix gh(dim, dim),grad (dim, dim);
  fes->GetElementVDofs(elem, udofs);
  DenseMatrix loc_data(dof, dim);
  for (int s=0 ; s<dim ; s++)
    {
      Array<int> udofs_tmp(dof);
      udofs_tmp = 0.;
      for(int j=0 ; j<dof ; j++){
 	udofs_tmp[j] = udofs[j+dof*s];}
      Vector loc_data_tmp;
      x.GetSubVector(udofs_tmp, loc_data_tmp);

      for (int j=0 ; j<dof ; j++)
 	loc_data(j,s) = loc_data_tmp(j);
    }
  double w = t.Weight() * ip.weight;
  fe->CalcDShape(ip, dshape);
  MultAtB(loc_data, dshape, gh);
  Mult(gh, t.InverseJacobian(), grad);
  stress(0)=grad(0,0);
  stress(1)=grad(1,1);
  if(dim==2){
    stress(2)=0.5*(grad(1,0)+grad(0,1));
  }
  else if(dim==3){
    stress(2)=grad(2,2);
    stress(3)=0.5*(grad(1,0)+grad(0,1));
    stress(4)=0.5*(grad(2,0)+grad(0,2));
    stress(5)=0.5*(grad(2,1)+grad(1,2));
  }
  else{
    cout<<"dimention not suported"<<endl;}
}

