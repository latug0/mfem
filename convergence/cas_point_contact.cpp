// Compile with: make cas_point_contact
//
//Cas test avec solution analytique. Plaque 2D avec point de contact
//d'intensité Force en (0,0).
//Calcul de la déformation en élasticité linéaire. 

//Ref: http://e6.ijs.si/medusa/wiki/index.php/Point_contact

//Calculs avec plusieurs maillages uniforméments raffinés possible. 
//Calcul des erreurs en normes H1 L2 et énergie
//Calcul des pentes de convergences avec ces trois normes.

#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <chrono>
#include <ctime>

using namespace std;
using namespace mfem;
static constexpr double Force = -1.;
static constexpr double E = 1000.;	//Module de young
static constexpr double nu = 0.25;	//Coef de poisson
static constexpr double lambda = E*nu/((1.+nu)*(1.-2.*nu));
static constexpr double mu = E/(2.*(1.+nu));	//coef de Lamé
//Solution exacte
void sol_exact(const Vector &, Vector &);
void Grad_Exact(const Vector &, DenseMatrix &);
void Stress_exacte(const Vector &, Vector &);
void Strain_exacte(const Vector &, Vector &);
void Mat_zero(const Vector &, DenseMatrix &);
void Vec_zero(const Vector &, Vector &);
double ux_exact(const Vector &);
double uy_exact(const Vector &);
void gradux_exact(const Vector &, Vector &);
void graduy_exact(const Vector &, Vector &);
//Chargement
double F(const Vector &);

//Erreur en norme H1
double ComputeH1Norm(GridFunction &, MatrixFunctionCoefficient &);

//Erreur en norme energy
double ComputeEnergyNorm(GridFunction &,
			 Coefficient &, Coefficient &,
			 VectorFunctionCoefficient &);

void Elasticy_mat(ElementTransformation &,const IntegrationPoint &, int,
		  Coefficient &, Coefficient &, DenseMatrix &);
void ComputeStrain(ElementTransformation &,const IntegrationPoint &,
		   GridFunction &, int,  Vector &);

int main(int argc, char *argv[])
{
  DenseMatrix slope_l2, slope_ener, slope_grad;
  double err_tmp_ener = 0., err_tmp_l2=0., err_tmp_grad=0.;
  double h_tmp = 0.;
  int iter = 0;

  // Parse command-line options.
  const char *mesh_file = "carre.msh";
  int order=1;
  bool static_cond = false;
  int rep=8;
  bool iterative=true;

  OptionsParser args(argc, argv);
  args.AddOption(&mesh_file, "-m", "--mesh",
		 "Mesh file to use.");
  args.AddOption(&order, "-o", "--order",
		 "Finite element order (polynomial degree).");
  args.AddOption(&rep, "-r", "--Repetition",
		 "Nombre de raffinement de maillage.");
  args.AddOption(&static_cond, "-sc", "--static-condensation", "-no-sc",
		 "--no-static-condensation", "Enable static condensation.");
  args.AddOption(&iterative, "-it", "--iterative", "-di",
		 "--direct", "Enable direct or iterative solver.");
  args.Parse();
  if (!args.Good())
    {
      args.PrintUsage(cout);
      return 1;
    }
  args.PrintOptions(cout);

  slope_ener.SetSize(rep-1,3); slope_l2.SetSize(rep-1,3),
				 slope_grad.SetSize(rep-1,3);
  string const err_energy("error_contact_ordretest.txt");
  ofstream err_energy_flux(err_energy.c_str());

  if (!err_energy_flux.is_open()) {
    cout << "Problem in openning file" << endl;
    exit(0);
  }
  else{
    cout <<"File openning"<<endl;
    //Lecture du malliage
    std::chrono::duration<double> time1;
    auto start1 = std::chrono::system_clock::now();

    Mesh *mesh = new Mesh(mesh_file, 1, 1);
    int dim = mesh->Dimension();
    for (int ref_levels=1; ref_levels<rep; ref_levels++){ 

      // 5. Define a finite element space on the mesh. Here we use vector finite
      //    elements, i.e. dim copies of a scalar finite element space. The vector
      //    dimension is specified by the last argument of the FiniteElementSpace
      //    constructor. For NURBS meshes, we use the (degree elevated) NURBS space
      //    associated with the mesh nodes.
      FiniteElementCollection *fec;
      FiniteElementSpace *fespace;
      fec = new H1_FECollection(order, dim);
      fespace = new FiniteElementSpace(mesh, fec, dim);
      cout << "Numbers of elements: " << mesh->GetNE() <<endl;
      cout << "Number of finite element unknowns: " << fespace->GetTrueVSize()
	   << endl << "Assembling: "<< flush;

      // 6. Determine the list of true (i.e. conforming) essential boundary dofs.
      //    In this example, the boundary conditions are defined by marking only
      //    boundary attribute 1 from the mesh as essential and converting it to a
      //    list of true dofs.

      // List of True DoFs : Define (here) Dirichlet conditions
      Array<int> ess_tdof_list;
      Array<int> ess_bdr(mesh->bdr_attributes.Max());
      ess_bdr = 1;		//On sélectionne tous les bords
      fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

      // 7. Set up the linear form b(.) which corresponds to the right-hand side of
      //    the FEM linear system. In this case, b_i equals the boundary integral
      //    of f*phi_i where f represents a "pull down" force on the Neumann part
      //    of the boundary and phi_i are the basis functions in the finite element
      //    fespace. The force is defined by the VectorArrayCoefficient object f,
      //    which is a vector of Coefficient objects. The fact that f is non-zero
      //    on boundary attribute 2 is indicated by the use of piece-wise constants
      //    coefficient for its last component.

      VectorArrayCoefficient f(dim);
      for (int i = 0; i < dim-1; i++)
	{
	  f.Set(i, new ConstantCoefficient(0.0));
	}
      LinearForm *b = new LinearForm(fespace);
      cout << "r.h.s. ... " << flush;


      // 8. Define the solution vector x as a finite element grid function
      //    corresponding to fespace. Initialize x with initial guess of zero,
      //    which satisfies the boundary conditions.
      GridFunction x(fespace);
      VectorFunctionCoefficient Boundary_Dirichlet_coef(dim, sol_exact);
      x = 0.0;
      // To use if there are different Dirichlet conditions.
      // Beware, the values of dirichlet boundary conditions are set here !
      //Projète solution exacte sur les bords
      x.ProjectBdrCoefficient(Boundary_Dirichlet_coef, ess_bdr);	

      // 9. Set up the bilinear form a(.,.) on the finite element space
      //    corresponding to the linear elasticity integrator with piece-wise
      //    constants coefficient lambda and mu.

      ConstantCoefficient mu_func(mu);
      ConstantCoefficient lambda_func(lambda);
      BilinearForm *a = new BilinearForm(fespace);
      BilinearFormIntegrator *integ = new ElasticityIntegrator(lambda_func, mu_func);
      a->AddDomainIntegrator(integ);

      // 10. Assemble the bilinear form and the corresponding linear system,
      //     applying any necessary transformations such as: eliminating boundary
      //     conditions, applying conforming constraints for non-conforming AMR,
      //     static condensation, etc.
      cout << "matrix ... " << flush;

      a->Assemble();
      b->Assemble();

      SparseMatrix A;
      Vector B, X;

      a->FormLinearSystem(ess_tdof_list, x, *b, A, X, B);
      cout << "done." << endl;

      cout << "Size of linear system: " << A.Height() << endl;

      if(iterative){
  	GSSmoother M(A);
  	PCG(A, M, B, X, 2, 50000, 1e-20, 0.0);
      } else {
#ifdef MFEM_USE_SUITESPARSE	  
  	// 11. If MFEM was compiled with SuiteSparse, use UMFPACK to solve the system.
  	UMFPackSolver umf_solver;
  	umf_solver.Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_METIS;
  	umf_solver.SetOperator(A);
  	umf_solver.Mult(B, X);
#else
  	cout<<"Direct solver not implemented" << endl;
  	exit(0);
#endif
      }
      // 12. Recover the solution as a finite element grid function.
      a->RecoverFEMSolution(X, *b, x);

      auto end1 = std::chrono::system_clock::now();
      time1 = end1 - start1;
      cout << "Time : " << time1.count()*1000.0 << " ms" << endl;

      // Compute norms of error
      int tdim = dim*(dim+1)/2; // num. entries in a symmetric tensor
      GridFunction zero(fespace);
      zero = 0.0;
      VectorFunctionCoefficient Strain_exacte_coef(tdim,Strain_exacte);
      double ener_reff = ComputeEnergyNorm(zero, lambda_func, mu_func,
					   Strain_exacte_coef);
      double ener_error = ComputeEnergyNorm(x, lambda_func, mu_func,
					    Strain_exacte_coef);
      MatrixFunctionCoefficient grad_exact_coef(dim, Grad_Exact);
      double H1_error = ComputeH1Norm(x,grad_exact_coef);
      VectorFunctionCoefficient sol_exact_coef(dim, sol_exact);
      double L2_error = x.ComputeL2Error(sol_exact_coef);
  	
      cout << "\nL2 norm of error: " << L2_error << endl;
      cout << "Energy norm of error: " << ener_error/ener_reff << endl;
      cout << "H1 norm of error: " << H1_error << endl;

      //Compute norme H1, ener
      VectorFunctionCoefficient Vec_zero_coef(tdim,Vec_zero);
      double ener = ComputeEnergyNorm(x, lambda_func, mu_func,
				      Vec_zero_coef);
      MatrixFunctionCoefficient Mat_zero_coef(dim,Mat_zero);
      double H1 = ComputeH1Norm(x,Mat_zero_coef);

      cout << "Rapport norme energie/H1: "<< ener_error/H1_error <<endl<<endl;

      //double h = mesh->GetElementSize(1);
      double h = fespace->GetNDofs();
      //Compute the slope
      slope_l2(iter,0) = log(err_tmp_l2/L2_error) / log(h_tmp/h);
      slope_l2(iter,1) = h;
      slope_l2(iter,2) = L2_error;
      err_tmp_l2 = L2_error;
      //Compute the slope
      slope_ener(iter,0) = time1.count()*1000.0;
      slope_ener(iter,1) = ener;
      slope_ener(iter,2) = ener_error;
      err_tmp_ener = ener_error;
      //Compute the slope
      slope_grad(iter,0) = log(err_tmp_grad/H1_error) / log(h_tmp/h);
      slope_grad(iter,1) = H1;
      slope_grad(iter,2) = H1_error;
      err_tmp_grad = H1_error;
      h_tmp = h;
      iter++;
      
      //Save in Praview format
      GridFunction diff(fespace);
      GridFunction ex1(fespace);
      diff.ProjectCoefficient(sol_exact_coef);
      ex1.ProjectCoefficient(sol_exact_coef);
      diff-= x;
 
      ParaViewDataCollection paraview_dc("Point_contact", mesh);
      paraview_dc.SetPrefixPath("ParaView");
      paraview_dc.SetLevelsOfDetail(order+1);
      paraview_dc.SetCycle(0);
      paraview_dc.SetDataFormat(VTKFormat::BINARY);
      paraview_dc.SetHighOrderOutput(true);
      paraview_dc.SetTime(0.0); // set the time
      paraview_dc.RegisterField("numerical_solution",&x);
      paraview_dc.RegisterField("diff-exact_solution",&diff);
      paraview_dc.RegisterField("exact_solution",&ex1);
      paraview_dc.Save();	
      

      //Free memory	
      delete a;
      delete b;
      if (fec) {
  	delete fespace;
  	delete fec;
      }
      //    refine the mesh to increase the resolution. In this example we do
      //    'ref_levels' of uniform refinement. We choose 'ref_levels' to be the
      //    largest number that gives a final mesh with no more than 5,000
      //    elements.
      mesh->UniformRefinement();
      cout<<endl;
    }      //end loop mesh

    //Affichage des normes et pentes.
    cout<<endl;
    cout<<"Erreur en norme:"<<endl;
    err_energy_flux << "Erreur en norme:"<<endl;
    for (int i=1; i<iter; i++){
      err_energy_flux << "L2: " << slope_l2(i,2)<<" H1: "<<slope_grad(i,2) 
		      << " Energie: " << slope_ener(i,2)<<" Taille de maille= "
		      <<slope_l2(i,1)<<endl;
      cout << "L2: " << slope_l2(i,2)<<" H1: "<<slope_grad(i,2) 
  	   << " Energie: " << slope_ener(i,2)<<" Taille de maille= "
  	   <<slope_l2(i,1)<<endl;}

    cout<<endl;
    err_energy_flux<<endl;
    cout<<"Pente de convergence:"<<endl;
    err_energy_flux << "Pente de convergence:"<<endl;
    for (int i=1; i<iter; i++){
      err_energy_flux<<"Pente L2: " << slope_l2(i,0)<<" H1: "<<slope_grad(i,0) 
		     << " Energie: " << slope_ener(i,0)<<" Taille de maille= "<<
  	slope_l2(i,1)<<endl;
      cout << "Pente L2: " << slope_l2(i,0)<<" H1: "<<slope_grad(i,0) 
  	   << " Energie: " << slope_ener(i,0)<<" Taille de maille= "<<
  	slope_l2(i,1)<<endl;}
    cout<<endl;

  }	//end flux .txt
  return 0;
}
void Mat_zero(const Vector &x, DenseMatrix &zero)
{
  zero(0,0) = 0.;
  zero(0,1) = 0.;
  zero(1,1) = 0.;
  zero(1,0) = 0.;
}
void Vec_zero(const Vector &x, Vector &zero)
{
  zero(0) = 0.;
  zero(1) = 0.;
  zero(2) = 0.;
}

//===================== Solution exacte =====================
void sol_exact(const Vector &x, Vector &u)
{
  double r2 = x(1)*x(1)+x(0)*x(0);
  double pi = M_PI;
  u(0) = -Force/(4*pi*mu)*(2*x(0)*x(1)/r2 + 2*mu/(lambda+mu)*atan2(x(1),x(0)));
  u(1) = -Force/(4*pi*mu)*((x(1)*x(1)-x(0)*x(0))/r2 - 
			   (lambda+2*mu)*log(r2)/(lambda+mu));
}
//===================== Grad exacte =====================
void Grad_Exact(const Vector &x, DenseMatrix &grad)
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

//===================== Erreur en norme H1 =====================
double ComputeH1Norm(GridFunction &x, MatrixFunctionCoefficient &grad_exact_coef)
{
  FiniteElementSpace *fes = x.FESpace();
  int dim = fes->GetMesh()->SpaceDimension();

  ElementTransformation *Trans;
  DenseMatrix grad, gradh;
  double error = 0.0;
  Array<int> udofs;
  for (int i = 0; i < fes->GetNE() ; i++)
    {
      const FiniteElement *fe = fes->GetFE(i);
      const int order = 2*fe->GetOrder() + 3;   //<----------
      const IntegrationRule *ir = &(IntRules.Get(fe->GetGeomType(), order));
      Trans = fes->GetElementTransformation(i);
      const int dof = fe->GetDof();

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
	  Trans->SetIntPoint(&ip);
	  double w = Trans->Weight() * ip.weight;
	  fe->CalcDShape(ip, dshape);
	  MultAtB(loc_data, dshape, gh);
	  Mult(gh, Trans->InverseJacobian(), gradh);
	  grad_exact_coef.Eval(grad,*Trans,ip);
	  grad -= gradh;
	  error += w * grad.FNorm2();
	}			
    }
  if(error>0.0){
    return sqrt(error);}
  else{
    cout<<"Negative H1 error"<<endl;
    exit(0);}
}

//==============Erreur en Norme energie ===================
double ComputeEnergyNorm(GridFunction &x, Coefficient &lambdah, Coefficient &muh,
			 VectorFunctionCoefficient &Strain_exacte_coef)
{
  FiniteElementSpace *fes = x.FESpace();
  int dim = fes->GetMesh()->SpaceDimension();
  const int tdim = dim*(dim+1)/2; // num. entries in a symmetric tensor
  ElementTransformation *Trans;

  double energy = 0.0;
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
	  ComputeStrain(*Trans, ip, x, i, strainh);
	  Strain_exacte_coef.Eval(strain,*Trans,ip);
	  //======= Stress vectors ========
	  Elasticy_mat(*Trans,ip,dim,lambdah,muh,Ch);
	  Ch.Mult(strainh,stressh);	//approx
	  Ch.Mult(strain,stress);	//exacte

	  strain -= strainh;
	  stress -= stressh;

	  double pdc=0.0;
	  for (int k = 0; k< dim; k++)
	    pdc += strain(k)*stress(k);

	  for (int k = dim; k < dim*(dim+1)/2; k++)
	    pdc += 2*strain(k)*stress(k);

	  energy += w * pdc;
	}
    }
  if(energy>0.0){
    return sqrt(energy);}
  else{
    return 0;}
}

//===================== Matrice élasticité =====================
void Elasticy_mat(ElementTransformation &T,const IntegrationPoint &ip, 
		  int dim, Coefficient &lambda, Coefficient &mu_func, DenseMatrix &C){
  double M = mu_func.Eval(T, ip);
  double L = lambda.Eval(T, ip);

  C.SetSize(dim*(dim+1)/2,dim*(dim+1)/2);
  C = 0.;
  for (int k = 0; k< dim; k++)
    {
      // Extra-diagonal terms
      for (int l = 0; l< dim; l++)
	C(k,l) = L;
      // Diagonal terms
      C(k,k) = L+2*M;
    }
  // Diagonal terms
  for (int k = dim; k < dim*(dim+1)/2; k++)
    C(k,k) = M;
}

//===================== Déformation =====================
void ComputeStrain(ElementTransformation &T,const IntegrationPoint &ip,
		   GridFunction &x, int elem,  Vector &strain){
  FiniteElementSpace *fes = x.FESpace();
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

  double w = T.Weight() * ip.weight;
  fe->CalcDShape(ip, dshape);
  MultAtB(loc_data, dshape, gh);
  Mult(gh, T.InverseJacobian(), grad);
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
    cout<<"dimention not suported"<<endl;}
}


