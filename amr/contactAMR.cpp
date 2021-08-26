// Compile with: make cas_point_contact
//
//Cas test avec solution analytique. Plaque 2D avec point de contact
//d'intensité Force en (0,0).
//Calcul de la déformation en élasticité linéaire. 

//Calcul des erreurs en normes H1 L2 et énergie
//Calcul des pentes de convergences avec ces trois normes.

//Ref: http://e6.ijs.si/medusa/wiki/index.php/Point_contact


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
			 VectorFunctionCoefficient &, Vector &);

void Elasticy_mat(ElementTransformation &,const IntegrationPoint &, int,
		  Coefficient &, Coefficient &, DenseMatrix &);
void ComputeStrain(ElementTransformation &,const IntegrationPoint &,
		   GridFunction &, int,  Vector &);

int main(int argc, char *argv[])
{

  // Parse command-line options.
  const char *mesh_file = "carre.msh";
  int order=1;
  bool static_cond = false;
  bool iterative=true;
  double err_goal = 1;
  int itermax = 100;
  OptionsParser args(argc, argv);
  args.AddOption(&mesh_file, "-m", "--mesh",
		 "Mesh file to use.");
  args.AddOption(&order, "-o", "--order",
		 "Finite element order (polynomial degree).");
  args.AddOption(&static_cond, "-sc", "--static-condensation", "-no-sc",
		 "--no-static-condensation", "Enable static condensation.");
  args.AddOption(&iterative, "-it", "--iterative", "-di",
		 "--direct", "Enable direct or iterative solver.");
  args.AddOption(&err_goal, "-e", "--error",
		 "Local relative Error goal.");
  args.AddOption(&itermax, "-i", "--itmax",
		 "Nombre iteration max de raffinement.");
  args.Parse();
  if (!args.Good())
    {
      args.PrintUsage(cout);
      return 1;
    }
  args.PrintOptions(cout);

  //Lecture du malliage
  std::chrono::duration<double> time1;
  auto start1 = std::chrono::system_clock::now();

  Mesh *mesh = new Mesh(mesh_file, 1, 1);
  mesh->UniformRefinement();
  int dim = mesh->Dimension();
  const int tdim = dim*(dim+1)/2; // num. entries in a symmetric tensor
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

  // 9. Set up the bilinear form a(.,.) on the finite element space
  //    corresponding to the linear elasticity integrator with piece-wise
  //    constants coefficient lambda and mu.

  ConstantCoefficient mu_func(mu);
  ConstantCoefficient lambda_func(lambda);
  BilinearForm *a = new BilinearForm(fespace);
  BilinearFormIntegrator *integ = new ElasticityIntegrator(lambda_func, mu_func);
  a->AddDomainIntegrator(integ);

  //===============Définition estimateur==========================
  FiniteElementSpace flux_fespace(mesh, fec, tdim);
  ZienkiewiczZhuEstimator estimator(*integ, x, flux_fespace);
  estimator.SetAnisotropic();

  ThresholdRefiner refiner(estimator);
  refiner.SetNCLimit(4);

  int it;
  for (it = 1; it<itermax ; it++){		// Boucle raffinement

    // 10. Assemble the bilinear form and the corresponding linear system,
    //     applying any necessary transformations such as: eliminating boundary
    //     conditions, applying conforming constraints for non-conforming AMR,
    //     static condensation, etc.

    // 6. Determine the list of true (i.e. conforming) essential boundary dofs.
    //    In this example, the boundary conditions are defined by marking only
    //    boundary attribute 1 from the mesh as essential and converting it to a
    //    list of true dofs.
    // List of True DoFs : Define (here) Dirichlet conditions
    Array<int> ess_tdof_list;
    Array<int> ess_bdr(mesh->bdr_attributes.Max());
    ess_bdr = 1;		//On sélectionne tous les bords
    fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
    VectorFunctionCoefficient Boundary_Dirichlet_coef(dim, sol_exact);
    x = 0.0;
    // To use if there are different Dirichlet conditions.
    // Beware, the values of dirichlet boundary conditions are set here !
    //Projète solution exacte sur les bords
    x.ProjectBdrCoefficient(Boundary_Dirichlet_coef, ess_bdr);	

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

    //Définition du seuil AMR
    VectorFunctionCoefficient Const_coef1(tdim,Strain_exacte,
					  new ConstantCoefficient(0.0));
    Vector errVec1;
    double ener_refer_ = ComputeEnergyNorm(x, lambda_func, mu_func,
					   Const_coef1,errVec1);
    refiner.SetLocalErrorGoal(ener_refer_*err_goal/(100*sqrt(fespace->GetNE())));
    refiner.SetTotalErrorFraction (0.);

    //===========Raffinement du maillage=================
    refiner.Apply(*mesh);
    fespace->Update();
    x.Update();
    a->Update();
    b->Update();

    if (refiner.Stop())
      {
	cout << "Stopping criterion satisfied. Stop." << endl;
	break;
      }

  }//end reffine loop

  auto end1 = std::chrono::system_clock::now();
  time1 = end1 - start1;
  cout << "Time : " << time1.count()*1000.0 << " ms" << endl;
  cout<<endl;
  // Compute norms of error
  VectorFunctionCoefficient Strain_exacte_coef(tdim,Strain_exacte);
  Vector err_exacte(fespace->GetNE());
  double ener_error = ComputeEnergyNorm(x, lambda_func, mu_func,
					Strain_exacte_coef,err_exacte);
  MatrixFunctionCoefficient grad_exact_coef(dim, Grad_Exact);
  double H1_error = ComputeH1Norm(x,grad_exact_coef);
  VectorFunctionCoefficient sol_exact_coef(dim, sol_exact);
  double L2_error = x.ComputeL2Error(sol_exact_coef);
  GridFunction zero(fespace);
  zero =0.0;
  Vector errVec1;
  double ener_refer = ComputeEnergyNorm(zero, lambda_func, mu_func,
					Strain_exacte_coef,errVec1);
  Vector tmp_estim = estimator.GetLocalErrors();
  cout << "Error Global estimée: "<< estimator.GetTotalError()/ener_refer << endl;
  cout << "Error Global référence: "<< ener_error/ener_refer << endl;
  cout << "Numbre of DOFS: "<< fespace->GetNDofs() << endl;
  cout << endl;

  cout << "Time : " << time1.count()*1000.0 << " ms" << endl;

  cout<<endl<<"Number of reffinement iterations: "<<it<<endl;

  //Save errors gnulplot format
  string const erreur("values_cont.txt");
  ofstream erreur_flux(erreur.c_str());
  if (!erreur_flux.is_open()) {
    cout << "Problem in openning file" << endl;
    exit(0);
  }
  else{
    int NE = mesh->GetNE();
    int count=0;
    for (int i = 0; i < NE; i++)
      {

	erreur_flux<<tmp_estim(i)*sqrt(NE)<<" "<<endl;
	if(tmp_estim(i)*sqrt(NE) > 
	   err_goal/100)
	  {
	    count++;
	  }
	/*
	  if(tmp_estim(i)*sqrt(NE) > 
	  err_goal/100)
	  {
	  count++;
	  erreur_flux<< tmp_estim(i)*sqrt(NE)-err_goal/100<<endl;
	  }
	  else{
	  erreur_flux<< 0.<<endl;
	  }
	*/
      }
    cout<<"Nombre d'élément au dessus du seuil: " << count <<
      " soit "<< (double)count/(double)mesh->GetNE()*100<<"%"<<endl;
  }
  //Save errors gnulplot format
  cout<<"Erreur max: "<<tmp_estim.Max()*sqrt(fespace->GetNE())/*-err_goal/100*/
      <<" Nombre éléments: "<<mesh->GetNE()<<endl;

  string const mesh_gnu("mesh_gnu_cont.txt");
  ofstream mesh_flux(mesh_gnu.c_str());
  if (!mesh_flux.is_open()) {
    cout << "Problem in openning file" << endl;
    exit(0);
  }
  else{
    int NE = mesh->GetNE();
    DenseMatrix coord_(2*NE+1,4);
    for (int i = 0; i < NE; i++)
      {
	Element *el = mesh->GetElement(i);
	int nv = el->GetNVertices();
	int *v = el->GetVertices();
	for (int j = 0; j < nv; j++)
	  {
	    const double *coord = mesh->GetVertex(v[j]);
	    coord_(2*i,j)=coord[0];
	    coord_(2*i+1,j)=coord[1];
	  }
      }

    for (int j = 0; j < 4; j++){
      for (int i = 0; i < NE; i++)
	{
	  //mesh_flux<<coord_(2*i,j)<<" "<<coord_(2*i+1,j)<<" "<<tmp_estim(i)*sqrt(fespace->GetNE())<<" ";
	  mesh_flux<<coord_(2*i,j)<<" "<<coord_(2*i+1,j)<<" ";	
	}
      mesh_flux<<endl;
    }
  }

  //Save in Praview format
  GridFunction ex1(fespace);
  ex1.ProjectCoefficient(sol_exact_coef);
 
  ParaViewDataCollection paraview_dc("ContactAMR", mesh);
  paraview_dc.SetPrefixPath("ParaView");
  paraview_dc.SetLevelsOfDetail(order+1);
  paraview_dc.SetCycle(0);
  paraview_dc.SetDataFormat(VTKFormat::BINARY);
  paraview_dc.SetHighOrderOutput(true);
  paraview_dc.SetTime(0.0); // set the time
  paraview_dc.RegisterField("numerical_solution",&x);
  paraview_dc.RegisterField("exact_solution",&ex1);
  paraview_dc.Save();	

  //Free memory	
  delete a;
  delete b;
  if (fec) {
    delete fespace;
    delete fec;
  }


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
			 VectorFunctionCoefficient &Strain_exacte_coef, Vector &errVec)
{
  FiniteElementSpace *fes = x.FESpace();
  int dim = fes->GetMesh()->SpaceDimension();
  const int tdim = dim*(dim+1)/2; // num. entries in a symmetric tensor
  ElementTransformation *Trans;
  errVec.SetSize(fes->GetNE());
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
      double energy_ = 0.0;
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
	  energy_ += w * pdc;
	}
      errVec(i) = 0.5*sqrt(energy_);
    }
  if(energy>0.0){
    return sqrt(energy);}
  else{
    cout<<"Negative Energy error"<<endl;
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


