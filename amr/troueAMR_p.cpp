//Compile with: make cas_troue
//
//Cas test elsatique linéaire d'une plaque trouée. Un 
//bord est en traction avec la contrainte Sigmainf imposée
//Par symétrie on peut considérer seulement un quart de la plaque
//
//Ref: https://www.researchgate.net/profile/Deghboudj-Samir/
//	publication/306058202_Etude_de_la_concentration_de_con
//	trainte_dans_une_plaque_trouee_sollicitee_en_traction/
//	links/57aceab808ae0932c974cb63/Etude-de-la-concentration-
//	de-contrainte-dans-une-plaque-trouee-sollicitee-en-tractio
//	n.pdf?origin=publication_detail

#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <chrono>
#include <ctime>

using namespace std;
using namespace mfem;

static constexpr double R2 = 1.*1.;	//Rayon troue au carré
static constexpr double Sigmainf = 1.;	//contrainte à l'infinie
static constexpr double E = 1e3;	//Module de Young
static constexpr double nu = 0.3;	//Coef de poisson
static constexpr double lambda = E*nu/((1.+nu)*(1.-2.*nu)); //Coef de lamé
static constexpr double mu = E/(2.*(1.+nu));
static constexpr double c = 10.;	//taille du domaine

//Solution exacte
void conversion(const double, const double , double &, double &);
void Stress_exacteCart(const Vector &x, Vector &stress);
void Stress_exactePol(const Vector &x, Vector &stress);
void Strain_(DenseMatrix &, Vector &);

//Chargement
double NeumannBdr_y(const Vector &);
double NeumannBdr_x(const Vector &);

//Erreur en norme energy
double ComputeEnergyNorm(ParGridFunction &,
			 Coefficient &, Coefficient &,
			 VectorFunctionCoefficient &, Vector &);
//Calcul erreur absolue
double ComputeAbsoluErr(ParGridFunction &, 
			double &,double &);
double ComputeRelativErr(ParGridFunction &, 
			 double &,double &);
Vector SaveError(const Vector &);
//Matrice d'élaticité
void Elasticy_mat(ElementTransformation &,const IntegrationPoint &,
		  double &, double &, DenseMatrix &);
void ComputeStrain(ElementTransformation &,const IntegrationPoint &,
		   ParGridFunction &,  Vector &);
void ComputeStress(ElementTransformation &,const IntegrationPoint &,
		   ParGridFunction &,Coefficient &, Coefficient &, Vector &, Vector &);

//Save stress QuadratureFunction format
void SaveStress(ParGridFunction &,Coefficient &, Coefficient &,
		Mesh &, QuadratureFunction );
//Changement de base
void ChangeBase(const double &,const double &, Vector &, Vector &);
//Recherche taille de maille max
double StepMax(Mesh &);
class DiagCoefficient : public Coefficient
{
protected:
  Coefficient &lambda, &mu;
  ParGridFunction *u; // displacement
  int si, sj; // component to evaluate, 0 <= si,sj < dim

  DenseMatrix grad; // auxiliary matrix, used in Eval

public:
  DiagCoefficient(Coefficient &lambda_, Coefficient &mu_)
    : lambda(lambda_), mu(mu_), u(NULL), si(0), sj(0) { }

  void SetDisplacement(ParGridFunction &u_) { u = &u_; }

  virtual double Eval(ElementTransformation &T, const IntegrationPoint &ip) = 0;
};

class StressCoefficient : public DiagCoefficient
{
public:
  using DiagCoefficient::DiagCoefficient;
  void SetComponent(int i, int j) { si = i; sj = j; }
  virtual double Eval(ElementTransformation &T, const IntegrationPoint &ip);

};

#define USE_PROFILER 1
#define LIB_PROFILER_IMPLEMENTATION
#define LIB_PROFILER_PRINTF MpiPrintf
#include <mpi.h>
#include "libProfiler.h"

int main(int argc, char *argv[])
{
  PROFILER_ENABLE;
  // initialize mpi.
  int num_procs, myid;
  MPI_Init(&argc, &argv);
  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  PROFILER_START(0_total);
  PROFILER_START(1_initialize_mesh);

  // Parse command-line options.
  bool static_cond = false;
  int order = 1;
  bool amg_elast = 0;
  bool iterative = true;
  bool reorder_space = false;
  int initial_mesh_refinement = 1;
  double err_goal = 1.;
  int itermax = 100;
  const char *mesh_file1 = "hole_mesh/quarter_phole1.msh";
  OptionsParser args(argc, argv);
  args.AddOption(&mesh_file1, "-m", "--mesh",
		 "Mesh file to use.");
  args.AddOption(&order, "-o", "--order",
		 "Finite element order (polynomial degree).");
  args.AddOption(&amg_elast, "-elast", "--amg-for-Elasticity", "-sys",
 		 "--amg-for-systems",
 		 "use the special amg Elasticity solver (gm/ln approaches), "
 		 "or standard amg for systems (unknown approach).");
  args.AddOption(&reorder_space, "-nodes", "--by-nodes", "-vdim", "--by-vdim",
 		 "use bynodes ordering of Vector space instead of byvdim");
  args.AddOption(&initial_mesh_refinement, "-r", "--initial_refinement",
		 "Number of inital uniform refinement.");
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

  Mesh *mesh = new Mesh(mesh_file1, 1, 1);
  for (int i = 0 ; i < initial_mesh_refinement ; i++)
    {
      // Perform refinement
      mesh->UniformRefinement();
    }
  mesh->EnsureNCMesh(true);

  ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);

  // Read the mesh from the given mesh file. We can handle triangular,
  //    quadrilateral, tetrahedral or hexahedral elements with the same code.
  // Ensure non-conforming mesh
  int dim = mesh->Dimension();
  const int tdim = dim*(dim+1)/2; // num. entries in a symmetric tensor
  PROFILER_END(); PROFILER_START(2_initialize_fem);
  // Define a finite element space on the mesh-> Here we use vector finite
  //    elements, i.e. dim copies of a scalar finite element space. The vector
  //    dimension is specified by the last argument of the ParFiniteElementSpace
  //    constructor. For NURBS meshes, we use the (degree elevated) NURBS space
  //    associated with the mesh nodes.
  FiniteElementCollection *fec;
  ParFiniteElementSpace *fespace;

  // 6. Define a finite element space on the mesh. The polynomial order is one
  //    (linear) by default, but this can be changed on the command line.
  fec = new H1_FECollection(order, dim);
  fespace =  new ParFiniteElementSpace(pmesh, fec, dim,  Ordering::byVDIM);;

  HYPRE_Int Size = fespace->GlobalTrueVSize();
  if (myid == 0)
    {
      cout << "Numbers of elements: " << pmesh->GetNE() <<endl;
	  cout << "assembling: " << flush;
    }


  //  Define the solution vector x as a finite element grid function
  //    corresponding to fespace. Initialize x with initial guess of zero,
  //    which satisfies the boundary conditions.
  ParGridFunction x(fespace);
  x = 0.0;

  //  Set up the linear form b(.) which corresponds to the right-hand side of
  //    the FEM linear system. In this case, b_i equals the boundary integral
  //    of f*phi_i where f represents a "pull down" force on the Neumann part
  //    of the boundary and phi_i are the basis functions in the finite element
  //    fespace. The force is defined by the VectorArrayCoefficient object f,
  //    which is a vector of Coefficient objects.
  //    The fact that f is non-zero
  //    on boundary attribute 2 is indicated by the use of piece-wise constants
  //    coefficient for its last component.
  // Add non zero Neumann conditions
  //Boundary conditions markers
  Array<int> bdr_attr_marker_neu(pmesh->bdr_attributes.Max());
  // values for neumann boundary conditions are set within boundary function
  bdr_attr_marker_neu = 0;
  bdr_attr_marker_neu[4] = 1;	//gauche
  bdr_attr_marker_neu[0] = 1;	//bas
  VectorArrayCoefficient fo(dim);
  fo.Set(1, new FunctionCoefficient(NeumannBdr_y));
  fo.Set(0, new FunctionCoefficient(NeumannBdr_x));
  ParLinearForm *b = new ParLinearForm(fespace);
  b->AddBoundaryIntegrator(new VectorBoundaryLFIntegrator(fo),bdr_attr_marker_neu);
  if (myid == 0)
    {
      cout << "r.h.s. ... " << endl;
    }
  
  //  Set up the bilinear form a(.,.) on the finite element space
  //    corresponding to the linear elasticity integrator with piece-wise
  //    constants coefficient lambda and mu.
  
  ConstantCoefficient mu_func(mu);
  ConstantCoefficient lambda_func(lambda);
  ParBilinearForm *a = new ParBilinearForm(fespace);
  BilinearFormIntegrator *integ = new ElasticityIntegrator(lambda_func, mu_func);
  a->AddDomainIntegrator(integ);

  ParGridFunction zero(fespace);
  zero =0.0;
  Vector errVec1;
  VectorFunctionCoefficient Stress_exactepol_coef(tdim,Stress_exactePol);
  double ener_refer = ComputeEnergyNorm(zero, lambda_func, mu_func,
					Stress_exactepol_coef,errVec1);

  //===============Définition estimateur==========================
  PROFILER_END(); PROFILER_START(3_amr_loop);
  L2_FECollection flux_fec(order, dim);
  RT_FECollection smooth_flux_fec(order-1, dim);
  ErrorEstimator* estimator;
  /*
    auto flux_fes = new ParFiniteElementSpace(pmesh, 
    &flux_fec, tdim);
    auto smooth_flux_fes = new ParFiniteElementSpace(pmesh,
    &smooth_flux_fec);
    estimator = new L2ZienkiewiczZhuEstimator(*integ, x, flux_fes, smooth_flux_fes);
  */
  auto flux_fes = new ParFiniteElementSpace(pmesh, fec, tdim);
  estimator = new ZienkiewiczZhuEstimator(*integ, x, flux_fes);

  ThresholdRefiner refiner(*estimator);
  refiner.SetNCLimit(4);

  int it;
  for (it = 0; it<itermax; it++){		// Boucle raffinement
    PROFILER_START(3_1_construct_system);
    HYPRE_Int global_dofs = fespace->GlobalTrueVSize();
    if (myid == 0)
      {
	cout << "Iteration: " << it << ", number of unknowns: "
	     << global_dofs << endl;
      }

    //  Determine the list of true (i.e. conforming) essential boundary dofs.
    //    In this example, the boundary conditions are defined by marking only
    //    boundary attribute 1 from the mesh as essential and converting it to a
    //    list of true dofs.
    // List of True DoFs : Define (here) Dirichlet conditions

    if (myid == 0) { cout << "matrix ... " << flush; }
    Array<int> ess_tdof_list, tmp_tdof, ess_bdr(pmesh->bdr_attributes.Max());
    ess_bdr = 0;
    ess_bdr[3] = 1;	//droite
    fespace->GetEssentialTrueDofs(ess_bdr, tmp_tdof, 1); 
    // ess_tof_list accumulates all needed dof
    ess_tdof_list.Append(tmp_tdof);
    ess_bdr = 0;
    ess_bdr[1] = 1;	//haut
    fespace->GetEssentialTrueDofs(ess_bdr, tmp_tdof, 0);
    // ess_tof_list accumulates all needed dof
    ess_tdof_list.Append(tmp_tdof);

    x = 0.0;
    //  Assemble the bilinear form and the corresponding linear system,
    //     applying any necessary transformations such as: eliminating boundary
    //     conditions, applying conforming constraints for non-conforming AMR,
    //     static condensation, etc.
    // 15. Recompute the field on the current mesh: assemble the stiffness
    //     matrix and the right-hand side.
    a->Assemble();
    b->Assemble();

    HypreParMatrix A;
    Vector B, X;
    a->FormLinearSystem(ess_tdof_list, x, *b, A, X, B);

    PROFILER_END(); PROFILER_START(3_2_solve_system);
    if (myid == 0) {cout<<"solve ... "<<flush; }
      if(iterative){
    HypreBoomerAMG amg(A);
    amg.SetPrintLevel(0);
    HyprePCG pcg(A);
    pcg.SetTol(1e-12);
    pcg.SetMaxIter(5000);
    pcg.SetPrintLevel(0);
    pcg.SetPreconditioner(amg);
    pcg.Mult(B, X);
      } else {
#ifdef MFEM_USE_SUITESPARSE	  
	// If MFEM was compiled with SuiteSparse, use UMFPACK to solve the system.
	UMFPackSolver umf_solver;
	umf_solver.Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_METIS;
	umf_solver.SetOperator(A);
	umf_solver.Mult(B, X);
#else
	cout<<"Direct solver not implemented" << endl;
	exit(0);
#endif
      }
    // 18. Extract the local solution on each processor.
    a->RecoverFEMSolution(X, *b, x);
    PROFILER_END(); PROFILER_START(3_3_amr_refine);
    //Définition du seuil AMR
    VectorFunctionCoefficient Const_coef1(tdim,Stress_exacteCart,
					  new ConstantCoefficient(0.0));
    Vector errVec1;
    double ener_refer_ = ComputeEnergyNorm(x, lambda_func, mu_func,
					   Const_coef1,errVec1);
    const int NE = pmesh->GetNE();
    int NE_tot = NE;
    MPI_Allreduce(&NE, &NE_tot, 1, MPI_INT,
		  MPI_SUM, MPI_COMM_WORLD);

    refiner.SetLocalErrorGoal(ener_refer_*err_goal/(100*sqrt(NE_tot)));
    refiner.SetTotalErrorFraction (0.);

    //===========Raffinement du maillage=================
    if (myid == 0) {cout<<"Refine ... "<<flush; }
    refiner.Apply(*pmesh);

    if (myid == 0) {cout<<"Update ... "<<endl<<endl; }
    PROFILER_END(); PROFILER_START(3_4_update);
    fespace->Update();
    x.Update();

    if (pmesh->Nonconforming())
      {

	// Load balance the mesh.
	pmesh->Rebalance();

	// Update the space again, this time a GridFunction redistribution matrix
	// is created. Apply it to the solution.
	fespace->Update();
	x.Update();
      }

    a->Update();
    b->Update();
    // Free any transformation matrices to save memory.
    fespace->UpdatesFinished();
    PROFILER_END();
    if (refiner.Stop())
      {
	if (myid == 0) {cout << "Stopping criterion satisfied. Stop." << endl;}
	break;
      }
  }//end reffine loop
  PROFILER_END();

  PROFILER_END();
  if (myid == 0)
    {
      cout << endl;
      LogProfiler();
    }
  PROFILER_DISABLE;
  Vector err_exacte(fespace->GetNE());
  double ener_error = ComputeEnergyNorm(x, lambda_func, mu_func,
					Stress_exactepol_coef,err_exacte);
  Vector tmp_estim = estimator->GetLocalErrors();
  HYPRE_Int global_dofs_ = fespace->GlobalTrueVSize();
  if (myid == 0)
    {
      cout<<endl;
      cout << "Error Global estimée: "<< estimator->GetTotalError()
	/ener_refer*100 << endl;
      cout << "Error Global référence: "<< ener_error/ener_refer*100 << endl;
      cout << "Number of DOFS: "<< global_dofs_ << endl;
	  cout<<endl<<"Number of reffinement iterations: "<<it<<endl;
    }

  //Save in Praview format
  if (!mesh->NURBSext)
    {
      pmesh->SetNodalFESpace(fespace);
    }
  ParaViewDataCollection paraview_dc("TroueAMR_p", pmesh);
  paraview_dc.SetPrefixPath("ParaView");

  //Stress exacte
  ParFiniteElementSpace *fieldspace_e;
  fieldspace_e = new ParFiniteElementSpace(pmesh, fec, tdim);
  ParGridFunction stress_exactepol(fieldspace_e);
  stress_exactepol.ProjectCoefficient(Stress_exactepol_coef);

  //  paraview_dc.SetLevelsOfDetail(1);
  paraview_dc.SetCycle(0);
  paraview_dc.SetDataFormat(VTKFormat::BINARY);
  //  paraview_dc.SetHighOrderOutput(true);
  paraview_dc.SetTime(0.0); // set the time
  paraview_dc.RegisterField("numerical_solution",&x);
  paraview_dc.RegisterField("Stress exacte polaire",&stress_exactepol);

  paraview_dc.Save();
  delete a;
  delete b;
  if (fec)
    {
      delete fespace;
      delete fec;
    }
  delete mesh;
  delete pmesh;

  MPI_Finalize();
  return 0;
}
void conversion(const double x, const double y , double &r, double &theta)
{
  r = sqrt(x*x + y*y);
  theta = atan2(y,x);
  if(x==0){
    theta = 3*M_PI/2;}
}
//===================== Stress exacte =====================
void Stress_exactePol(const Vector &x, Vector &stress_cart)
{
  double r, theta;
  conversion(x(0), x(1), r, theta);
  int dim = x.Size();
  double r2 = r*r;
  Vector stress(dim*(dim+1)/2);
  stress(0) = Sigmainf*0.5*((1.-R2/r2)+(1.+3*R2*R2/(r2*r2)-
					4*R2/r2)*cos(2*theta));
  stress(1) = Sigmainf*0.5*((1.+R2/r2) - (1.+3*R2*R2/(r2*r2))*
			    cos(2*theta));
  stress(2) = -Sigmainf*0.5*(1. - 3*R2*R2/(r2*r2) + 2*R2/r2)*
    sin(2*theta);
  ChangeBase(x(0),x(1),stress,stress_cart);
}
void Stress_exacteCart(const Vector &x, Vector &stress)
{
  double r, theta;
  conversion(x(0), x(1), r, theta);
  int dim = x.Size();
  double r2 = r*r;
  stress(0) = Sigmainf*(1-R2/r2*(3/2*cos(2*theta)+cos(4*theta))
			+3/2*R2*R2/(r2*r2)*cos(4*theta));
  stress(1) = Sigmainf*(-R2/r2*(0.5*cos(2*theta)-cos(4*theta))
			-3/2*R2*R2/(r2*r2)*cos(4*theta));
  stress(2) = Sigmainf*(-R2/r2*(0.5*sin(2*theta)+sin(4*theta))
			+3/2*R2*R2/(r2*r2)*sin(4*theta));
}
//===================== Condition Neumann exacte =====================
double NeumannBdr_x(const Vector &x)
{
  int dim = x.Size(), tdim = dim*(dim+1)/2;
  Vector stress(tdim);
  Stress_exactePol(x,stress);
  if(x(0) == -c){
    return -stress(0);//composante suivant x sur le bord gauche
  }
  if(x(1) == -c){
    return -stress(2);//composante suivant x sur le bord bas
  }
  else{
    return 0.;
  }
}
double NeumannBdr_y(const Vector &x)
{
  int dim = x.Size(), tdim = dim*(dim+1)/2;
  Vector stress(tdim);
  Stress_exactePol(x,stress);
  if(x(0) == -c){
    return -stress(2);//composante suivant y sur le bord gauche
  }
  if(x(1) == -c){
    return -stress(1);//composante suivant y sur le bord bas
  }
  else{
    return 0.;
  }
}

double StressCoefficient::Eval(ElementTransformation &T,
                               const IntegrationPoint &ip)
{
  MFEM_ASSERT(u != NULL, "displacement field is not set");
  double L = lambda.Eval(T, ip);
  double M = mu.Eval(T, ip);
  u->GetVectorGradient(T, grad);
  if (si == sj)
    {
      double div_u = grad.Trace();
      return (L*div_u + 2*M*grad(si,si));
    }
  else
    {
      return M*(grad(si,sj) + grad(sj,si));
    }
}
//===================== Strain exacte =====================
void Strain_(DenseMatrix &grad, Vector &strain)
{
  int dim = grad.Size();
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
//===================== Calcul erreur absolue =====================
double ComputeRelativErr(ParGridFunction &x, double &ener_error, double &absolu_err)
{
  return (absolu_err)/ener_error;
}

double ComputeAbsoluErr(ParGridFunction &x, double &ener_error,double &relative_err)
{
  return relative_err*ener_error/100;
}

//==============Erreur en Norme energie ===================
double ComputeEnergyNorm(ParGridFunction &x, Coefficient &lambdah, Coefficient &muh,
			 VectorFunctionCoefficient &Stress_exacte_coef, Vector &errVec)
{
  ParFiniteElementSpace *fes = x.ParFESpace();
  int dim = fes->GetMesh()->SpaceDimension();
  const int tdim = dim*(dim+1)/2; // num. entries in a symmetric tensor
  ElementTransformation *Trans;
  errVec.SetSize(fes->GetNE());
  double error_local = 0.0, error_global = 0.0;
  for (int i = 0; i < fes->GetNE() ; i++)
    {
      const FiniteElement *fe = fes->GetFE(i);
      const int order = 2*fe->GetOrder()+3; // <----------
      const IntegrationRule *ir = &IntRules.Get(fe->GetGeomType(), order);
      const int dof = fe->GetDof();
      Trans = fes->GetElementTransformation(i);
      Vector stressh(tdim), strainh(tdim), strainh_pol(tdim);	//approché
      Vector stress(tdim), strain(tdim);	//exacte
      DenseMatrix C;
      double energy_ = 0.0;
      for (int j = 0; j < ir->GetNPoints(); j++)
	{
	  const IntegrationPoint &ip = ir->IntPoint(j);
	  Trans->SetIntPoint(&ip);
	  double w = Trans->Weight() * ip.weight;
	  ComputeStress(*Trans, ip, x, lambdah, muh, stressh, strainh);	
	  Stress_exacte_coef.Eval(stress,*Trans,ip);//exact
	  double M = muh.Eval(*Trans, ip);
	  double L = lambdah.Eval(*Trans, ip);
	  Elasticy_mat(*Trans,ip,L,M,C);
	  C.Invert();
	  C.Mult(stress,strain);

	  strainh -= strain;
	  stressh -= stress;

	  double pdc=0.0;
	  for (int k = 0; k< dim; k++)
	    pdc += stressh(k)*strainh(k);
	  for (int k = dim; k < dim*(dim+1)/2; k++)
	    pdc += 2*stressh(k)*strainh(k); 

	  error_local += w * pdc;
	  energy_ += w * pdc;
	}
    }

  MPI_Allreduce(&error_local, &error_global, 1, MPI_DOUBLE, MPI_SUM,
		MPI_COMM_WORLD);
  if(error_global>0.0){
    return sqrt(error_global);}
  else{
    cout<<"Negative energy Norme "<<endl;
    MPI_Finalize();
    exit(0);
  }
}


//===================== Matrice élasticité =====================
void Elasticy_mat(ElementTransformation &T,const IntegrationPoint &ip, 
		  double &L, double &M, DenseMatrix &C){
  int dim = T.GetSpaceDim();
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
    C(k,k) = 2*M;
}

//===================== Déformation =====================
void ComputeStrain(ElementTransformation &T,const IntegrationPoint &ip,
		   ParGridFunction &x, Vector &strain){
  ParFiniteElementSpace *fes = x.ParFESpace();
  Array<int> udofs;
  const FiniteElement *fe = fes->GetFE(T.ElementNo);
  const int dof = fe->GetDof();
  const int dim = fe->GetDim();
  const int tdim = dim*(dim+1)/2; // num. entries in a symmetric tensor
  DenseMatrix dshape(dof, dim);
  DenseMatrix gh(dim, dim),grad (dim, dim);
  fes->GetElementVDofs(T.ElementNo, udofs);
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
  fe->CalcDShape(ip, dshape);
  MultAtB(loc_data, dshape, gh);
  Mult(gh, T.InverseJacobian(), grad);
  Strain_(grad, strain);
}
//===================== Contrainte =====================
void ComputeStress(ElementTransformation &T,const IntegrationPoint &ip,
		   ParGridFunction &x,Coefficient &lambda, Coefficient &mu, 
		   Vector &stress_cart, Vector &strainh_cart){
  double L = lambda.Eval(T, ip);
  double M = mu.Eval(T, ip);
  int dim = T.GetSpaceDim();
  int tdim = dim*(dim+1)/2;
  Vector  stress(tdim), strainh(tdim);
  stress.SetSize(tdim);
  ComputeStrain(T, ip, x, strainh_cart);	//calcul de la déformation
  DenseMatrix C;
  Elasticy_mat(T,ip,L,M,C);	//matrice d'élasticité
  C.Mult(strainh_cart,stress_cart); //calcul de la contrainte
}
//===================== Changement de base =====================
void ChangeBase(const double &x,const double &y, 
		Vector &stress_cart, Vector &Stress_pol){
  int dim = 2;
  double r, theta;
  conversion(x,y,r,theta);
  Stress_pol(0)=cos(theta)*cos(theta)*stress_cart(0)+sin(theta)*sin(theta)*stress_cart(1)
    -sin(2*theta)*stress_cart(2);  Stress_pol(1)=cos(theta)*cos(theta)*stress_cart(1)+sin(theta)*sin(theta)*stress_cart(0)
				     +sin(2*theta)*stress_cart(2);
  Stress_pol(2)=cos(theta)*sin(theta)*(stress_cart(0)-stress_cart(1))
    +cos(2*theta)*stress_cart(2);
}
