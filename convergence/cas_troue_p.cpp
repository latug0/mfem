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
//
//Version Parallèle
#include "mfem.hpp"
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <time.h>

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
			 VectorFunctionCoefficient &);
//Matrice d'élaticité
void Elasticy_mat(ElementTransformation &,const IntegrationPoint &,
		  double &, double &, DenseMatrix &);
void ComputeStrain(ElementTransformation &,const IntegrationPoint &,
		   ParGridFunction &,  Vector &);
void ComputeStress(ElementTransformation &,const IntegrationPoint &,
		   ParGridFunction &,Coefficient &, Coefficient &, Vector &, Vector &);

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
  int time = 0;
  PROFILER_ENABLE;
  // initialize mpi.
  int num_procs, myid;
  MPI_Init(&argc, &argv);
  MPI_Barrier(MPI_COMM_WORLD);
  PROFILER_START(Total);
  PROFILER_START(0_initialize);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  // Parse command-line options.
  const char *mesh_file = "hole_mesh/quarter_phole00.msh";
  bool static_cond = false;
  int order = 1;
  bool amg_elast = 0;
  bool reorder_space = false;
  bool iterative = true;
  int ref_levels = 0;
  OptionsParser args(argc, argv);
  args.AddOption(&mesh_file, "-m", "--mesh",
 		 "mesh file to use.");
  args.AddOption(&order, "-o", "--order",
		 "Finite element order (polynomial degree).");
  args.AddOption(&amg_elast, "-elast", "--amg-for-Elasticity", "-sys",
 		 "--amg-for-systems",
 		 "use the special amg Elasticity solver (gm/ln approaches), "
 		 "or standard amg for systems (unknown approach).");
  args.AddOption(&static_cond, "-sc", "--static-condensation", "-no-sc",
		 "--no-static-condensation", "Enable static condensation.");
  args.AddOption(&iterative, "-it", "--iterative", "-di",
		 "--direct", "Enable direct or iterative solver.");
  args.AddOption(&reorder_space, "-nodes", "--by-nodes", "-vdim", "--by-vdim",
 		 "use bynodes ordering of Vector space instead of byvdim");
  args.AddOption(&ref_levels, "-r", "--num_ref", "nombre de rafinement de maillage");
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
   PROFILER_END(); PROFILER_START(1_read_mesh);
  //lecture du malliage
  Mesh *mesh = new Mesh(mesh_file, 1, 1);

  for (int r = 0; r < ref_levels; r++)
    {
      mesh->UniformRefinement();
    }

  ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
  // Read the mesh from the given mesh file. We can handle triangular,
  //    quadrilateral, tetrahedral or hexahedral elements with the same code.
  int dim = mesh->Dimension();
  PROFILER_END(); PROFILER_START(2_initialize_fem);
  // Define a finite element space on the mesh. Here we use vector finite
  //    elements, i.e. dim copies of a scalar finite element space. The vector
  //    dimension is specified by the last argument of the FiniteElementSpace
  //    constructor. For NURBS meshes, we use the (degree elevated) NURBS space
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
  
  //  Determine the list of true (i.e. conforming) essential boundary dofs.
  //    In this example, the boundary conditions are defined by marking only
  //    boundary attribute 1 from the mesh as essential and converting it to a
  //    list of true dofs.
  
  // List of True DoFs : Define (here) Dirichlet conditions
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
      cout << "r.h.s. ... " << flush;
    }

  //  Set up the bilinear form a(.,.) on the finite element space
  //    corresponding to the linear elasticity integrator with piece-wise
  //    constants coefficient lambda and mu.
  
  ConstantCoefficient mu_func(mu);
  ConstantCoefficient lambda_func(lambda);
  ParBilinearForm *a = new ParBilinearForm(fespace);
  BilinearFormIntegrator *integ = new ElasticityIntegrator(lambda_func, mu_func);
  a->AddDomainIntegrator(integ);
  PROFILER_END(); PROFILER_START(3_assembly);
  //  Assemble the bilinear form and the corresponding linear system,
  //     applying any necessary transformations such as: eliminating boundary
  //     conditions, applying conforming constraints for non-conforming AMR,
  //     static condensation, etc.
  if (myid == 0) { cout << "matrix ... " << flush; }
  if (static_cond) { a->EnableStaticCondensation(); }
  a->Assemble();
  b->Assemble();

  PROFILER_END(); PROFILER_START(4_solve);

  HypreParMatrix A;
  Vector B, X;
  a->FormLinearSystem(ess_tdof_list, x, *b, A, X, B);  
  if (myid == 0)
    {
      cout << "done." << endl;
      cout << "Size of linear system: " << A.GetGlobalNumRows() << endl;
    }
  
  if(iterative){
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
  }
  else{
    cout<<"solver direct non implémenté"<<endl;
  }   
  // Recover the solution as a finite element grid function.
  a->RecoverFEMSolution(X, *b, x);
  if (!use_nodal_fespace)
    {
      pmesh->SetNodalFESpace(fespace);
    }
   PROFILER_END();
   PROFILER_END();
  if (myid == 0)
    {
   LogProfiler();
    }
   PROFILER_DISABLE;

  // Compute norms of error
  int tdim = dim*(dim+1)/2; // num. entries in a symmetric tensor
  VectorFunctionCoefficient Stress_exactecart_coef(tdim,Stress_exacteCart);
  VectorFunctionCoefficient Stress_exactepol_coef(tdim,Stress_exactePol);
  double ener_error = ComputeEnergyNorm(x, lambda_func, mu_func,
					Stress_exactepol_coef);
  double h = StepMax(*mesh);	//recherche pas max
  if (myid == 0)
    {
      cout << "Energy norm of error: " << ener_error <<" Taille de maille= "
	   <<h<< endl<<endl;
    }


  //  free the used memory.
  delete a;
  delete b;
  if (fec)
    {
      delete fespace;
      delete fec;
    }
  delete mesh;
  delete pmesh;

  time = clock();
  if (myid == 0)
    {
  printf("Temps d'execution = %d 10^-5s \n", time);
	}
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

//==============Erreur en Norme energie ===================
double ComputeEnergyNorm(ParGridFunction &x, Coefficient &lambdah, 
			 Coefficient &muh,VectorFunctionCoefficient &Stress_exacte_coef)
{
  ParFiniteElementSpace *fes = x.ParFESpace();
  int dim = fes->GetMesh()->SpaceDimension();
  const int tdim = dim*(dim+1)/2; // num. entries in a symmetric tensor
  ElementTransformation *Trans;
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
	}
    }
  MPI_Reduce(&error_local, &error_global, 1, MPI_DOUBLE, MPI_SUM, 0,
 	     MPI_COMM_WORLD);
  if(error_global>0.0){
    return sqrt(error_global);}
  else{
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
//===================== Recherche pas max =====================
double StepMax(Mesh &mesh){
  double h = 0., h_tmp;
  for (int i = 0; i < mesh.GetNE() ; i++)
    {
      h_tmp = mesh.GetElementSize(i);
      if(h_tmp > h){
	h = h_tmp;
      }
    }
  return h;
}
