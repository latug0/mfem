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
double ComputeEnergyNorm(GridFunction &,
			 Coefficient &, Coefficient &,
			 VectorFunctionCoefficient &);
//Calcul erreur absolue
double ComputeAbsoluErr(GridFunction &, 
			Coefficient &, Coefficient &, 
			double &);
//Matrice d'élaticité
void Elasticy_mat(ElementTransformation &,const IntegrationPoint &,
		  double &, double &, DenseMatrix &);
void ComputeStrain(ElementTransformation &,const IntegrationPoint &,
		   GridFunction &,  Vector &);
void ComputeStress(ElementTransformation &,const IntegrationPoint &,
		   GridFunction &,Coefficient &, Coefficient &, Vector &, Vector &);

//Save stress QuadratureFunction format
void SaveStress(GridFunction &,Coefficient &, Coefficient &,
		Mesh &, QuadratureFunction );
//Changement de base
void ChangeBase(const double &,const double &, Vector &, Vector &);
//Recherche taille de maille max
double StepMax(Mesh &);
class DiagCoefficient : public Coefficient
{
protected:
  Coefficient &lambda, &mu;
  GridFunction *u; // displacement
  int si, sj; // component to evaluate, 0 <= si,sj < dim

  DenseMatrix grad; // auxiliary matrix, used in Eval

public:
  DiagCoefficient(Coefficient &lambda_, Coefficient &mu_)
    : lambda(lambda_), mu(mu_), u(NULL), si(0), sj(0) { }

  void SetDisplacement(GridFunction &u_) { u = &u_; }

  virtual double Eval(ElementTransformation &T, const IntegrationPoint &ip) = 0;
};

class StressCoefficient : public DiagCoefficient
{
public:
  using DiagCoefficient::DiagCoefficient;
  void SetComponent(int i, int j) { si = i; sj = j; }
  virtual double Eval(ElementTransformation &T, const IntegrationPoint &ip);

};

int main(int argc, char *argv[])
{

  DenseMatrix slope_ener;
  double err_tmp_ener = 0;
  double h_tmp = 0.;
  int iter = 0;
  // Parse command-line options.
  bool static_cond = false;
  int order = 1;
  bool iterative = true;
  int rep = 1;
  const char *mesh_file1 = "hole_mesh/quarter_phole1.msh";
  OptionsParser args(argc, argv);
  args.AddOption(&mesh_file1, "-m", "--mesh",
		 "Mesh file to use.");
  args.AddOption(&order, "-o", "--order",
		 "Finite element order (polynomial degree).");
  args.AddOption(&rep, "-r", "--repet",
		 "Repetition de malliage.");
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
  slope_ener.SetSize(rep,3);
  Mesh *mesh = new Mesh(mesh_file1, 1, 1);
  for (int r = 0; r < rep; r++){
    /*
      string  mesh_file = "hole_mesh/quarter_phole";
      string buf(mesh_file);
      string X_(to_string(r));
      buf.append(X_);
      string buff(buf);
      string mesh_file2 = ".msh";
      buff.append(mesh_file2);

      char* mesh_file_ = &buff[0];
      cout<<"Nom du maillage: "<<mesh_file_<<endl;
      Mesh *mesh = new Mesh(mesh_file_, 1, 1);
    */
    // Read the mesh from the given mesh file. We can handle triangular,
    //    quadrilateral, tetrahedral or hexahedral elements with the same code.
    int dim = mesh->Dimension();
  
    // Define a finite element space on the mesh. Here we use vector finite
    //    elements, i.e. dim copies of a scalar finite element space. The vector
    //    dimension is specified by the last argument of the FiniteElementSpace
    //    constructor. For NURBS meshes, we use the (degree elevated) NURBS space
    //    associated with the mesh nodes.
    FiniteElementCollection *fec;
    FiniteElementSpace *fespace;
    fec = new H1_FECollection(order, dim);
    fespace = new FiniteElementSpace(mesh, fec, dim, Ordering::byVDIM);
    cout << "Numbers of elements: " << mesh->GetNE() <<endl;
    cout << "Number of finite element unknowns: " << fespace->GetTrueVSize()
	 << endl << "Assembling: " << flush;
  
    //  Determine the list of true (i.e. conforming) essential boundary dofs.
    //    In this example, the boundary conditions are defined by marking only
    //    boundary attribute 1 from the mesh as essential and converting it to a
    //    list of true dofs.
  
    // List of True DoFs : Define (here) Dirichlet conditions
    Array<int> ess_tdof_list, tmp_tdof, ess_bdr(mesh->bdr_attributes.Max());
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
    GridFunction x(fespace);
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
    Array<int> bdr_attr_marker_neu(mesh->bdr_attributes.Max());
    // values for neumann boundary conditions are set within boundary function
    bdr_attr_marker_neu = 0;
    bdr_attr_marker_neu[4] = 1;	//gauche
    bdr_attr_marker_neu[0] = 1;	//bas
    VectorArrayCoefficient fo(dim);
    fo.Set(1, new FunctionCoefficient(NeumannBdr_y));
    fo.Set(0, new FunctionCoefficient(NeumannBdr_x));
    LinearForm *b = new LinearForm(fespace);
    b->AddBoundaryIntegrator(new
			     VectorBoundaryLFIntegrator(fo),bdr_attr_marker_neu);
    cout << "r.h.s. ... " << flush;
  
    //  Set up the bilinear form a(.,.) on the finite element space
    //    corresponding to the linear elasticity integrator with piece-wise
    //    constants coefficient lambda and mu.
  
    ConstantCoefficient mu_func(mu);
    ConstantCoefficient lambda_func(lambda);
    BilinearForm *a = new BilinearForm(fespace);
    BilinearFormIntegrator *integ = new ElasticityIntegrator(lambda_func, mu_func);
    a->AddDomainIntegrator(integ);

    std::chrono::duration<double> time1;
    auto start1 = std::chrono::system_clock::now();
    //  Assemble the bilinear form and the corresponding linear system,
    //     applying any necessary transformations such as: eliminating boundary
    //     conditions, applying conforming constraints for non-conforming AMR,
    //     static condensation, etc.
    cout << "matrix ... " << flush;
    a->Assemble();
    b->Assemble();
  
    SparseMatrix A;
    Vector B, X;
  
    a->FormLinearSystem(ess_tdof_list, x, *b, A, X, B);  
    cout << "done." << endl << "Size of linear system: " <<
      A.Height() << endl;
  
    if(iterative){
      GSSmoother M(A);
      PCG(A, M, B, X, 2, 10000, 1e-25, 0.0);
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
    // Recover the solution as a finite element grid function.
    a->RecoverFEMSolution(X, *b, x);

    auto end1 = std::chrono::system_clock::now();
    time1 = end1 - start1;
    cout << "Time : " << time1.count()*1000.0 << " ms" << endl;

    // Compute norms of error
    int tdim = dim*(dim+1)/2; // num. entries in a symmetric tensor
    VectorFunctionCoefficient Stress_exactecart_coef(tdim,Stress_exacteCart);
    GridFunction zero(fespace);
    zero = 0.0;
    VectorFunctionCoefficient Stress_exactepol_coef(tdim,Stress_exactePol);
    double ener_refer = ComputeEnergyNorm(zero, lambda_func, mu_func,
					  Stress_exactepol_coef);
    double ener_error = ComputeEnergyNorm(x, lambda_func, mu_func,
					  Stress_exactepol_coef);
    ener_error = ener_error/ener_refer;
    cout << "Energy norm of error: " << ener_error << endl;
    cout << "Numbre of DOFS: "<< fespace->GetNDofs() << endl;
    cout << endl;
    double h = StepMax(*mesh);	//recherche pas max
    //Compute the slope
    slope_ener(iter,0) = time1.count()*1000.0;
    slope_ener(iter,1) = fespace->GetNDofs();
    slope_ener(iter,2) = ener_error;
    err_tmp_ener = ener_error;
    h_tmp = h;
    iter++;
    //Save in Praview format
    //QuadratureFunction Stress;
    //SaveStress(x, lambda_func, mu_func, *mesh, Stress);
    if (!mesh->NURBSext)
      {
        mesh->SetNodalFESpace(fespace);
      }
    //GridFunction *nodes = mesh->GetNodes();
    //*nodes += x;
    /*
      FiniteElementSpace *fieldspace;
      fieldspace = new FiniteElementSpace(mesh, fec, 1, Ordering::byVDIM);
      ParaViewDataCollection paraview_dc("Troue", mesh);
      paraview_dc.SetPrefixPath("ParaView");
      std::string letters = "xyz";
      Array<GridFunction *> stress(dim*(dim+1)/2);
      StressCoefficient stress_c(lambda_func, mu_func);
      stress_c.SetDisplacement(x);
      int c = 0;
      for (int si = 0; si < dim; si++)
      {
      for (int sj = si; sj < dim; sj++)
      {
      std::string stressname= "S" + letters.substr(si,1) +
      letters.substr(sj,1);
      stress[c] = new GridFunction(fieldspace);
      stress_c.SetComponent(si, sj);
      stress[c]->ProjectDiscCoefficient(stress_c, GridFunction::ARITHMETIC);
      paraview_dc.RegisterField(stressname,stress[c]);
      std::string strainname= (si != sj ? "G" : "E") + letters.substr(si,1) +
      letters.substr(sj,1);
      }}
      //Stress exacte
      FiniteElementSpace *fieldspace_e;
      fieldspace_e = new FiniteElementSpace(mesh, fec, tdim, Ordering::byVDIM);
      GridFunction stress_exactepol(fieldspace_e);
      stress_exactepol.ProjectCoefficient(Stress_exactepol_coef);
      GridFunction stress_exactecart(fieldspace_e);
      stress_exactecart.ProjectCoefficient(Stress_exactecart_coef);
      paraview_dc.SetLevelsOfDetail(order+1);
      paraview_dc.SetCycle(0);
      paraview_dc.SetDataFormat(VTKFormat::BINARY);
      paraview_dc.SetHighOrderOutput(true);
      paraview_dc.SetTime(0.0); // set the time
      paraview_dc.RegisterField("numerical_solution",&x);
      paraview_dc.RegisterField("Stress exacte cartésien",&stress_exactecart);
      paraview_dc.RegisterField("Stress exacte polaire",&stress_exactepol);
      paraview_dc.Save();*/

    mesh->UniformRefinement();
  }//Fin de boucle maillage

  //Affichage des normes et pentes.
  cout<<endl;
  cout<<"Erreur en norme:"<<endl;
  for (int i=1; i<iter; i++){
    cout << " Energie: " << slope_ener(i,2)<<" Taille de maille= "
	 <<slope_ener(i,1)<<endl;}

  cout<<endl;
  cout<<"Pente de convergence:"<<endl;
  for (int i=1; i<iter; i++){
    cout << " Energie: " << slope_ener(i,0)<<" Taille de maille= "<<
      slope_ener(i,1)<<endl;}
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
double ComputeAbsoluErr(GridFunction &x, Coefficient &lambdah, Coefficient &muh, double &relative_goal){
  const int dim = x.VectorDim();
  const int tdim = dim*(dim+1)/2;
  VectorFunctionCoefficient Stress_exactecart_coef(tdim,Stress_exacteCart,
						   new ConstantCoefficient(0.0));
  double ener_error = ComputeEnergyNorm(x, lambdah, muh,
					Stress_exactecart_coef);
  return relative_goal*ener_error;
}

//==============Erreur en Norme energie ===================
double ComputeEnergyNorm(GridFunction &x, Coefficient &lambdah, Coefficient &muh,
			 VectorFunctionCoefficient &Stress_exacte_coef)
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

	  energy += w * pdc;
	}
    }
  if(energy>0.0){
    return sqrt(energy);}
  else{
    cout<<"Negative Energy error"<<endl;
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
		   GridFunction &x, Vector &strain){
  FiniteElementSpace *fes = x.FESpace();
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
		   GridFunction &x,Coefficient &lambda, Coefficient &mu, 
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
