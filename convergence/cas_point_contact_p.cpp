// Compile with: make cas_point_contact
//
//Cas test avec solution analytique. Plaque 2D avec point de contact en (0,0).
//Calcul de la déformation en élasticité linéaire. 
//Version parallèle

//Ref: http://e6.ijs.si/medusa/wiki/index.php/Point_contact


#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;
static double Force = -1.;
static double E = 1000.;
static double nu = 0.25;

void sol_exact(const Vector &, Vector &);

double ComputeEnergyNorm(Mesh &, GridFunction &,
			 Coefficient &, Coefficient &);

void Elasticy_mat(ElementTransformation &,const IntegrationPoint &, int, Coefficient &,
		  Coefficient &, DenseMatrix &);


void Grad(ElementTransformation &,const IntegrationPoint &,
	  const GridFunction &, DenseMatrix &);

int main(int argc, char *argv[])
{
  // Initialize MPI.
  int num_procs, myid;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
   
  // Parse command-line options.
  int order=1;
  const char *mesh_file = "carre.msh";
  bool static_cond = false;
  bool amg_elast = 0;
  bool reorder_space = false;
  bool solver=true;
  OptionsParser args(argc, argv);
  args.AddOption(&mesh_file, "-m", "--mesh",
		 "Mesh file to use.");
  args.AddOption(&order, "-o", "--order",
		 "Finite element order (polynomial degree).");
  args.AddOption(&amg_elast, "-elast", "--amg-for-elasticity", "-sys",
		 "--amg-for-systems",
		 "Use the special AMG elasticity solver (GM/LN approaches), "
		 "or standard AMG for systems (unknown approach).");
  args.AddOption(&static_cond, "-sc", "--static-condensation", "-no-sc",
		 "--no-static-condensation", "Enable static condensation.");
  args.AddOption(&reorder_space, "-nodes", "--by-nodes", "-vdim", "--by-vdim",
		 "Use byNODES ordering of vector space instead of byVDIM");
  args.AddOption(&solver, "-sol", "--Itératif", "-Direct",
		 "--Solver Itératif", "Solver direct.");

  int ref_levels = 9;
  //cout << "Combien de rafinement uniforme : "; cin >> ref_levels;

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

  //Lecture du malliage
  Mesh *mesh = new Mesh(mesh_file, 1, 1);
  int dim = mesh->Dimension();

  //  Select the order of the finite element discretization space. For NURBS
  //    meshes, we increase the order by degree elevation.
  if (mesh->NURBSext)
    {
      mesh->DegreeElevate(order, order);
    }

  //  refine the mesh to increase the resolution. In this example we do
  //    'ref_levels' of uniform refinement. We choose 'ref_levels' to be the
  //    largest number that gives a final mesh with no more than 5,000
  //    elements.
  //   int ref_levels = (int)floor(log(5000./mesh->GetNE())/log(2.)/dim);

  for (int l = 0; l < ref_levels; l++)
    {
      mesh->UniformRefinement();
    }

  //Define parallel mesh by a partitioning of the serial mesh.
  ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
  {
    int par_ref_levels = 1;
    for (int l = 0; l < par_ref_levels; l++)
      {
	pmesh->UniformRefinement();
      }
  }

  //  Define a finite element space on the mesh. Here we use vector finite
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
	  fespace = new ParFiniteElementSpace(pmesh, fec, dim, Ordering::byNODES);
	}
      else
	{
	  fespace = new ParFiniteElementSpace(pmesh, fec, dim, Ordering::byVDIM);
	}
    }
  HYPRE_Int size = fespace->GlobalTrueVSize();
  if (myid == 0)
    {
      cout << "Number of finite element unknowns: " << size << endl
           << "Assembling: " << flush;
    }

  // Determine the list of true (i.e. conforming) essential boundary dofs.
  //    In this example, the boundary conditions are defined by marking only
  //    boundary attribute 1 from the mesh as essential and converting it to a
  //    list of true dofs.


  // List of True DoFs : Define (here) Dirichlet conditions
  Array<int> ess_tdof_list;
  Array<int> ess_bdr(pmesh->bdr_attributes.Max());
  ess_bdr = 1;		//On sélectionne tous les bords
  fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

  // Set up the linear form b(.) which corresponds to the right-hand side of
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

  ParLinearForm *b = new ParLinearForm(fespace);
  b->Assemble();
  if (myid == 0)
    {
      cout << "r.h.s. ... " << flush;
    }

  // Define the solution vector x as a finite element grid function
  //    corresponding to fespace. Initialize x with initial guess of zero,
  //    which satisfies the boundary conditions.
  ParGridFunction x(fespace);
  VectorFunctionCoefficient Boundary_Dirichlet_coef(dim, sol_exact);
  x = 0.0;
  // To use if there are different Dirichlet conditions.
  // Beware, the values of dirichlet boundary conditions are set here !
  //Projète solution exacte sur les bords
  x.ProjectBdrCoefficient(Boundary_Dirichlet_coef, ess_bdr);	

  // Set up the bilinear form a(.,.) on the finite element space
  //    corresponding to the linear elasticity integrator with piece-wise
  //    constants coefficient lambda and mu.

  double lambda;
  lambda = E*nu/((1.+nu)*(1.-2.*nu));
  double mu;
  mu = E/(2.*(1.+nu));
  ConstantCoefficient mu_func(mu);
  ConstantCoefficient lambda_func(lambda);

  ParBilinearForm *a = new ParBilinearForm(fespace);
  BilinearFormIntegrator *integ = new ElasticityIntegrator(lambda_func, mu_func);
  a->AddDomainIntegrator(integ);

  // Assemble the parallel bilinear form and the corresponding linear
  //     system, applying any necessary transformations such as: parallel
  //     assembly, eliminating boundary conditions, applying conforming
  //     constraints for non-conforming AMR, static condensation, etc.
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

  // 13. Define and apply a parallel PCG solver for A X = B with the BoomerAMG
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
    pcg->SetTol(1e-15);
    pcg->SetMaxIter(5000);
    pcg->SetPrintLevel(2);
    pcg->SetPreconditioner(*amg);
    pcg->Mult(B, X);
  }
  else{
    cout<<"Solver direct non implémenté<"endl;
  }   


  // 12. Recover the solution as a finite element grid function.
  a->RecoverFEMSolution(X, *b, x);


  // 13. For non-NURBS meshes, make the mesh curved based on the finite element
  //     space. This means that we define the mesh elements through a fespace
  //     based transformation of the reference element. This allows us to save
  //     the displaced mesh as a curved mesh when using high-order finite
  //     element displacement field. We assume that the initial mesh (read from
  //     the file) is not higher order curved mesh c<ompared to the chosen FE
  //     space.
  if (!use_nodal_fespace)
    {
      pmesh->SetNodalFESpace(fespace);
    }

  // Compute errors
  //double ener_error = ComputeEnergyNorm(*mesh, x, lambda_func, mu_func);
  VectorFunctionCoefficient sol_exact_coef(dim, sol_exact);
  double L2_error = x.ComputeL2Error(sol_exact_coef);
  if (myid == 0)
    {
      cout << "\nL2 norm of error: " << L2_error << endl;
      //cout << "Energy norm of error: " << ener_error << endl;
      cout << "numbers of elements: " << pmesh->GetNE() <<endl;
    }
  //  Free the used memory.
  delete pcg;
  delete amg;
  delete a;
  delete b;
  if (fec)
    {
      delete fec;
    }
  delete mesh;
  delete pmesh;
  MPI_Finalize();
  /*
  //Save in Praview format
  if (ref1==0){
  GridFunction diff(fespace);
  GridFunction ex1(fespace);
  diff.ProjectCoefficient(sol_exact_coef);
  ex1.ProjectCoefficient(sol_exact_coef);
  diff-= x;
 
  ParaViewDataCollection paraview_dc("Example2", mesh);
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
  }
  */

  return 0;
}

//===================== Solution exacte =====================
void sol_exact(const Vector &x, Vector &u)
{
  double r2 = x(1)*x(1)+x(0)*x(0);
  double pi = M_PI;
  double lambda;
  lambda = E*nu/((1.+nu)*(1.-2.*nu));
  double mu;
  mu = E/(2.*(1.+nu));
  u(0) = -Force/(4*pi*mu)*(2*x(0)*x(1)/r2 + 2*mu/(lambda+mu)*atan2(x(1),x(0)));
  u(1) = -Force/(4*pi*mu)*((x(1)*x(1)-x(0)*x(0))/r2 - (lambda+2*mu)*log(r2)/(lambda+mu));
}

//===================== Erreur Norme Energie =====================
double ComputeEnergyNorm(Mesh &mesh, GridFunction &x,
			 Coefficient &lambdah, Coefficient &muh)
{
  FiniteElementSpace *fes = x.FESpace();
  int dim = fes->GetMesh()->SpaceDimension();
  VectorFunctionCoefficient sol_exact_coef (dim, sol_exact);
  GridFunction ex(fes);
  ex.ProjectCoefficient(sol_exact_coef);
  
  BilinearFormIntegrator *integh = new ElasticityIntegrator(lambdah, muh);

  ElementTransformation *Trans;
  Array<int> vdofs;

  // compute sigma(ex) and sigma(x)
  // then compute ||sigma(ex) - sigma(x)||
  FiniteElementSpace *ufes = ex.FESpace();
  FiniteElementSpace *uhfes = x.FESpace();
  Vector sigma, sigma_h;
  L2_FECollection flux_fec(fes->GetOrder(1), dim);
  FiniteElementSpace flux_fes(&mesh, &flux_fec, dim);
  double energy = 0.0;
  for (int i = 0; i < ufes->GetNE() ; i++)
    {
      const FiniteElement *fe = ufes->GetFE(i);
      const int order = 2*fe->GetOrder()+3; // <----------
      const IntegrationRule *ir = &IntRules.Get(fe->GetGeomType(), order);

      Vector u, uh;
      // compute sigma(ex)
      ufes->GetElementVDofs(i, vdofs);
      ex.GetSubVector(vdofs, u);
      Trans = ufes->GetElementTransformation(i);
      integ->ComputeElementFlux(*ufes->GetFE(i), *Trans, u, *flux_fes.GetFE(i), sigma);

      // compute sigma(x)
      uhfes->GetElementVDofs(i, vdofs);
      x.GetSubVector(vdofs, uh);
      Trans = uhfes->GetElementTransformation(i);
      integh->ComputeElementFlux(*uhfes->GetFE(i), *Trans, uh, *flux_fes.GetFE(i), sigma_h);
		
      const int dof = fe->GetDof();
      const int dim = fe->GetDim();
      const int tdim = dim*(dim+1)/2; // num. entries in a symmetric tensor

      // View of the 'flux' vector as a (dof x tdim) matrix
      DenseMatrix flux_mat(sigma.GetData(), dof, tdim);
      DenseMatrix flux_math(sigma_h.GetData(), dof, tdim);

      Vector shape(dof),  shapeh(dof);
      Vector stressh(tdim), strainh(tdim);	//approché
      Vector stress(tdim), strain(tdim);	//exacte
      DenseMatrix C,Ch;
      for (int i = 0; i < ir->GetNPoints(); i++)
	{
	  const IntegrationPoint &ip = ir->IntPoint(i);
	  Trans->SetIntPoint(&ip);
	  double w = Trans->Weight() * ip.weight;

	  //======= Strains vectors ========
	  fe->CalcShape(ip, shape);
	  fe->CalcShape(ip, shapeh);

	  flux_mat.MultTranspose(shape, strain);
	  flux_math.MultTranspose(shapeh, strainh);

	  //======= Stress vectors ========
	  Elasticy_mat(*Trans,ip,dim,lambdah,muh,Ch);
	  Ch.Invert();
	  Ch.Mult(strainh,stressh);	//approx
	  Ch.Mult(strain,stress);	//exacte

	  strainh -= strain;
	  stressh -= stress;

	  double pdc=0.0;
	  for (int k = 0; k< dim; k++)
	    pdc += strainh(k)*stressh(k);

	  for (int k = dim; k < dim*(dim+1)/2; k++)
	    pdc += 2*strainh(k)*stressh(k);

	  energy += w * pdc;
	}
    }
  return (energy < 0.0) ? -sqrt(-energy) : sqrt(energy);
}


//===================== Matrice élasticité =====================
void Elasticy_mat(ElementTransformation &T,const IntegrationPoint &ip, 
		  int dim, Coefficient &lambda, Coefficient &mu, DenseMatrix &C){
  double M = mu.Eval(T, ip);
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

