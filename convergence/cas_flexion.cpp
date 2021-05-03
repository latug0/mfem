//Cas test d'une poutre en flexion avec solution analytique.
//Calculs avec plusieurs maillage possible

//GL: ajouter le papier de reference et une description plus detaillée
//GL: revoir l'indentation du fichier (en utilisant emacs par exemple).

#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

void sol_exact(const Vector &, Vector &);
void grad_exact(const Vector &, DenseMatrix &);
double F(const Vector &);

double ComputeGradNorm(Mesh &, GridFunction &);

double ComputeEnergyNorm(Mesh &, GridFunction &,
						 Coefficient &, Coefficient &);

//GL: sur deux lignes, évite d'avoir plus de 90 caractère sur une seul ligne
void Elasticy_mat(ElementTransformation &,const IntegrationPoint &, int, Coefficient &, 					  Coefficient &, DenseMatrix &);

//GL: coquille à Exact ...
double Norm_Energie_Eaxct();

int main(int argc, char *argv[])
{
int ref1;
DenseMatrix slope_l2, slope_ener, slope_grad;
double err_tmp_ener = 0., err_tmp_l2=0., err_tmp_grad=0.;
double h_tmp = 0.;
int iter = 0;

   // 1. Parse command-line options.
   const char *mesh_file = "../data/beam-quad.mesh";
   int order;
   bool static_cond = false;
   int rep=5;
   //GL: utiliser args.AddOption pour les paramètres d'entrée, pas "cin"
	cout << "Ordre de la méthode: ";  cin >> order;
	cout << "Pour plusieurs maillages tapez 1: ";  cin >> ref1;
	if (ref1==1){
	cout << "Combien de maillage : ";cin >> rep;
	slope_ener.SetSize(rep-1,3); slope_l2.SetSize(rep-1,3), slope_grad.SetSize(rep-1,3);}

string const err_energy("err_flexion.txt");
    ofstream err_energy_flux(err_energy.c_str());

    if (err_energy_flux)    
    {

for (int ref_levels=1; ref_levels<((rep-2)*ref1+2); ref_levels++){ 

   //GL: utiliser args.AddOption pour les paramètres d'entrée, pas "cin"
if (ref1==0){
	cout << "Combien de rafinement uniforme : "; cin >> ref_levels;}

   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to use.");
   args.AddOption(&order, "-o", "--order",
                  "Finite element order (polynomial degree).");
   args.AddOption(&static_cond, "-sc", "--static-condensation", "-no-sc",
                  "--no-static-condensation", "Enable static condensation.");
   args.Parse();
   if (!args.Good())
   {
      args.PrintUsage(cout);
      return 1;
   }
   args.PrintOptions(cout);

   // 2. Read the mesh from the given mesh file. We can handle triangular,
   //    quadrilateral, tetrahedral or hexahedral elements with the same code.
   Mesh *mesh = new Mesh(mesh_file, 1, 1);
   int dim = mesh->Dimension();

   // 3. Select the order of the finite element discretization space. For NURBS
   //    meshes, we increase the order by degree elevation.
   if (mesh->NURBSext)
   {
      mesh->DegreeElevate(order, order);
   }

   // 4. refine the mesh to increase the resolution. In this example we do
   //    'ref_levels' of uniform refinement. We choose 'ref_levels' to be the
   //    largest number that gives a final mesh with no more than 5,000
   //    elements.
   {
   //   int ref_levels = (int)floor(log(5000./mesh->GetNE())/log(2.)/dim);

      for (int l = 0; l < ref_levels; l++)
      {
         mesh->UniformRefinement();
      }
   }

   // 5. Define a finite element space on the mesh. Here we use vector finite
   //    elements, i.e. dim copies of a scalar finite element space. The vector
   //    dimension is specified by the last argument of the FiniteElementSpace
   //    constructor. For NURBS meshes, we use the (degree elevated) NURBS space
   //    associated with the mesh nodes.
   FiniteElementCollection *fec;
   FiniteElementSpace *fespace;
   if (mesh->NURBSext)
   {
      fec = NULL;
      fespace = mesh->GetNodes()->FESpace();
   }
   else
   {
      fec = new H1_FECollection(order, dim);
      fespace = new FiniteElementSpace(mesh, fec, dim);
   }
   cout << "Number of finite element unknowns: " << fespace->GetTrueVSize()
        << endl << "Assembling: " <<dim<< flush;

   // 6. Determine the list of true (i.e. conforming) essential boundary dofs.
   //    In this example, the boundary conditions are defined by marking only
   //    boundary attribute 1 from the mesh as essential and converting it to a
   //    list of true dofs.

/*
   Array<int> ess_tdof_list, tmp_tdof, ess_bdr(mesh->bdr_attributes.Max());
   cout << "bdr_attributes.Max()" << std::endl;
   ess_bdr = 0;
   // ess_bdr[0] refers to dof on the "xz0" plane
   ess_bdr[0] = 1;
   // grab all dof on "xz0" and put them into tmp_dof.
   // last parameter is 1 to set "y" direction for dirichlet condition.
   fespace->GetEssentialTrueDofs(ess_bdr, tmp_tdof); 
   // ess_tof_list accumulates all needed dof
   ess_tdof_list.Append(tmp_tdof);
*/

  // List of True DoFs : Define (here) Dirichlet conditions
  Array<int> ess_tdof_list;
  Array<int> ess_bdr(mesh->bdr_attributes.Max());
  ess_bdr = 0;
  ess_bdr[0] =1;
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
{
	f.Set(dim-1, new FunctionCoefficient(F));
}

   LinearForm *b = new LinearForm(fespace);
   b->AddBoundaryIntegrator(new VectorBoundaryLFIntegrator(f));
   cout << "r.h.s. ... " << flush;


   // 8. Define the solution vector x as a finite element grid function
   //    corresponding to fespace. Initialize x with initial guess of zero,
   //    which satisfies the boundary conditions.
   GridFunction x(fespace);
   VectorFunctionCoefficient Boundary_Dirichlet_coef(dim, sol_exact);
   x = 0.0;
   // To use if there are different Dirichlet conditions.
   // Beware, the values of dirichlet boundary conditions are set here !
   x.ProjectBdrCoefficient(Boundary_Dirichlet_coef, ess_bdr);

   // 9. Set up the bilinear form a(.,.) on the finite element space
   //    corresponding to the linear elasticity integrator with piece-wise
   //    constants coefficient lambda and mu.

   double lambda;
  double E = 1000.;
  double nu = 0.25;

	lambda = E*nu/((1.+nu)*(1.-2.*nu));

   double mu;
   mu = E/(2.*(1.+nu));
   ConstantCoefficient mu_func(mu);
//lambda = 2*lambda*mu/(lambda+2*mu);
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
   a->RecoverFEMSolution(X, *b, x);


  //GL: Ajouter une option en paramètre d'entrée qui permette de faire itératif ou direct
//#ifndef MFEM_USE_SUITESPARSE
   // 11. Define a simple symmetric Gauss-Seidel preconditioner and use it to
   //     solve the system Ax=b with PCG.
   GSSmoother M(A);
   PCG(A, M, B, X, 1, 50000, 1e-25, 0.0);
/*
#else
   // 11. If MFEM was compiled with SuiteSparse, use UMFPACK to solve the system.
   UMFPackSolver umf_solver;
   umf_solver.Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_METIS;
   umf_solver.SetOperator(A);
   umf_solver.Mult(B, X);
#endif
*/
   // 12. Recover the solution as a finite element grid function.
   a->RecoverFEMSolution(X, *b, x);


   // 13. For non-NURBS meshes, make the mesh curved based on the finite element
   //     space. This means that we define the mesh elements through a fespace
   //     based transformation of the reference element. This allows us to save
   //     the displaced mesh as a curved mesh when using high-order finite
   //     element displacement field. We assume that the initial mesh (read from
   //     the file) is not higher order curved mesh c<ompared to the chosen FE
   //     space.
   if (!mesh->NURBSext)
   {
      mesh->SetNodalFESpace(fespace);
   }


// Compute error
	double ener_error = ComputeEnergyNorm(*mesh, x, lambda_func, mu_func);
	double pdc = Norm_Energie_Eaxct();
	double err_grad = ComputeGradNorm(*mesh, x);
	//cout << "ener exact: "<< abs(pdc - ener_error)<<endl;
	VectorFunctionCoefficient sol_exact_coef(dim, sol_exact);
	double L2_error = x.ComputeL2Error(sol_exact_coef);
	
	cout << "\nL2 norm of error: " << L2_error << endl;
	cout << "Energy norm of error: " << ener_error << endl;
	cout << "Grad norm of error: " << err_grad << endl;
	cout << "numbers of elements: " << mesh->GetNE() <<endl;

if (ref1==1){
	double h = mesh->GetElementSize(0);
	//Compute the slope
	slope_l2(iter,0) = log(err_tmp_l2/L2_error) / log(h_tmp/h);
	slope_l2(iter,1) = h;
	slope_l2(iter,2) = L2_error;
	err_tmp_l2 = L2_error;
	//Compute the slope
	slope_ener(iter,0) = log(err_tmp_ener/ener_error) / log(h_tmp/h);
	slope_ener(iter,1) = h;
	slope_ener(iter,2) = ener_error;
	err_tmp_ener = ener_error;
	//Compute the slope
	slope_grad(iter,0) = log(err_tmp_grad/err_grad) / log(h_tmp/h);
	slope_grad(iter,1) = h;
	slope_grad(iter,2) = err_grad;
	err_tmp_grad = err_grad;
	h_tmp = h;

	//Save in errors
	//col1: elmentSize		  col2: L2 norm error	col3: ener norm error	
	//col6: slope L2 error 	  col7: slope Ener error	col8: écart	
	 err_energy_flux <<h<<" "<<L2_error<<" "<<ener_error <<" "<<err_grad<<" "<<slope_l2(iter,0)<<" "<<slope_ener(iter,0)<<" "<< pdc <<endl;
	iter++;
}

//Save in Praview format
if (ref1==0){
	GridFunction ex(fespace);
	ex.ProjectCoefficient(sol_exact_coef);
	GridFunction diff(fespace);
	diff.ProjectCoefficient(sol_exact_coef);
	diff -= x;
 
	ParaViewDataCollection paraview_dc("Flexion", mesh);
	paraview_dc.SetPrefixPath("ParaView");
	paraview_dc.SetLevelsOfDetail(order);
	paraview_dc.SetCycle(0);
	paraview_dc.SetDataFormat(VTKFormat::BINARY);
	paraview_dc.SetHighOrderOutput(true);
	paraview_dc.SetTime(0.0); // set the time
	paraview_dc.RegisterField("numerical_solution",&x);
	paraview_dc.RegisterField("diff-exact_solution",&diff);
	paraview_dc.RegisterField("exact_solution",&ex);
	paraview_dc.Save();	
}

//end loop mesh
}

if (ref1==1){
	cout<<endl;
	cout<<"Erreur en norme:"<<endl;
	for (int i=1; i<iter; i++)
			cout << "L2: " << slope_l2(i,2)<<" Grad: "<<slope_grad(i,2) << " Energie: " << slope_ener(i,2)<<" Taille de maille= "<<slope_l2(i,1)<<endl;
	cout<<endl;
	cout<<"Pente de convergence:"<<endl;
	for (int i=1; i<iter; i++)
			cout << "Pente L2: " << slope_l2(i,0)<<" Grad: "<<slope_grad(i,0) << " Energie: " << slope_ener(i,0)<<" Taille de maille= "<<slope_l2(i,1)<<endl;
	cout<<endl;
}

}
    else
    {
        cout << "ERREUR: Impossible d'ouvrir le fichier." << endl;
    }
return 0;
}

// Definition of exact solution
//GL: commenter les formules et donner la référence/papier ou les trouver
void sol_exact(const Vector &x, Vector &u)
{
//GL: Il faut que ces données pull_force, L, D, E, nu, I soit communes à toutes les routines.
//GL: Une solution simple (mais pas idéale, on verra cela dans un second temps) est de faire des variables globales.
//GL: Voir ex22p.cpp par exemple, les constantes sont définies au début comme des "static double".
//GL: Si elles sont connues à la compilation, tu peux les définir comme des "constexpr double".
  double pull_force = -1.;
  double L = 8.0;
  double D =1.0;
  double E = 1000.;
  double nu = 0.25;
  double I = D*D*D*D/12.;
	double y = x(1)-D*0.5;
//GL: expliquer ces transformations de manière adéquate.
E=E/(1.-nu*nu); nu = nu/(1.-nu);
 //u(0) = -pull_force*y/(6.*E*I) * ((6.*L-3.*x(0))*x(0) + y*y*(2.+nu)- 1.5*D*D*(1.+nu));
 //u(1) = pull_force/(6.*E*I) * (3.*nu*y*y*(L-x(0)) + (3*L-x(0))*x(0)*x(0));

	u(0) = -pull_force*y/(6.*E*I) * ((6.*L-3.*x(0))*x(0) + (2.+nu)*(y*y - D*D/4.));
	u(1) = pull_force/(6.*E*I) * (3.*nu*y*y*(L-x(0)) + (3.*L-x(0))*x(0)*x(0)+ (4.+5.*nu)*D*D*x(0)/4.);
}

//GL: commenter les formules et donner la référence/papier ou les trouver
void grad_exact(const Vector &x, DenseMatrix &grad)
{
	double pull_force = -1;
	double L = 8.0;
	double D =1.0;
	double E = 1000.;
	double nu = 0.25;
	double I = D*D*D/12.;
	double y = x(1)-D*0.5;

	E=E/(1.-nu*nu); nu = nu/(1.-nu);

	grad(0,0) = -pull_force*y/(E*I)*(L-x(0));
	grad(0,1) = -pull_force/(6.*E*I)*((6*L-3*x(0))*x(0) + (2+nu)*(y*y-D*D/4) + 2*y*y*(2+nu));
	grad(1,1) = pull_force*y*nu/(E*I)*(L-x(0));
	grad(1,0) = pull_force/(6.*E*I)*(2*x(0)*(3*L-x(0))-x(0)*x(0)-3*nu*y*y+(4+5*nu)*D*D/4);
}
// load at the end of the beam  
double F(const Vector &x)
{
  double pull_force = -1.;
  double L = 8.0;
  double D =1.0;
	double I = D*D*D*D/12.;
	double force;
	if(x(0) >= 8.-1.e-6){
		force = pull_force/(2.*I)*(pow(D*0.5,2) - pow((x(1)-D*0.5),2));}
	else{
		force = 0.;}

	return force;
}


double ComputeGradNorm(Mesh &mesh, GridFunction &x){
	FiniteElementSpace *fes = x.FESpace();
	int dim = fes->GetMesh()->SpaceDimension();

	MatrixFunctionCoefficient grad_exact_coef (dim, grad_exact);
 
	ElementTransformation *Trans;
	Array<int> vdofs;
	DenseMatrix grad, gradh;
	L2_FECollection flux_fec(fes->GetOrder(1), dim);
	FiniteElementSpace flux_fes(&mesh, &flux_fec, dim);
	double error = 0.0;
	Array<int> udofs;
	Vector ul;
	for (int i = 0; i < fes->GetNE() ; i++)
	{
		const FiniteElement *fe = fes->GetFE(i);
		const FiniteElement *flux_fe = flux_fes.GetFE(i);
		const int order = 2*fe->GetOrder() + 3;   //<----------
		const IntegrationRule *ir = &(IntRules.Get(fe->GetGeomType(), order));
		Trans = fes->GetElementTransformation(i);
		const int dof = fe->GetDof();
		const int dim = fe->GetDim();
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

   		//DenseMatrix loc_data_mat(ul.GetData(), dof, dim);
	for (int j = 0; j < ir->GetNPoints(); j++)
	{
		const IntegrationPoint &ip = ir->IntPoint(j);
		Trans->SetIntPoint(&ip);
		double w = Trans->Weight() * ip.weight;
		fe->CalcDShape(ip, dshape);
		MultAtB(loc_data, dshape, gh);
		Mult(gh, Trans->InverseJacobian(), gradh);
		grad_exact_coef.Eval(grad,*Trans,ip);
/*
	cout<<endl;
	cout<<endl;
	cout<<"Intergration poin: "<<j<<" element: "<<i<<endl;
	cout<<endl;
	cout<<"Coord x,y: "<<ip.x<<" "<<ip.y<<endl;
	cout<<grad(0,0)<<" "<<gradh(0,0)<<" écart entre les composantes "<<grad(0,0) - gradh(0,0)<<endl;
	cout<<grad(1,1)<<" "<<gradh(1,1)<<" écart entre les composantes "<<grad(1,1) - gradh(1,1)<<endl;
	cout<<grad(0,1)<<" "<<gradh(0,1)<<" écart entre les composantes "<<grad(0,1) - gradh(0,1)<<endl;
	cout<<grad(1,0)<<" "<<gradh(1,0)<<" écart entre les composantes "<<grad(1,0) - gradh(1,0)<<endl;
*/
		grad -= gradh;
		error += w * grad(0,0)*grad(0,0) + grad(1,1)*grad(1,1) + grad(0,1)*grad(0,1) + grad(1,0)*grad(1,0);
	}			
	}
return (error < 0.0) ? -sqrt(-error) : sqrt(error);
}

double ComputeEnergyNorm(Mesh &mesh, GridFunction &x,
						 Coefficient &lambdah, Coefficient &muh)
{
	FiniteElementSpace *fes = x.FESpace();
	int dim = fes->GetMesh()->SpaceDimension();
	VectorFunctionCoefficient sol_exact_coef (dim, sol_exact);
	GridFunction ex(fes);
	ex.ProjectCoefficient(sol_exact_coef);
  
	double E = 1000.;
	double nu = 0.25;
	E=E/(1.-nu*nu); nu = nu/(1.-nu);
	double lambda_d;
	lambda_d = E*nu/((1.+nu)*(1.-2.*nu));
	ConstantCoefficient lambda(lambda_d);
	double mu_d;
	mu_d = E/(2.*(1.+nu));
	ConstantCoefficient mu(mu_d);
	BilinearFormIntegrator *integ = new ElasticityIntegrator(lambda, mu);
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
	for (int j = 0; j < ir->GetNPoints(); j++)
    {
		const IntegrationPoint &ip = ir->IntPoint(j);
		Trans->SetIntPoint(&ip);
		double w = Trans->Weight() * ip.weight;

	  //======= Strains vectors ========
		fe->CalcShape(ip, shape);
		fe->CalcShape(ip, shapeh);

		flux_mat.MultTranspose(shape, strain);
		flux_math.MultTranspose(shapeh, strainh);

	  //======= Stress vectors ========
		Elasticy_mat(*Trans,ip,dim,lambda,mu,C);
		Elasticy_mat(*Trans,ip,dim,lambdah,muh,Ch);

		C.Invert();	Ch.Invert();
		Ch.Mult(strainh,stressh);	//approx
		C.Mult(strain,stress);	//exacte

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

// Definition of exact solution
double Norm_Energie_Eaxct()
{
	double P = 2.;
  	double L = 8.0;
	double D =1.0;
	double E = 1000.;
	double nu = 0.25;
	double I = D*D*D/12.;
	//E=E/(1.-nu*nu); nu = nu/(1.-nu);
	double lambda;
	lambda = E*nu/((1.+nu)*(1.-2.*nu));
	double mu;
	mu = E/(2.*(1.+nu));

	double pdc;
	pdc = pow(P/I,2) * ((lambda+2.*mu)/(mu*(4.*lambda +mu))*L*L*L*D*D*D/36. 
					 - D*D*D*L/(12.*mu));
return (pdc < 0.0) ? -sqrt(-pdc) : sqrt(pdc);

}

