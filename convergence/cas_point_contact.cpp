//
// Compile with: make cas_point_contact
//
//
//GL: ajouter commentaire descriptif


//GL: indenter correctement le fichier. Sous emacs
//GL: https://unix.stackexchange.com/questions/162057/can-i-have-emacs-automatically-indent-my-whole-code-after-it-is-all-written

#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

void sol_exact(const Vector &, Vector &);

double ComputeEnergyNorm(Mesh &, GridFunction &,
			 Coefficient &, Coefficient &);

void Elasticy_mat(ElementTransformation &,const IntegrationPoint &, int, Coefficient &,
		  Coefficient &, DenseMatrix &);


void Grad(ElementTransformation &,const IntegrationPoint &,
	  const GridFunction &, DenseMatrix &);

int main(int argc, char *argv[])
{
  //Variables pour l'affichage
  int ref1;
  DenseMatrix slope_l2, slope_ener;
  double err_tmp_ener = 0., err_tmp_l2=0.;
  double h_tmp = 0.;
  int iter = 0;
   
	const char *mesh_file = "carre.msh";
	bool static_cond = false;
	int rep=5;

	int order;
	//GL: utiliser "args.AddOption" utilisé ci-après pour spécifier 'order' et 'rep'
	//GL: une valeur sera donnée par défaut, par exemple order=1, ref1=0, rep=4
	cout << "Ordre de la méthode: ";  cin >> order;
	cout << "Pour plusieurs maillages tapez 1: ";  cin >> ref1;
	if (ref1==1){
	cout << "Combien de maillage : ";cin >> rep;
	slope_ener.SetSize(rep-1,3); slope_l2.SetSize(rep-1,3);}

	string const err_energy("err_contact.txt");		//Fichier avec les erreurs
    ofstream err_energy_flux(err_energy.c_str());
    if (err_energy_flux)    
    {
//Boucle pour maillage plus fin
for (int ref_levels=1; ref_levels<((rep-2)*ref1+2); ref_levels++){ 

if (ref1==0){
	//GL: utiliser "args.AddOption" utilisé ci-après pour spécifier les valeurs,
        //GL: n'utilise pas de "cin" STP
	cout << "Combien de rafinement uniforme : "; cin >> ref_levels;}

	//GL: Le traitement des options doit être remonté au-dessus de la boucle ref_level
	// Parse command-line options.
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

	//Lecture du malliage
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


	double E = 1000.;
	double nu = 0.25;
	double lambda;
	lambda = E*nu/((1.+nu)*(1.-2.*nu));
	double mu;
	mu = E/(2.*(1.+nu));
	ConstantCoefficient mu_func(mu);
	//GL: pas nécessaire de conserver cette formule de lambda, non ?
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

   //GL: utiliser une option en ligne de commande pour activer soit suitesparse, soit PCG
   //GL: il faudrait ne pas avoir à recompiler pour changer de solveur 
   //#ifndef MFEM_USE_SUITESPARSE
   // 11. Define a simple symmetric Gauss-Seidel preconditioner and use it to
   //     solve the system Ax=b with PCG.
   GSSmoother M(A);
   PCG(A, M, B, X, 1, 5000, 1e-11, 0.0);
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


// Compute errors
	double ener_error = ComputeEnergyNorm(*mesh, x, lambda_func, mu_func);
	VectorFunctionCoefficient sol_exact_coef(dim, sol_exact);
	double L2_error = x.ComputeL2Error(sol_exact_coef);
	
	cout << "\nL2 norm of error: " << L2_error << endl;
	cout << "Energy norm of error: " << ener_error << endl;
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
	h_tmp = h;

	//Save in errors
	//col1: elmentSize		  col2: L2 norm error	col3: ener norm error	
	// col4: logelementSize	  col5: logL2Error		col6: log(enererror)
	//col6: slope L2 error 	  col7: slope Ener error	col8: écart	
	 err_energy_flux <<h<<" "<<L2_error<<" "<<ener_error <<" "<< log(h)<<" "<<log(L2_error)<<" "<<log(ener_error)<<" "<<slope_l2(iter,0)<<" "<<slope_ener(iter,0)<<endl;
	iter++;
   //  Free the used memory.
   delete a;
   delete b;
   if (fec)
   {
      delete fec;
   }
   delete mesh;
}

//Save in Praview format
if (ref1==0){
	GridFunction diff(fespace);
	GridFunction ex1(fespace);
	diff.ProjectCoefficient(sol_exact_coef);
	ex1.ProjectCoefficient(sol_exact_coef);
	//GL: Regarder si dans les exemples MFEM une solution plus classique ne serait pas utilisée pour cette différence
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

//end loop mesh
}

//Affichage
if (ref1==1){
	cout<<endl;
	cout<<"Erreur en norme:"<<endl;
	for (int i=1; i<iter; i++)
			cout << "L2: " << slope_l2(i,2) << " Energie: " << slope_ener(i,2)<<" Taille de maille= "<<slope_l2(i,1)<<endl;
	cout<<endl;
	cout<<"Pente de convergence:"<<endl;
	for (int i=1; i<iter; i++)
			cout << "Pente L2: " << slope_l2(i,0) << " Energie: " << slope_ener(i,0)<<" Taille de maille= "<<slope_l2(i,1)<<endl;
	cout<<endl;
}

}
    else
    {
        cout << "ERREUR: Impossible d'ouvrir le fichier." << endl;
    }
return 0;
}

//===================== Solution exacte =====================
//GL: Ajouter la référence
//GL: Test/vérifier ce que donne sol_exact pour x={0, 0}.
void sol_exact(const Vector &x, Vector &u)
{
  double P = 2.;
  double E = 1000.;
  double nu = 0.25;
  double r2 = x(1)*x(1)+x(0)*x(0);
  double pi = M_PI;
	double lambda;
    lambda = E*nu/((1.+nu)*(1.-2.*nu));
	double mu;
	mu = E/(2.*(1.+nu));
	u(0) = -P/(4*pi*mu)*(2*x(0)*x(1)/r2 + 2*mu/(lambda+mu)*atan2(x(1),x(0)));
	u(1) = -P/(4*pi*mu)*((x(1)*x(1)-x(0)*x(0))/r2 - (lambda+2*mu)*log(r2)/(lambda+mu));
}

//===================== Erreur Norme Energie =====================
//GL: Ajouter de commentaires STP
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
	//E=E/(1.-nu*nu); nu = nu/(1.-nu);
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


//GL: indiquer comment tu as trouvé cette matrice élasticité, quel algo utilises-tu
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

