//                                      
// Compile with: make test_full_periodic
//
// Author: Guillaume Latu
//
// Description:  This code shows how to export stress diagnostics.
//

#include "mfem.hpp"
#include "fem/bilininteg.hpp"
#include <fstream>
#include <iostream>
#include "general/forall.hpp"
#include "fem/bilininteg.hpp"
#include "fem/gridfunc.hpp"
#include "fem/qfunction.hpp"

using namespace std;
using namespace mfem;

// central frontier in-between the two materials
const double xthr = 0.5;

void sol_exact1(const Vector &x, Vector &u);
void sol_exact2(const Vector &x, Vector &u);
void sol_exact3(const Vector &x, Vector &u);
void sol_exact4(const Vector &x, Vector &u);
void sol_exact5(const Vector &x, Vector &u);
void sol_exact6(const Vector &x, Vector &u);
void ComputeAvg(const GridFunction &gf, Vector &avg);
void ComputeAvgPerAtt(const GridFunction &gf, Vector &avg, Vector &vol);
void ComputeAvgStressStrain(const GridFunction &u, Coefficient &lambda, Coefficient &mu,
			    int si, int sj, Vector &strain, Vector &stress, Vector &vol);



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

class StrainCoefficient : public DiagCoefficient
{
public:
  using DiagCoefficient::DiagCoefficient;
  void SetComponent(int i, int j) { si = i; sj = j; }
  virtual double Eval(ElementTransformation &T, const IntegrationPoint &ip);
};


/// Class for domain integrator L(v) := (f, grad v)
class DomainLFGrad2Integrator : public LinearFormIntegrator
{
protected:
  Coefficient *lambda, *mu;
private:
  Vector shape, Qvec, Svec;
  DenseMatrix dshape;
  int tcase;

public:
   /// Constructs the domain integrator (Q, grad v)
  DomainLFGrad2Integrator(Coefficient &l, Coefficient &m, int _tcase)
  {
    lambda = &l; mu = &m; tcase = _tcase;
  }

  void AssembleRHSElementVect(
    const FiniteElement &el, ElementTransformation &Tr, Vector &elvect)
  {
    int dof = el.GetDof(); 
    int spaceDim = Tr.GetSpaceDim();
    
    dshape.SetSize(dof, spaceDim);
    Qvec.SetSize(3);
    Svec.SetSize(3);
    elvect.SetSize(dof*spaceDim);
    elvect = 0.0;
   const IntegrationRule *ir = IntRule ? IntRule :  
     &IntRules.Get(el.GetGeomType(), 2 * Tr.OrderGrad(&el));
    
    for (int i = 0; i < ir->GetNPoints(); i++)
      {
	const IntegrationPoint &ip = ir->IntPoint(i);
	Tr.SetIntPoint(&ip);
	el.CalcPhysDShape(Tr, dshape);

	double M = mu->Eval(Tr, ip);
	double L = lambda->Eval(Tr, ip);
	switch(tcase) {
	case 1:
	  Svec[0]=0;Qvec[0] = -(L+2*M);
	  Svec[1]=1;Qvec[1] = -L;
	  Svec[2]=2;Qvec[2] = -L;
	  break;
	case 2:
	  Svec[0]=0;Qvec[0] = -L;
	  Svec[1]=1;Qvec[1] = -(L+2*M);
	  Svec[2]=2;Qvec[2] = -L;
	  break;
	case 3:
	  Svec[0]=0;Qvec[0] = -L;
	  Svec[1]=1;Qvec[1] = -L;
	  Svec[2]=2;Qvec[2] = -(L+2*M);
	  break;
	case 4:
	  Svec[0]=1;Qvec[0] = -M;
	  Svec[1]=0;Qvec[1] = -M;
	  Svec[2]=2;Qvec[2] = 0;
	  break;
	case 5:
	  Svec[0]=0;Qvec[0] = 0;
	  Svec[1]=2;Qvec[1] = -M;
	  Svec[2]=1;Qvec[2] = -M;
	  break;
	case 6:
	  Svec[0]=2;Qvec[0] = -M;
	  Svec[1]=1;Qvec[1] = 0;
	  Svec[2]=0;Qvec[2] = -M;
	  break;
	}
	Qvec *= ip.weight * Tr.Weight();
	for (int k = 0; k < spaceDim; k++) 
	  {	
	    for (int s = 0; s < dof; s++)
	      {
		elvect(dof*k+s)	+= Qvec[k]*dshape(s,Svec[k]);
	      }
	  }
      }
  }
  
  
};



int main(int argc, char *argv[])
{
   //    Parse command-line options.
   const char *mesh_file = "square_2mat_per.msh";
   bool pa = false;
   bool ea = false;
   bool postproc = true;
   int ref = 0;
   int order = 1;
   int tcase = 1;
   bool static_cond = false;
   const char *device_config = "cpu";
   bool visualization = 1;

   OptionsParser args(argc, argv);
   args.AddOption(&pa, "-pa", "--partial-assembly", "-no-pa",
                  "--no-partial-assembly", "Enable Partial Assembly.");
   args.AddOption(&ea, "-ea", "--element-assembly", "-no-ea",
                  "--no-element-assembly", "Enable Element Assembly.");
   args.AddOption(&postproc, "-po", "--postproc", "-no-po",
                  "--no-postproc", "Enable prosprocessing.");
   args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to use.");
   args.AddOption(&order, "-o", "--order",
                  "Finite element order (polynomial degree).");
   args.AddOption(&tcase, "-t", "--tcase",
                  "identifier of the case : Exx->1, Eyy->2, Ezz->3, Exy->4, Eyz->5, Exz->6");
   args.AddOption(&ref, "-r", "--refine",
                  "Nb of refine steps");
   args.AddOption(&device_config, "-d", "--device",
                  "Device configuration string, see Device::Configure().");
   args.Parse();
   if (!args.Good())
   {
      args.PrintUsage(cout);
      return 1;
   }
   args.PrintOptions(cout);

   Device device(device_config);
   device.Print();
   if (!strcmp(device_config, "cuda") && (!pa && !ea))
     MFEM_ABORT("CUDA device requires partial or element assembly")
       
   //   Read the mesh from the given mesh file. We can handle triangular,
   //    quadrilateral, tetrahedral or hexahedral elements with the same code.
   Mesh *mesh;
   std::string mesh_name(mesh_file);
   size_t lastindex = mesh_name.rfind(".med");
   const bool medconfig = lastindex != string::npos;
   if (medconfig){
#ifdef MFEM_USE_MED
     std::string per_name = mesh_name.substr(0, lastindex).append("_per.txt"); 
     mesh = new Mesh();
     mesh->ImportMED(mesh_file, "Box", 0, per_name);
#else /* MFEM_USE_MED */
     MFEM_VERIFY(0, "MED support not activated\n");
#endif /* MFEM_USE_MED */
   } else {
     mesh = new Mesh(mesh_file, 1, 1);
   }
   mesh->CheckElementOrientation(1);

   int dim = mesh->Dimension();

   if ( dim == 2 && tcase != 1 && tcase != 2 && tcase != 4 ) {
     cout << "This test case is undefined in 2D" << endl;
     exit(1);
   }
   if ( dim == 3 && ((tcase < 1) || (tcase > 6))) {
     cout << "This test case is undefined in 3D" << endl;
     exit(1);
   }

   //  Refine the mesh.
   int ref_levels = ref;
   for (int l = 0; l < ref_levels; l++)
     {
       mesh->UniformRefinement();
     }

   //    Define a finite element space on the mesh. Here we use vector finite
   //    elements, i.e. dim copies of a scalar finite element space. The vector
   //    dimension is specified by the last argument of the FiniteElementSpace
   //    constructor. 
   FiniteElementCollection *fec;
   FiniteElementSpace *fespace;
   FiniteElementSpace *fieldspace;
   fec = new H1_FECollection(order, dim);
   fespace = new FiniteElementSpace(mesh, fec, dim, Ordering::byVDIM);
   fieldspace = new FiniteElementSpace(mesh, fec, 1, Ordering::byVDIM);
   cout << "Number of unknowns: " << fespace->GetVSize() << endl;
  
   int ndof = fespace->GetTrueVSize() / dim;
   cout << "Number of finite element unknowns: " << fespace->GetTrueVSize() << " ndof: " << ndof << " dim: " << dim
        << endl << "Assembling: " << flush;

   //    Define the solution vector x as a finite element grid function
   //    corresponding to fespace. Initialize x with initial guess of zero,
   //    which satisfies the boundary conditions.
   GridFunction x(fespace);
   x = 0.; 

   // Set up the bilinear form a(.,.) on the finite element space
   //    corresponding to the linear elasticity integrator with piece-wise
   //    constants coefficient lambda and mu.
   
   Vector lambda(mesh->attributes.Max());
   lambda = 100.0;
   if (mesh->attributes.Max() > 1)
     lambda(1) = lambda(0)*2;
   PWConstCoefficient lambda_func(lambda); 
   Vector mu(mesh->attributes.Max());
   mu = 75.0;
   if (mesh->attributes.Max() > 1)
     mu(1) = mu(0)*2;
   PWConstCoefficient mu_func(mu);


   // Impose no displacement on the first node
   // which needs to be on x=xmin or x=xmax axis.
   // ux=0, uy=0, uz=0 on this point.
   Array<int> ess_tdof_list;
   ess_tdof_list.SetSize(0);
   {
     GridFunction nodes(fespace);
     mesh->GetNodes(nodes);
     int size = nodes.Size()/dim;
     int found = 0;

     // Traversal of all dofs to detect which one is (0,0,0)
     for (int i = 0; i < size; ++i) {
       double coord[dim]; // coordinates of a node
       double dist = 0.;
       for (int j = 0; j < dim; ++j) {
	 coord[j] = (nodes)[i * dim + j];
	 // because of periodic BC, 0. is also 1.
        if (abs(coord[j] - 1.) < 1e-7)
          coord[j] = 0.;
        dist += coord[j] * coord[j];
       }
       // If distance is close to zero, we have our reference point
       if (dist < 1.e-16) {
	 cout << "coord " << coord[0] << " " << coord[1] << " " << coord[2] << " " << endl;
	 for (int j = 0; j < dim; ++j) {
	   int id_unk;
	   id_unk = i * dim + j;
	   if (id_unk >= 0)
	     {
	       found = 1;
	       ess_tdof_list.Append(id_unk);
	       x(id_unk) = 0.;
	       cout << "blocked unknown: " << id_unk << endl;
	     }
	 }
       }
     }
   }

   for (int i = 0; i < mesh->GetNE(); i++)
   {
      Element *el = mesh->GetElement(i);
      int nv = el->GetNVertices();
      int *v = el->GetVertices();
      for (int j = 0; j < nv; j++) {
	const double *coord = mesh->GetVertex(v[j]);
      }
   }
   BilinearForm *a = new BilinearForm(fespace);
   if (pa) { a->SetAssemblyLevel(AssemblyLevel::PARTIAL); }
   else if (ea) { a->SetAssemblyLevel(AssemblyLevel::ELEMENT); }
   auto ei = new ElasticityIntegrator(lambda_func,mu_func);
   a->AddDomainIntegrator(ei);
   const FiniteElementSpace &fes = *fespace;
   if (!pa && !ea) ei->AssemblePA(*fespace);
   a->Assemble();
  
   // Set up the right-hand side of the FEM linear system.
   LinearForm rhs(fespace);
   rhs.AddDomainIntegrator(new DomainLFGrad2Integrator(lambda_func,mu_func,tcase));
   rhs.Assemble();

   // Assemble the bilinear form and the corresponding linear system,
   // applying any necessary transformations such as: eliminating boundary
   // conditions, applying conforming constraints for non-conforming AMR,
   // static condensation, etc.
   OperatorPtr A;
   Vector B, X;
   a->FormLinearSystem(ess_tdof_list, x, rhs, A, X, B);

   // Define a simple symmetric Gauss-Seidel preconditioner and use it to
   // solve the system Ax=b with PCG.
//   GSSmoother M(A);
//   PCG(A, M, B, X, 1, 500, 1e-24, 0.0);
   CG(*A, B, X, 1, 2100, 1e-24, 0.0);
   //CG(A, B, X, 1, 500, 1e-24, 0.0);

   //  Recover the solution as a finite element grid function.
   a->RecoverFEMSolution(X, rhs, x);
   std::cout << "syst solve nbunkowns= " << B.Size() << " " << x.Size() << " " << X.Size() << "\n" ;
   
   //  Save the results
   if (postproc) {
     ParaViewDataCollection paraview_dc("per", mesh);
     std::string letters = "xyz";
     Array<GridFunction *> stress(dim*(dim+1)/2);
     StressCoefficient stress_c(lambda_func, mu_func);
     stress_c.SetDisplacement(x);
     Array<GridFunction *> strain(dim*(dim+1)/2);
     StrainCoefficient strain_c(lambda_func, mu_func);
     strain_c.SetDisplacement(x);
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
	     strain[c] = new GridFunction(fieldspace);
	     strain_c.SetComponent(si, sj);
	     strain[c]->ProjectDiscCoefficient(strain_c, GridFunction::ARITHMETIC);
	     paraview_dc.RegisterField(strainname,strain[c]);
	     //	     Vector avg, vol;
	     //	     ComputeAvg(*(stress[c]),avg);
	     //	     cout << setprecision(10) << stressname << " " << avg(0) << endl;
	     //	     ComputeAvg(*(strain[c]),avg);
	     //	     cout << setprecision(10) << strainname << " " << avg(0) << endl;
//	     const int nattr = fespace->GetMesh()->attributes.Max();
//	     ComputeAvgPerAtt(*(stress[c]),avg,vol);
//	     for (int attr = 0; attr < nattr; attr++)
//	     {
//	       cout << setprecision(10) << stressname << " attr: " << attr << " avg: "
//		    <<  setw(20) << avg(attr) << ", vol: " << vol(attr) << endl;
//	     }
//	     ComputeAvgPerAtt(*(strain[c]),avg,vol);
//	     for (int attr = 0; attr < nattr; attr++)
//	     {
//	       cout << setprecision(10) << strainname << " attr: " << attr << " avg: "
//		    <<  setw(20) <<  avg(attr) << ", vol: " << vol(attr) << endl;
//	     }
	     const int nattr = fespace->GetMesh()->attributes.Max();
	     Vector vstrain, vstress, vol;
	     ComputeAvgStressStrain(x, lambda_func, mu_func, si, sj, vstrain, vstress, vol);
	     for (int attr = 0; attr < nattr; attr++)
	       {
		 cout << setprecision(10) << strainname << " attr: " << attr << " avg: "
		      <<  setw(20) <<  vstrain(attr) << ", vol: " << vol(attr) << endl;
		 cout << setprecision(10) << stressname << " attr: " << attr << " avg: "
		      <<  setw(20) <<  vstress(attr) << ", vol: " << vol(attr) << endl;
	       }
	     c++;
	   }
       }
     MFEM_VERIFY(c==dim*(dim+1)/2,"");
     paraview_dc.RegisterField("u",&x);
     paraview_dc.SetLevelsOfDetail(order+1);
     paraview_dc.SetDataFormat(VTKFormat::BINARY);
     paraview_dc.SetHighOrderOutput(true);
     paraview_dc.SetCycle(0);
     paraview_dc.SetTime(0.0);
     paraview_dc.Save();
     while (c) { delete stress[--c]; } 
   }

   void (*sol_exact)(const Vector &x, Vector &u);
   switch(tcase) {
   case 1:
     sol_exact=sol_exact1;
     break;
   case 2:
     sol_exact=sol_exact2;
     break;
   case 3:
     sol_exact=sol_exact3;
     break;
   case 4:
     sol_exact=sol_exact4;
     break;
   case 5:
     sol_exact=sol_exact5;
     break;
   case 6:
     sol_exact=sol_exact6;
     break;
   default:
     cout << "Test undefined" << endl;
     exit (1);
     break;
   }
   double errorL2 = 0.;
   if (postproc) {
     VectorFunctionCoefficient sol_coef (dim, sol_exact);
     errorL2= x.ComputeL2Error(sol_coef);
     cerr<<"\ntcase " << tcase << " -- L2 norm: " << errorL2 << endl;
   }
   delete a;
   delete fespace;
   delete fec;
   delete fieldspace;
   delete mesh;
   if (postproc) {
     if (errorL2 < 1e-10)
       {
	 cerr << "OK" << endl;
	 return 0;
       }
     else
       {
	 cerr << "Fail" << endl;
	 return 1;
       }
   }
}


void ComputeAvg(const GridFunction &gf, Vector &avg)
{
   double vol;
   const FiniteElementSpace *fes = gf.FESpace();
   ElementTransformation *T;
   DenseMatrix vals;
   const int dim = fes->GetVDim();
   avg.SetSize(dim);
   avg = 0.;
   for (int i = 0; i < fes->GetNE(); i++)
   {
      const FiniteElement *fe = fes->GetFE(i);
      const IntegrationRule *ir;
      int intorder = 2*fe->GetOrder() + 3; // <----------
      ir = &(IntRules.Get(fe->GetGeomType(), intorder));
      T = fes->GetElementTransformation(i);
      gf.GetVectorValues(*T, *ir, vals);
      MFEM_VERIFY(vals.Width() == ir->GetNPoints(),"");
      MFEM_VERIFY(vals.Height() == dim,"");
      for (int j = 0; j < ir->GetNPoints(); j++)
      {
         const IntegrationPoint &ip = ir->IntPoint(j);
         T->SetIntPoint(&ip);
	 for (int d = 0; d < dim; d++)
         {
	   avg(d) += ip.weight * T->Weight() * vals(d,j);
	 }
	 vol += ip.weight * T->Weight();
      }
   }
   for (int d = 0; d < dim; d++)
   {
     avg(d) /= vol;
   }
}

void ComputeAvgPerAtt(const GridFunction &gf, Vector &avg, Vector &vol)
{
   const FiniteElementSpace *fes = gf.FESpace();
   Vector vals;
   const int dim = fes->GetVDim();
   const int nattr = fes->GetMesh()->attributes.Max();
   MFEM_VERIFY(dim==1,"");
   avg.SetSize(nattr);
   vol.SetSize(nattr);
   avg = 0.;
   vol = 0.;
   for (int i = 0; i < fes->GetNE(); i++)
   {
      const FiniteElement *fe = fes->GetFE(i);
      const int intorder = 2*fe->GetOrder() + 3; // <----------
      const IntegrationRule *ir = &(IntRules.Get(fe->GetGeomType(), intorder));
      ElementTransformation *T = fes->GetElementTransformation(i);
      int attr = fes->GetAttribute(i)-1;
      gf.GetValues(*T, *ir, vals);
      MFEM_VERIFY((attr >= 0) && (attr < nattr),"");
      MFEM_VERIFY(vals.Size() == ir->GetNPoints(),"");
      for (int j = 0; j < ir->GetNPoints(); j++)
      {
         const IntegrationPoint &ip = ir->IntPoint(j);
         T->SetIntPoint(&ip);
	 avg(attr) += ip.weight * T->Weight() * vals(j);
	 vol(attr) += ip.weight * T->Weight();
      }
   }
   for (int attr = 0; attr < nattr; attr++)
   {
     avg(attr) /= vol(attr);
   }
}

void ComputeAvgStressStrain(const GridFunction &u, Coefficient &lambda, Coefficient &mu,
			    int si, int sj, Vector &strain, Vector &stress, Vector &vol)
{
   const FiniteElementSpace *fes = u.FESpace();
   Vector vals;
   DenseMatrix grad;
   double fstrain, fstress;
   const int dim = fes->GetVDim();
   const int nattr = fes->GetMesh()->attributes.Max();

   stress.SetSize(nattr);
   strain.SetSize(nattr);
   vol.SetSize(nattr);
   stress = 0.;
   strain = 0.;
   vol = 0.;
   for (int i = 0; i < fes->GetNE(); i++)
   {
      const FiniteElement *fe = fes->GetFE(i);
      const int intorder = 2*fe->GetOrder() + 3; // <----------
      const IntegrationRule *ir = &(IntRules.Get(fe->GetGeomType(), intorder));
      ElementTransformation *T = fes->GetElementTransformation(i);
      int attr = fes->GetAttribute(i)-1;
      //      u.GetValues(*T, *ir, vals);
      MFEM_VERIFY((attr >= 0) && (attr < nattr),"");
      for (int j = 0; j < ir->GetNPoints(); j++)
      {
        const IntegrationPoint &ip = ir->IntPoint(j);
	const double L = lambda.Eval(*T, ip);
	const double M = mu.Eval(*T, ip);
        T->SetIntPoint(&ip);
        u.GetVectorGradient(*T, grad);
        if (si == sj)
        {
          double div_u = grad.Trace();
          fstress = L*div_u + 2*M*grad(si,si);
	  fstrain = grad(si,si);
        }
        else
        {
          fstress = M*(grad(si,sj) + grad(sj,si));
	  fstrain = .5*(grad(si,sj) + grad(sj,si));
        }
        strain(attr) += ip.weight * T->Weight() * fstrain;
        stress(attr) += ip.weight * T->Weight() * fstress;
        vol(attr)    += ip.weight * T->Weight();
      }
   }
   for (int attr = 0; attr < nattr; attr++)
   {
     strain(attr) /= vol(attr);
     stress(attr) /= vol(attr);
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
      return L*div_u + 2*M*grad(si,si);
   }
   else
   {
      return M*(grad(si,sj) + grad(sj,si));
   }
}


double StrainCoefficient::Eval(ElementTransformation &T,
                               const IntegrationPoint &ip)
{
   MFEM_ASSERT(u != NULL, "displacement field is not set");

   u->GetVectorGradient(T, grad);
   return (grad(si,sj)+grad(sj,si));
}

void sol_exact1(const Vector &x, Vector &u)
{
  const double gradx = 1./3.;
  u =  0.;
  if (x(0) < xthr)
    {
      u(0) = gradx*x(0);
    }
  else
    {
      u(0) = gradx*xthr - gradx*(x(0)-xthr);
    }
}

void sol_exact2(const Vector &x, Vector &u)
{
  const double gradx = 4./30.;
  u =  0.;
  if (x(0) < xthr)
    {
      u(0) = gradx*x(0);
    }
  else
    {
      u(0) = gradx*xthr - gradx*(x(0)-xthr);
    }
}

void sol_exact3(const Vector &x, Vector &u)
{
  const double gradx = 4./30.;
  u =  0.;
  if (x(0) < xthr)
    {
      u(0) = gradx*x(0);
    }
  else
    {
      u(0) = gradx*xthr - gradx*(x(0)-xthr);
    }
}

void sol_exact4(const Vector &x, Vector &u)
{
  const double gradx = 1./3.;
  u =  0.;
  if (x(0) < xthr)
    {
      u(1) = gradx*x(0);
    }
  else
    {
      u(1) = gradx*xthr - gradx*(x(0)-xthr);
    }
}

void sol_exact5(const Vector &x, Vector &u)
{
  u =  0.;
}

void sol_exact6(const Vector &x, Vector &u)
{
  const double gradx = 1./3.;
  u =  0.;
  if (x(0) < xthr)
    {
      u(2) = gradx*x(0);
    }
  else
    {
      u(2) = gradx*xthr - gradx*(x(0)-xthr);
    }
}
