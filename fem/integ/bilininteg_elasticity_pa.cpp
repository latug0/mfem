// Copyright (c) 2010-2023, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#include "../../general/forall.hpp"
#include "../bilininteg.hpp"
#include "../gridfunc.hpp"
#include "../qfunction.hpp"
#include "../fespace.hpp"

namespace mfem
{

  constexpr int MAXNDOF = 45;
  constexpr int MAXNQ = 44;
  static void PAElasticitySetup_bis(const int dim,
				    const int nq,
				    const int ne,
				    const int ndof,
				    const Vector &pJ,
				    const Array<double> &pGt,
				    const Array<double> &pW,
				    const Vector &pMU,
				    const Vector &pLAMBDA,
				    Vector &op) {

    auto J = Reshape(pJ.Read(), nq, dim, dim, ne);
    auto W = pW.Read();
    auto MU = Reshape(pMU.Read(),nq,ne);
    auto LAMBDA = Reshape(pLAMBDA.Read(),nq,ne);
    auto LM = Reshape(op.Write(), 2+dim*dim, nq, ne);
    mfem::forall(ne, [=] MFEM_HOST_DEVICE (int e)
    {
      for (int i = 0; i < nq; i++) {
	if (dim == 3) {
	  const double J11 = J(i,0,0,e);
	  const double J21 = J(i,1,0,e);
	  const double J31 = J(i,2,0,e);
	  const double J12 = J(i,0,1,e);
	  const double J22 = J(i,1,1,e);
	  const double J32 = J(i,2,1,e);
	  const double J13 = J(i,0,2,e);
	  const double J23 = J(i,1,2,e);
	  const double J33 = J(i,2,2,e);
	  const double detJ =
	    J11 * (J22 * J33 - J32 * J23) -
	     J21 * (J12 * J33 - J32 * J13) +
	    J31 * (J12 * J23 - J22 * J13);
	  const double i_detJ = 1/detJ;
	  LM(0,i,e) = W[i] * i_detJ * MU(i,e);
	  LM(1,i,e) = W[i] * i_detJ * LAMBDA(i,e);
	  // adj(J)
	  const double A11 = (J22 * J33) - (J23 * J32);
	  const double A12 = (J32 * J13) - (J12 * J33);
	  const double A13 = (J12 * J23) - (J22 * J13);
	  const double A21 = (J31 * J23) - (J21 * J33);
	  const double A22 = (J11 * J33) - (J13 * J31);
	  const double A23 = (J21 * J13) - (J11 * J23);
	  const double A31 = (J21 * J32) - (J31 * J22);
	  const double A32 = (J31 * J12) - (J11 * J32);
	  const double A33 = (J11 * J22) - (J12 * J21);
	  LM( 2,i,e) = A11; // 1,1
	  LM( 3,i,e) = A21; // 1,2
	  LM( 4,i,e) = A31; // 1,3
	  LM( 5,i,e) = A12; // 2,1
	  LM( 6,i,e) = A22; // 2,2
	  LM( 7,i,e) = A32; // 2,3
	  LM( 8,i,e) = A13; // 3,1
	  LM( 9,i,e) = A23; // 3,2
	  LM(10,i,e) = A33; // 3,3
	} else if (dim == 2) {
	  const double J11 = J(i,0,0,e);
	  const double J21 = J(i,1,0,e);
	  const double J12 = J(i,0,1,e);
	  const double J22 = J(i,1,1,e);
	  const double detJ = ((J11*J22)-(J21*J12));
	  const double i_detJ = 1/detJ;
	  LM(0,i,e) = W[i] * i_detJ * MU(i,e);
	  LM(1,i,e) = W[i] * i_detJ * LAMBDA(i,e);
	  const double iJ11    =  J22; 
	  const double iJ12    = -J12; 
	  const double iJ21    = -J21; 
	  const double iJ22    =  J11;
	  LM(2,i,e) = iJ11;
	  LM(3,i,e) = iJ21;
	  LM(4,i,e) = iJ12;
	  LM(5,i,e) = iJ22;
	}
      }
    });
  }

  static void PAElasticitySetup(const int dim,
				const int nq,
				const int ne,
				const int ndof,
				const Vector &pJ,
				const Array<double> &pGt,
				const Array<double> &pW,
				const Vector &pMU,
				const Vector &pLAMBDA,
				Vector &op) {

    if (dim == 1) { MFEM_ABORT("dim==1 not supported in PAElasticitySetup"); }
    if (dim == 2)
      {
	PAElasticitySetup_bis(dim, nq, ne, ndof, pJ, pGt, pW, pMU, pLAMBDA, op);
      }
    if (dim == 3)
      {
	PAElasticitySetup_bis(dim, nq, ne, ndof, pJ, pGt, pW, pMU, pLAMBDA, op);
      }
  }

void ElasticityIntegrator::AssemblePA(const FiniteElementSpace &fes)
{
   const MemoryType mt = (pa_mt == MemoryType::DEFAULT) ?
                         Device::GetDeviceMemoryType() : pa_mt;
   Mesh *mesh = fes.GetMesh();
   const FiniteElement &el = *fes.GetFE(0);
   ElementTransformation &TransRef = *fes.GetElementTransformation(0);
   const IntegrationRule *ir = IntRule ? IntRule :  
     &IntRules.Get(el.GetGeomType(), 2 * TransRef.OrderGrad(&el));
   maps = &el.GetDofToQuad(*ir, DofToQuad::FULL);
   auto  geom = mesh->GetGeometricFactors(*ir, GeometricFactors::JACOBIANS, mt);

   QuadratureSpace qs(*mesh, *ir);
   CoefficientVector coeffmu(qs, CoefficientStorage::COMPRESSED);
   CoefficientVector coefflambda(qs, CoefficientStorage::COMPRESSED);
   coeffmu.Project(*mu);
   coefflambda.Project(*lambda);
   dim = el.GetDim();
   ndof = el.GetDof();
   nq = ir->GetNPoints();
   ne = fes.GetNE();
   MFEM_VERIFY(dim == 2 || dim == 3, "");
   pa_data.SetSize((2+dim*dim) * nq * ne, mt);
   if (ndof > MAXNDOF)
     MFEM_ABORT("MAXNDOF too low");
   if (nq > MAXNQ)
     MFEM_ABORT("MAXNQ too low");
   
   std::cout << "optim ne=" << ne << " nq=" << nq << " ndof " << ndof << "\n" ;
  
   PAElasticitySetup(dim, nq, ne, ndof, geom->J, maps->Gt, ir->GetWeights(),
		     coeffmu, coefflambda, pa_data);
}

void  ElasticityIntegrator::AssembleDiagonalPA(Vector& diag)
{
  MFEM_ABORT("ElasticityIntegrator::AssembleDiagonalPA not defined");
}


void PAElasticityApply2D(const int dim,
			 const int nq,
			 const int ne,
			 const int ndof,
			 const Array<double> &pGt,
			 const Vector &op,
			 const Vector &px,
			 Vector &py)
{
  constexpr int MDIM = 2;
  auto LM = Reshape(op.Read(), 2+dim*dim, nq, ne);
  auto Gt = Reshape(pGt.Read(),nq*ndof*dim);
  auto X = Reshape(px.Read(), ndof, dim, ne);
  auto Y = Reshape(py.ReadWrite(), ndof, dim, ne);
  
  mfem::forall(ne, [=] MFEM_HOST_DEVICE (int e)
      {
//	int e = q_global / nq;
//	int i = q_global % nq;
	for (int i=0; i<nq; i++) {
	  double gshape[MAXNDOF][MDIM];

	  double LW = LM(1,i,e);
	  double MW = LM(0,i,e);
	  
	for (int kk = 0; kk < MDIM; kk++) 
	  for (int ll = 0; ll < ndof; ll++) {
	    gshape[ll][kk] = 0.;
	    for (int ii = 0; ii < MDIM; ii++)
	      gshape[ll][kk] += Gt[ll+ndof*(i+nq*ii)] * LM(2+ii+kk*MDIM,i,e);
	  }
	
	for (int ii = 0; ii < MDIM; ii++)
	  for (int jj = 0; jj < MDIM; jj++) {
	    double contribB = 0.;
	    double contribA = 0;
	    double contribC = 0;
	    for (int ll = 0; ll < ndof; ll++) {
	      contribA += X(ll,ii,e) * gshape[ll][jj];
	      contribB += X(ll,jj,e) * gshape[ll][jj];
	      contribC += X(ll,jj,e) * gshape[ll][ii]; 
	    }
	    for (int kk = 0; kk < ndof; kk++) {
	      // conflict in writing if parallelized on dofs
	      Y(kk,ii,e) += MW * gshape[kk][jj] * (contribA + contribC) +
		LW * gshape[kk][ii] * contribB ;
	    }
	  }
	}
      });
}

void PAElasticityApply2D_optim(const int dim,
			       const int nq,
			       const int ne,
			       const int ndof,
			       const Array<double> &pGt,
			       const Vector &op,
			       const Vector &px,
			       Vector &py)
{
  constexpr int MDIM = 2;
  auto LM = Reshape(op.Read(), 2+dim*dim, nq, ne);
  auto Gt = Reshape(pGt.Read(),nq*ndof*dim);
  auto X = Reshape(px.Read(), ndof, dim, ne);
  auto Y = Reshape(py.ReadWrite(), ndof, dim, ne);
  int maxndofnq = (ndof>nq?ndof:nq);
  mfem::forall_2D(ne, maxndofnq, MDIM, [=] MFEM_HOST_DEVICE (int e)
      {
	const int tidx = MFEM_THREAD_ID(x);
	const int tidy = MFEM_THREAD_ID(y);
	MFEM_SHARED double gshape[MAXNQ*MAXNDOF][MDIM];
	MFEM_FOREACH_THREAD(ll,x,ndof) {  
	  MFEM_FOREACH_THREAD(jj,y,MDIM) {
	    for (int i=0; i<nq; i++) {
	      gshape[ll+ndof*i][jj] = 0.;
	      for (int ii = 0; ii < MDIM; ii++)
		gshape[ll+ndof*i][jj] += Gt[ll+ndof*(i+nq*ii)] * LM(2+ii+jj*MDIM,i,e);
	    }
	  }
	}
	MFEM_SYNC_THREAD;
 
	MFEM_SHARED double contribA[MAXNQ][MDIM*MDIM];
	MFEM_SHARED double contribB[MAXNQ][MDIM*MDIM];
	MFEM_SHARED double contribC[MAXNQ][MDIM*MDIM];
	MFEM_FOREACH_THREAD(jj,y,MDIM) {
	  MFEM_FOREACH_THREAD(i,x,nq) {
	    for (int ii = 0; ii < MDIM; ii++) {
	      //	      for (int i = 0; i<nq; i++) {
	      contribB[i][jj+MDIM*ii] = 0.;
	      contribA[i][jj+MDIM*ii] = 0;
	      contribC[i][jj+MDIM*ii] = 0;
	      for (int ll = 0; ll < ndof; ll++) {
		contribA[i][jj+MDIM*ii] += X(ll,ii,e) * gshape[ll+ndof*i][jj];
		contribB[i][jj+MDIM*ii] += X(ll,jj,e) * gshape[ll+ndof*i][jj];
		contribC[i][jj+MDIM*ii] += X(ll,jj,e) * gshape[ll+ndof*i][ii]; 
	      }
	    }
	  }
	}
	MFEM_SYNC_THREAD;
	
	MFEM_FOREACH_THREAD(kk,x,ndof) {
	  MFEM_FOREACH_THREAD(ii,y,MDIM) {
	    for (int jj = 0; jj < MDIM; jj++) {
	      for (int i = 0; i < nq; i++) {	    
		const double LW = LM(1,i,e);
		const double MW = LM(0,i,e);
		Y(kk,ii,e) +=
		  MW * gshape[kk+ndof*i][jj] * (contribA[i][jj+MDIM*ii]+ contribC[i][jj+MDIM*ii]) +
		  LW * gshape[kk+ndof*i][ii] * contribB[i][jj+MDIM*ii] ;
	      }
	    }
	  }
	}
	MFEM_SYNC_THREAD;
      });
}

void PAElasticityApply3D(const int dim,
			 const int nq,
			 const int ne,
			 const int ndof,
			 const Array<double> &pGt,
			 const Vector &op,
			 const Vector &px,
			 Vector &py)
{
  constexpr int MDIM = 3;
  auto LM = Reshape(op.Read(), 2+dim*dim, nq, ne);
  auto Gt = Reshape(pGt.Read(),nq*ndof*dim);
  auto X = Reshape(px.Read(), ndof, dim, ne);
  auto Y = Reshape(py.ReadWrite(), ndof, dim, ne);
  mfem::forall(ne, [=] MFEM_HOST_DEVICE (int e)
      {
	for (int i=0; i<nq; i++) {
	double gshape[MAXNDOF][MDIM];
	double LW = LM(1,i,e);
	double MW = LM(0,i,e);
	
	for (int kk = 0; kk < MDIM; kk++) 
	  for (int ll = 0; ll < ndof; ll++) {
	    gshape[ll][kk] = 0.;
	    for (int mm = 0; mm < MDIM; mm++)
	      gshape[ll][kk] += Gt[ll+ndof*(i+nq*mm)] * LM(2+mm+kk*MDIM,i,e);
	  }
	
	for (int ii = 0; ii < MDIM; ii++)
	  for (int jj = 0; jj < MDIM; jj++) {
	    double contribB = 0.;
	    double contribA = 0;
	    double contribC = 0;
	    for (int ll = 0; ll < ndof; ll++) {
	      contribA += X(ll,ii,e) * gshape[ll][jj];
	      contribB += X(ll,jj,e) * gshape[ll][jj];
	      contribC += X(ll,jj,e) * gshape[ll][ii]; 
	    }
	    for (int kk = 0; kk < ndof; kk++) {
	      Y(kk,ii,e) += MW * gshape[kk][jj] * (contribA + contribC) +
		LW * gshape[kk][ii] * contribB ;
	    }
	  }
      }
      });
}

void ElasticityIntegrator::AddMultPA(const Vector &x, Vector &y) const
{
    if (dim == 1) { MFEM_ABORT("dim==1 not supported in PAElasticitySetup"); }
    if (dim == 2)
      {
	PAElasticityApply2D_optim(dim, nq, ne, ndof, maps->Gt, pa_data, x, y);
      }
    if (dim == 3)
      {
	PAElasticityApply3D(dim, nq, ne, ndof, maps->Gt, pa_data, x, y);
      }
}

} // namespace mfem

