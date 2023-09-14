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

void ElasticityIntegrator::AssemblePA(const FiniteElementSpace &fes)
{
   const MemoryType mt = (pa_mt == MemoryType::DEFAULT) ?
                         Device::GetDeviceMemoryType() : pa_mt;
   Mesh *mesh = fes.GetMesh();
   const FiniteElement &el = *fes.GetFE(0);
   ElementTransformation &TransRef = *fes.GetElementTransformation(0);
   const IntegrationRule *ir = IntRule ? IntRule :  
     &IntRules.Get(el.GetGeomType(), 2 * TransRef.OrderGrad(&el));
   if (geom == NULL) 
     geom = mesh->GetGeometricFactors(*ir, GeometricFactors::JACOBIANS);
   if (maps == NULL) 
     maps = &el.GetDofToQuad(*ir, DofToQuad::FULL);

   QuadratureSpace qs(*mesh, *ir);
   CoefficientVector coeffmu(qs, CoefficientStorage::COMPRESSED);
   CoefficientVector coefflambda(qs, CoefficientStorage::COMPRESSED);
   coeffmu.Project(*mu);
   coefflambda.Project(*lambda);
   dim = el.GetDim();
   dof = el.GetDof();
   nq = ir->GetNPoints();
   MFEM_VERIFY(dim == 2 || dim == 3, "");

   ne = fes.GetNE();
   geom = mesh->GetGeometricFactors(*ir, GeometricFactors::JACOBIANS, mt);

   pa_data.SetSize((2+dim*dim) * nq * ne, mt);


   auto J = Reshape(geom->J.Read(), nq, dim, dim, ne);
   auto MU = Reshape(coeffmu.Read(),nq,ne);
   auto LAMBDA = Reshape(coefflambda.Read(),nq,ne);
   auto LM = Reshape(pa_data.Write(), 2+dim*dim, nq, ne);
   auto W = ir->GetWeights().Read();
   mfem::forall(ne, [=] MFEM_HOST_DEVICE (int e)
     {
       for (int i = 0; i < ir -> GetNPoints(); i++) {
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

void  ElasticityIntegrator::AssembleDiagonalPA(Vector& diag)
{
  MFEM_ABORT("ElasticityIntegrator::AssembleDiagonalPA not defined");
}


void ElasticityIntegrator::AddMultPA(const Vector &x_, Vector &y_) const
{
  auto LM = Reshape(pa_data.Read(), 2+dim*dim, nq, ne);
  auto Gt = maps->Gt.Read();
  auto X = Reshape(x_.Read(), dof, dim, ne);
  auto Y = Reshape(y_.ReadWrite(), dof, dim, ne);

  mfem::forall(ne*nq, [=] MFEM_HOST_DEVICE (int q_global)
      {
	const int e = q_global / nq;
	const int i = q_global % nq;
	double gshape[dof][dim];

	const double LW = LM(1,i,e);
	const double MW = LM(0,i,e);
	
	for (int kk = 0; kk < dim; kk++) 
	  for (int ll = 0; ll < dof; ll++) {
	    gshape[ll][kk] = 0.;
	    for (int ii = 0; ii < dim; ii++)
	      gshape[ll][kk] += Gt[ll+dof*(i+nq*ii)] * LM(2+ii+kk*dim,i,e);
	  }
	
	for (int ii = 0; ii < dim; ii++)
	  for (int jj = 0; jj < dim; jj++) {
	    double contribB = 0.;
	    double contribA = 0;
	    double contribC = 0;
	    for (int ll = 0; ll < dof; ll++) {
	      contribA += X(ll,ii,e) * gshape[ll][jj];
	      contribB += X(ll,jj,e) * gshape[ll][jj];
	      contribC += X(ll,jj,e) * gshape[ll][ii]; 
	    }
	    for (int kk = 0; kk < dof; kk++) {
	      Y(kk,ii,e) += MW * gshape[kk][jj] * (contribA + contribC) +
		LW * gshape[kk][ii] * contribB ;
	    }
	  }
      });
}

} // namespace mfem

