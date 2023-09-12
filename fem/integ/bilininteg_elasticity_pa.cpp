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

   dim = el.GetDim();
   dof = el.GetDof();
   nq = ir->GetNPoints();
   MFEM_VERIFY(dim == 2 || dim == 3, "");

   ne = fes.GetNE();
   geom = mesh->GetGeometricFactors(*ir, GeometricFactors::JACOBIANS, mt);

   pa_data.SetSize((2+dim*dim) * nq * ne, mt);
   QuadratureSpace qs(*mesh, *ir);

   auto LM = Reshape(pa_data.Write(), 2+dim*dim, nq, ne);
   mfem::forall(ne, [=] MFEM_HOST_DEVICE (int e)
     {
       for (int i = 0; i < ir -> GetNPoints(); i++) {
	 ElementTransformation *Trans = fes.GetElementTransformation(e);
	 const IntegrationPoint &ip = ir->IntPoint(i);
	 Trans->SetIntPoint(&ip);
	 double w = ip.weight * Trans->Weight();
	 LM(0,i,e) = w * mu->Eval(*Trans, ip);
	 LM(1,i,e) = w * lambda->Eval(*Trans, ip);
	 const DenseMatrix &invJ = Trans->InverseJacobian();
	 for (int d1=0; d1<dim; d1++)
	   for (int d2=0; d2<dim; d2++) {
	     LM(2+d2+d1*dim,i,e) = invJ(d2,d1);
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
  auto X = Reshape(x_.Read(), dof, dim, ne);
  auto Y = Reshape(y_.ReadWrite(), dof, dim, ne);
  std::cout << "nq : " << nq << "\n";
  mfem::forall(ne, [=] MFEM_HOST_DEVICE (int e)
  {
    for (int i = 0; i < nq; i++)
      {
	DenseMatrix gshape(dof, dim);

	const double LW = LM(1,i,e);
	const double MW = LM(0,i,e);
	
	for (int kk = 0; kk < dim; kk++) 
	  for (int ll = 0; ll < dof; ll++) {
	    register double accu = 0.;
	    for (int ii = 0; ii < dim; ii++)
	      accu += maps->Gt[ll+dof*(i+nq*ii)] * LM(2+ii+kk*dim,i,e);
	    gshape(ll,kk) = accu;
	  } 
	
	for (int ii = 0; ii < dim; ii++)
	  for (int kk = 0; kk < dof; kk++) {
	    for (int jj = 0; jj < dim; jj++) {
	      double contribB = 0.;
	      double contribA = 0;
	      double contribC = 0;
	      for (int ll = 0; ll < dof; ll++) {
		contribA += X(ll,ii,e) * gshape(ll, jj);
		contribB += X(ll,jj,e) * gshape(ll, jj);
		contribC += X(ll,jj,e) * gshape(ll, ii); 
	      }
	      Y(kk,ii,e) += MW * gshape(kk, jj) * (contribA + contribC) + LW * gshape(kk, ii) * contribB ;
	    }
	  }
      }  
  });
}

} // namespace mfem

