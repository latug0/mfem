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

namespace mfem
{

const int MAXNDOF = ElasticityIntegrator::MAXNDOF;
const int MAXNQ = ElasticityIntegrator::MAXNQ;

template<int T_D1D = 0, int T_Q1D = 0>
static void EAElasticityAssemble2D(const int NE,
				   const int NDOFS,
				   const int NQ,
				   const Array<double> &pGt,
				   const Vector &padata,
				   Vector &eadata,
				   const bool add,
				   const int d1d = 0,
				   const int q1d = 0)
{
   constexpr int MDIM = 2;
   auto Gt = Reshape(pGt.Read(),NQ*NDOFS*MDIM);
   auto LM = Reshape(padata.Read(), 2+MDIM*MDIM, NQ, NE);

   // TODO: Check ndof value
   auto A = Reshape(eadata.ReadWrite(), NDOFS, NDOFS, NE);
   mfem::forall(NE, [=] MFEM_HOST_DEVICE (int e)
   {
     for (int i=0; i<NQ; i++) {
       double gshape[MAXNDOF][MDIM];
       double LW = LM(1,i,e);
       double MW = LM(0,i,e);
       
       for (int kk = 0; kk < MDIM; kk++) 
	 for (int ll = 0; ll < NDOFS; ll++) {
	   gshape[ll][kk] = 0.;
	   for (int mm = 0; mm < MDIM; mm++)
	     gshape[ll][kk] +=
	       Gt[ll+NDOFS*(i+NQ*mm)] * LM(2+mm+kk*MDIM,i,e);
	 }
       
       for (int ii = 0; ii < MDIM; ii++)
	 for (int jj = 0; jj < MDIM; jj++) {
	   for (int kk = 0; kk < NDOFS; kk++) {
	     for (int ll = 0; ll < NDOFS; ll++) {
//	       A(kk,ii,ll,ii,e) +=
//		 MW * gshape[kk][jj] * gshape[ll][jj];
	     }
	   }
	 }

       for (int ii = 0; ii < MDIM; ii++)
	 for (int jj = 0; jj < MDIM; jj++) {
	   for (int kk = 0; kk < NDOFS; kk++) {
	     for (int ll = 0; ll < NDOFS; ll++) {
//	       A(kk,ii,ll,jj,e) +=
//		 LW * gshape[kk][ii] * gshape[ll][jj] +
//		 MW * gshape[kk][jj] * gshape[ll][ii];
	     }
	   }
	 }
     }
   });
}

  
template<int T_D1D = 0, int T_Q1D = 0>
static void EAElasticityAssemble3D(const int NE,
				   const int NDOFS,
				   const int NQ,
				   const Array<double> &pGt,
				   const Vector &padata,
				   Vector &eadata,
				   const bool add,
				   const int d1d = 0,
				   const int q1d = 0)
{
   constexpr int MDIM = 3;
   auto Gt = Reshape(pGt.Read(),NQ*NDOFS*MDIM);
   auto LM = Reshape(padata.Read(), 2+MDIM*MDIM, NQ, NE);

   // TODO: Check ndof value
   auto A = Reshape(eadata.ReadWrite(), NDOFS, NDOFS, NE);
   mfem::forall(NE, [=] MFEM_HOST_DEVICE (int e)
   {
     for (int i=0; i<NQ; i++) {
       double gshape[MAXNDOF][MDIM];
       double LW = LM(1,i,e);
       double MW = LM(0,i,e);
       
       for (int kk = 0; kk < MDIM; kk++) 
	 for (int ll = 0; ll < NDOFS; ll++) {
	   gshape[ll][kk] = 0.;
	   for (int mm = 0; mm < MDIM; mm++)
	     gshape[ll][kk] +=
	       Gt[ll+NDOFS*(i+NQ*mm)] * LM(2+mm+kk*MDIM,i,e);
	 }
       
       for (int ii = 0; ii < MDIM; ii++)
	 for (int jj = 0; jj < MDIM; jj++) {
	   for (int kk = 0; kk < NDOFS; kk++) {
	     for (int ll = 0; ll < NDOFS; ll++) {
//	       A(kk,ii,ll,ii,e) +=
//		 MW * gshape[kk][jj] * gshape[ll][jj];
	     }
	   }
	 }

       for (int ii = 0; ii < MDIM; ii++)
	 for (int jj = 0; jj < MDIM; jj++) {
	   for (int kk = 0; kk < NDOFS; kk++) {
	     for (int ll = 0; ll < NDOFS; ll++) {
//	       A(kk,ii,ll,jj,e) +=
//		 LW * gshape[kk][ii] * gshape[ll][jj] +
//		 MW * gshape[kk][jj] * gshape[ll][ii];
	     }
	   }
	 }
     }
   });
}
  
  
void ElasticityIntegrator::AssembleEA(const FiniteElementSpace &fes,
				      Vector &ea_data,
				      const bool add)
{
   AssemblePA(fes);
   //   ne = fes.GetMesh()->GetNE();
   const Array<double> &Gt = maps->Gt;

   std::cout << "param_ea predict ne*dim= " << ne*dim <<  "\n" ;
   std::cout << "param_ea ndof= " << ndof  << "  nq=" << nq  << " ne=" << ne << "\n" ;
   std::cout << "          d1d=" << dofs1D << " q1d=" << quad1D << "\n" ;

   if (dim == 2)
   {
     return EAElasticityAssemble2D(ne,ndof,nq,Gt,pa_data,ea_data,
				   add,dofs1D,quad1D);
   }
   else if (dim == 3)
   {
     return EAElasticityAssemble3D(ne,ndof,nq,Gt,pa_data,ea_data,
				   add,dofs1D,quad1D);
   }
   MFEM_ABORT("Unknown kernel.");
}

  

} // namespace mfem

