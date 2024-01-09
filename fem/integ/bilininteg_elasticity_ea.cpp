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


template<int T_D1D = 0, int T_Q1D = 0>
static void EAElasticityAssemble2D(const int NE,
				   const int NDOFS,
				   const int NQ,
				   const Array<double> &b,
				   const Array<double> &g,
				   const Vector &padata,
				   Vector &eadata,
				   const bool add,
				   const int d1d = 0,
				   const int q1d = 0)
{
   constexpr int MDIM = 2;
   auto LM = Reshape(padata.Read(), 2+MDIM*MDIM, NQ, NE);

   // TODO: Check ndof value
   auto A = Reshape(eadata.ReadWrite(), NDOFS, NDOFS, NE);
   mfem::forall(NE*NDOFS, [=] MFEM_HOST_DEVICE (int e)
   {
   });
}

  
template<int T_D1D = 0, int T_Q1D = 0>
static void EAElasticityAssemble3D(const int NE,
				   const int NDOFS,
				   const int NQ,
				   const Array<double> &b,
				   const Array<double> &g,
				   const Vector &padata,
				   Vector &eadata,
				   const bool add,
				   const int d1d = 0,
				   const int q1d = 0)
{
   constexpr int MDIM = 3;
   auto LM = Reshape(padata.Read(), 2+MDIM*MDIM, NQ, NE);
   auto A = Reshape(eadata.ReadWrite(), NDOFS, NDOFS, NE);
   mfem::forall(NE*NDOFS, [=] MFEM_HOST_DEVICE (int e)
   {
      const int D1D = T_D1D ? T_D1D : d1d;
      const int Q1D = T_Q1D ? T_Q1D : q1d;
      constexpr int MD1 = T_D1D ? T_D1D : MAX_D1D;
      constexpr int MQ1 = T_Q1D ? T_Q1D : MAX_Q1D;
   });
}
  
  
void ElasticityIntegrator::AssembleEA(const FiniteElementSpace &fes,
				      Vector &ea_data,
				      const bool add)
{
   AssemblePA(fes);
   //   ne = fes.GetMesh()->GetNE();
   const Array<double> &B = maps->B;
   const Array<double> &G = maps->G;
   if (dim == 2)
   {
     return EAElasticityAssemble2D(ne,ndof,nq,B,G,pa_data,ea_data,
				   add,dofs1D,quad1D);
   }
   else if (dim == 3)
   {
     return EAElasticityAssemble3D(ne,ndof,nq,B,G,pa_data,ea_data,
				   add,dofs1D,quad1D);
   }
   MFEM_ABORT("Unknown kernel.");
}

  

} // namespace mfem

