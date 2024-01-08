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
static void EADiffusionAssemble2D(const int NE,
                                  const Array<double> &b,
                                  const Array<double> &g,
                                  const Vector &padata,
                                  Vector &eadata,
                                  const bool add,
                                  const int d1d = 0,
                                  const int q1d = 0)
{
   const int D1D = T_D1D ? T_D1D : d1d;
   const int Q1D = T_Q1D ? T_Q1D : q1d;
   MFEM_VERIFY(D1D <= MAX_D1D, "");
   MFEM_VERIFY(Q1D <= MAX_Q1D, "");
   auto B = Reshape(b.Read(), Q1D, D1D);
   auto G = Reshape(g.Read(), Q1D, D1D);
   auto D = Reshape(padata.Read(), Q1D, Q1D, 3, NE);
   auto A = Reshape(eadata.ReadWrite(), D1D, D1D, D1D, D1D, NE);
   mfem::forall_2D(NE, D1D, D1D, [=] MFEM_HOST_DEVICE (int e)
   {
      const int D1D = T_D1D ? T_D1D : d1d;
      const int Q1D = T_Q1D ? T_Q1D : q1d;
      constexpr int MD1 = T_D1D ? T_D1D : MAX_D1D;
      constexpr int MQ1 = T_Q1D ? T_Q1D : MAX_Q1D;
      double r_B[MQ1][MD1];
      double r_G[MQ1][MD1];
   });
}

  
template<int T_D1D = 0, int T_Q1D = 0>
static void EADiffusionAssemble3D(const int NE,
                                  const Array<double> &b,
                                  const Array<double> &g,
                                  const Vector &padata,
                                  Vector &eadata,
                                  const bool add,
                                  const int d1d = 0,
                                  const int q1d = 0)
{
   const int D1D = T_D1D ? T_D1D : d1d;
   const int Q1D = T_Q1D ? T_Q1D : q1d;
   MFEM_VERIFY(D1D <= MAX_D1D, "");
   MFEM_VERIFY(Q1D <= MAX_Q1D, "");
   auto B = Reshape(b.Read(), Q1D, D1D);
   auto G = Reshape(g.Read(), Q1D, D1D);
   auto D = Reshape(padata.Read(), Q1D, Q1D, Q1D, 6, NE);
   auto A = Reshape(eadata.ReadWrite(), D1D, D1D, D1D, D1D, D1D, D1D, NE);
   mfem::forall_3D(NE, D1D, D1D, D1D, [=] MFEM_HOST_DEVICE (int e)
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
   ne = fes.GetMesh()->GetNE();
   const Array<double> &B = maps->B;
   const Array<double> &G = maps->G;
   if (dim == 2)
   {
      return EADiffusionAssemble2D(ne,B,G,pa_data,ea_data,add,
				   dofs1D,quad1D);
   }
   else if (dim == 3)
   {
     return EADiffusionAssemble3D(ne,B,G,pa_data,ea_data,add,
				  dofs1D,quad1D);
   }
   MFEM_ABORT("Unknown kernel.");
}

  

} // namespace mfem

