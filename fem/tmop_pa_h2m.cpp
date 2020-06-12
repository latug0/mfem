// Copyright (c) 2010-2020, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#include "tmop.hpp"
#include "tmop_pa.hpp"
#include "linearform.hpp"
#include "pgridfunc.hpp"
#include "tmop_tools.hpp"
#include "../general/forall.hpp"
#include "../linalg/kernels.hpp"
#include "../linalg/dtensor.hpp"

namespace mfem
{

template<int T_D1D = 0, int T_Q1D = 0, int T_NBZ = 0, int T_MAX = 0>
static void AddMultGradPA_Kernel_2D(const int NE,
                                    const Array<double> &b_,
                                    const Array<double> &g_,
                                    const DenseTensor &j_,
                                    const Vector &p_,
                                    const Vector &x_,
                                    Vector &y_,
                                    const int d1d = 0,
                                    const int q1d = 0)
{
   constexpr int DIM = 2;
   constexpr int NBZ = T_NBZ ? T_NBZ : 1;

   const int D1D = T_D1D ? T_D1D : d1d;
   const int Q1D = T_Q1D ? T_Q1D : q1d;

   const auto b = Reshape(b_.Read(), Q1D, D1D);
   const auto g = Reshape(g_.Read(), Q1D, D1D);
   const auto J = Reshape(j_.Read(), DIM, DIM, Q1D, Q1D, NE);
   const auto X = Reshape(x_.Read(), D1D, D1D, DIM, NE);
   const auto dP = Reshape(p_.Read(), DIM, DIM, DIM, DIM, Q1D, Q1D, NE);
   auto Y = Reshape(y_.ReadWrite(), D1D, D1D, DIM, NE);

   MFEM_FORALL_2D(e, NE, Q1D, Q1D, NBZ,
   {
      const int D1D = T_D1D ? T_D1D : d1d;
      const int Q1D = T_Q1D ? T_Q1D : q1d;
      constexpr int NBZ = T_NBZ ? T_NBZ : 1;
      constexpr int MQ1 = T_Q1D ? T_Q1D : T_MAX;
      constexpr int MD1 = T_D1D ? T_D1D : T_MAX;

      MFEM_SHARED double BG[2][MQ1*MD1];
      MFEM_SHARED double XY[2][NBZ][MD1*MD1];
      MFEM_SHARED double DQ[4][NBZ][MD1*MQ1];
      MFEM_SHARED double QQ[4][NBZ][MQ1*MQ1];

      kernels::LoadX<MD1,NBZ>(e,D1D,X,XY);
      kernels::LoadBG<MD1,MQ1>(D1D,Q1D,b,g,BG);

      kernels::GradX<MD1,MQ1,NBZ>(D1D,Q1D,BG,XY,DQ);
      kernels::GradY<MD1,MQ1,NBZ>(D1D,Q1D,BG,DQ,QQ);

      MFEM_FOREACH_THREAD(qy,y,Q1D)
      {
         MFEM_FOREACH_THREAD(qx,x,Q1D)
         {
            double A[4], B[4], C[4];
            const double *Jtr = &J(0,0,qx,qy,e);

            // Jrt = Jtr^{-1}
            double Jrt[4];
            kernels::CalcInverse<2>(Jtr, Jrt);

            double hX[4];
            kernels::PullGradXY<MQ1,NBZ>(qx,qy,QQ,hX);

            // A = X^T . Jrt
            kernels::Mult(2,2,2, hX, Jrt, A);

            // B = A : dP
            for (int r = 0; r < DIM; r++)
            {
               for (int c = 0; c < DIM; c++)
               {
                  B[r+2*c] = 0.0;
                  for (int i = 0; i < DIM; i++)
                  {
                     for (int j = 0; j < DIM; j++)
                     {
                        B[r+2*c] += dP(i,j,r,c,qx,qy,e) * A[i+2*j];
                     }
                  }
               }
            }

            // C = Jrt . B
            kernels::MultABt(2,2,2, Jrt, B, C);
            kernels::PushGradXY<MQ1,NBZ>(qx,qy,C,QQ);
         }
      }
      MFEM_SYNC_THREAD;
      kernels::LoadBGt<MD1,MQ1>(D1D,Q1D,b,g,BG);
      kernels::GradYt<MD1,MQ1,NBZ>(D1D,Q1D,BG,QQ,DQ);
      kernels::GradXt<MD1,MQ1,NBZ>(D1D,Q1D,BG,DQ,Y,e);
   });
}

void TMOP_Integrator::AddMultGradPA_2D(const Vector &X, const Vector &R,
                                       Vector &C) const
{
   const int N = PA.ne;
   const int D1D = PA.maps->ndof;
   const int Q1D = PA.maps->nqpt;
   const int id = (D1D << 4 ) | Q1D;
   const DenseTensor &J = PA.Jtr;
   const Array<double> &B = PA.maps->B;
   const Array<double> &G = PA.maps->G;
   const Vector &A = PA.A;

   if (!PA.setup)
   {
      PA.setup = true;
      AssembleGradPA_2D(X);
   }

   switch (id)
   {
      case 0x21: return AddMultGradPA_Kernel_2D<2,1,1>(N,B,G,J,A,R,C);
      case 0x22: return AddMultGradPA_Kernel_2D<2,2,1>(N,B,G,J,A,R,C);
      case 0x23: return AddMultGradPA_Kernel_2D<2,3,1>(N,B,G,J,A,R,C);
      case 0x24: return AddMultGradPA_Kernel_2D<2,4,1>(N,B,G,J,A,R,C);
      case 0x25: return AddMultGradPA_Kernel_2D<2,5,1>(N,B,G,J,A,R,C);
      case 0x26: return AddMultGradPA_Kernel_2D<2,6,1>(N,B,G,J,A,R,C);

      case 0x31: return AddMultGradPA_Kernel_2D<3,1,1>(N,B,G,J,A,R,C);
      case 0x32: return AddMultGradPA_Kernel_2D<3,2,1>(N,B,G,J,A,R,C);
      case 0x33: return AddMultGradPA_Kernel_2D<3,3,1>(N,B,G,J,A,R,C);
      case 0x34: return AddMultGradPA_Kernel_2D<3,4,1>(N,B,G,J,A,R,C);
      case 0x35: return AddMultGradPA_Kernel_2D<3,5,1>(N,B,G,J,A,R,C);
      case 0x36: return AddMultGradPA_Kernel_2D<3,6,1>(N,B,G,J,A,R,C);

      case 0x41: return AddMultGradPA_Kernel_2D<4,1,1>(N,B,G,J,A,R,C);
      case 0x42: return AddMultGradPA_Kernel_2D<4,2,1>(N,B,G,J,A,R,C);
      case 0x43: return AddMultGradPA_Kernel_2D<4,3,1>(N,B,G,J,A,R,C);
      case 0x44: return AddMultGradPA_Kernel_2D<4,4,1>(N,B,G,J,A,R,C);
      case 0x45: return AddMultGradPA_Kernel_2D<4,5,1>(N,B,G,J,A,R,C);
      case 0x46: return AddMultGradPA_Kernel_2D<4,6,1>(N,B,G,J,A,R,C);

      case 0x51: return AddMultGradPA_Kernel_2D<5,1,1>(N,B,G,J,A,R,C);
      case 0x52: return AddMultGradPA_Kernel_2D<5,2,1>(N,B,G,J,A,R,C);
      case 0x53: return AddMultGradPA_Kernel_2D<5,3,1>(N,B,G,J,A,R,C);
      case 0x54: return AddMultGradPA_Kernel_2D<5,4,1>(N,B,G,J,A,R,C);
      case 0x55: return AddMultGradPA_Kernel_2D<5,5,1>(N,B,G,J,A,R,C);
      case 0x56: return AddMultGradPA_Kernel_2D<5,6,1>(N,B,G,J,A,R,C);

      default:
      {
         constexpr int T_MAX = 8;
         MFEM_VERIFY(D1D <= MAX_D1D && Q1D <= MAX_Q1D, "Max size error!");
         return AddMultGradPA_Kernel_2D<0,0,0,T_MAX>(N,B,G,J,A,R,C,D1D,Q1D);
      }
   }
}

} // namespace mfem