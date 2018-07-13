// Copyright (c) 2010, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-443211. All Rights
// reserved. See file COPYRIGHT for details.
//
// This file is part of the MFEM library. For more information and source code
// availability see http://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free
// Software Foundation) version 2.1 dated February 1999.

#include "../../config/config.hpp"
#if defined(MFEM_USE_BACKENDS) && defined(MFEM_USE_OMP)

#include "backend.hpp"
#include "bilinearform.hpp"
#include "adiffusioninteg.hpp"

namespace mfem
{

namespace omp
{

BilinearForm::~BilinearForm()
{
   // Make sure all integrators free their data
   for (int i = 0; i < tbfi.Size(); i++) delete tbfi[i];

   delete element_matrices;
}

void BilinearForm::TransferIntegrators()
{
   mfem::Array<mfem::BilinearFormIntegrator*> &dbfi = *bform->GetDBFI();
   for (int i = 0; i < dbfi.Size(); i++)
   {
      std::string integ_name(dbfi[i]->Name());
      Coefficient *scal_coeff = dbfi[i]->GetScalarCoefficient();
      // ConstantCoefficient *const_coeff =
      //    dynamic_cast<ConstantCoefficient*>(scal_coeff);
      // // TODO: other types of coefficients ...
      // double val = const_coeff ? const_coeff->constant : 1.0;

      if (integ_name == "(undefined)")
      {
         MFEM_ABORT("BilinearFormIntegrator does not define Name()");
      }
      else if (integ_name == "diffusion")
      {
         switch (OmpEngine().IntegType())
         {
         case Acrotensor:
            tbfi.Append(new AcroDiffusionIntegrator(*scal_coeff, bform->FESpace()->Get_PFESpace()->As<FiniteElementSpace>()));
            break;
         default:
            mfem_error("integrator is not supported for any MultType");
            break;
         }
      }
      else
      {
         MFEM_ABORT("BilinearFormIntegrator [Name() = " << integ_name
                    << "] is not supported");
      }
   }
}

void BilinearForm::InitRHS(const mfem::Array<int> &ess_tdof_list,
                           mfem::Vector &mfem_x, mfem::Vector &mfem_b,
                           mfem::OperatorHandle &A,
                           mfem::Vector &mfem_X, mfem::Vector &mfem_B,
                           int copy_interior) const
{
   const mfem::Operator *P = GetProlongation();
   const mfem::Operator *R = GetRestriction();

   if (P)
   {
      // Variational restriction with P
      mfem_B.Resize(P->InLayout());
      P->MultTranspose(mfem_b, mfem_B);
      mfem_X.Resize(R->OutLayout());
      R->Mult(mfem_x, mfem_X);
   }
   else
   {
      // rap, X and B point to the same data as this, x and b
      mfem_X.MakeRef(mfem_x);
      mfem_B.MakeRef(mfem_b);
   }

   if (A.Type() != mfem::Operator::ANY_TYPE)
   {
      A.EliminateBC(mat_e, ess_tdof_list, mfem_X, mfem_B);
   }

   if (!copy_interior && ess_tdof_list.Size() > 0)
   {
      Vector &X = mfem_X.Get_PVector()->As<Vector>();
      const Array &constraint_list = ess_tdof_list.Get_PArray()->As<Array>();

      double *X_data = X.GetData<double>();
      const int* constraint_data = constraint_list.GetData<int>();

      Vector subvec(constraint_list.OmpLayout());
      double *subvec_data = subvec.GetData<double>();

      const std::size_t num_constraint = constraint_list.Size();
      const bool use_target = constraint_list.ComputeOnDevice();
      const bool use_parallel = (use_target || num_constraint > 1000);

      // This operation is a general version of mfem::Vector::SetSubVectorComplement()
      // {
#pragma omp target teams distribute parallel for        \
   map(to: subvec_data, constraint_data, X_data)        \
   if (target: use_target)                              \
   if (parallel: use_parallel)
      for (std::size_t i = 0; i < num_constraint; i++) subvec_data[i] = X_data[constraint_data[i]];

      X.Fill(0.0);

#pragma omp target teams distribute parallel for        \
   map(to: X_data, constraint_data, subvec_data)        \
   if (target: use_target)                              \
   if (parallel: use_parallel)
      for (std::size_t i = 0; i < num_constraint; i++) X_data[constraint_data[i]] = subvec_data[i];
      // }
   }

   if (A.Type() == mfem::Operator::ANY_TYPE)
   {
      ConstrainedOperator *A_constrained = static_cast<ConstrainedOperator*>(A.Ptr());
      A_constrained->EliminateRHS(mfem_X, mfem_B);
   }
}


bool BilinearForm::Assemble()
{
   if (!has_assembled)
   {
      TransferIntegrators();
      has_assembled = true;
   }

   return true;
}

void BilinearForm::ComputeElementMatrices()
{
   // Only called if performing full assembly
   const int nelements = trial_fes->GetFESpace()->GetNE();
   const int trial_ndofs = trial_fes->GetFESpace()->GetFE(0)->GetDof() * trial_fes->GetFESpace()->GetVDim();
   const int test_ndofs = test_fes->GetFESpace()->GetFE(0)->GetDof() * test_fes->GetFESpace()->GetVDim();
   const std::size_t length = nelements * trial_ndofs * test_ndofs;

   if (!element_matrices) element_matrices = new mfem::Vector(*(new Layout(OmpEngine(), length)));
   else element_matrices->Push();

   element_matrices->Fill(0.0);
   Vector &elmats = element_matrices->Get_PVector()->As<Vector>();

   tbfi[0]->ComputeElementMatrices(elmats);

   if (tbfi.Size() > 1)
   {
      for (int k = 1; k < tbfi.Size(); k++)
      {
         tbfi[k]->ComputeElementMatrices(elmats);
      }
   }
}

void BilinearForm::FormSystemMatrix(const mfem::Array<int> &ess_tdof_list,
                                    mfem::OperatorHandle &A)
{
   if (A.Type() == mfem::Operator::ANY_TYPE)
   {
      // FIXME: Support different test and trial spaces (MixedBilinearForm)
      const mfem::Operator *P = GetProlongation();

      mfem::Operator *rap = this;
      if (P != NULL) rap = new mfem::RAPOperator(*P, *this, *P);

      A.Reset(new ConstrainedOperator(rap, ess_tdof_list, (rap != this)));

      return;
   }
   else
   {
      // ASSUMPTION: some sort of sparse matrix
      // Compute the local matrices (stored in bform->element_matrices
      ComputeElementMatrices();
      bform->AllocateMatrix();
      mfem::SparseMatrix &mat = bform->SpMat();

      element_matrices->Pull();
      double *data = element_matrices->GetData();

      const bool skip_zeros = true;
      mfem::Array<int> tr_vdofs, te_vdofs;
      for (int i = 0; i < trial_fes->GetFESpace()->GetNE(); i++)
      {
         trial_fes->GetFESpace()->GetElementVDofs(i, tr_vdofs);
         test_fes->GetFESpace()->GetElementVDofs(i, te_vdofs);
         const mfem::DenseMatrix elmat(data, te_vdofs.Size(), tr_vdofs.Size());
         mat.AddSubMatrix(te_vdofs, tr_vdofs, elmat, skip_zeros);
         data += tr_vdofs.Size() * te_vdofs.Size();
      }
   }

   if (A.Type() == mfem::Operator::MFEM_SPARSEMAT)
   {
      // This works because the FormSystemMatrix call with an explicit
      // SparseMatrix doesnt call the backend version... This might
      // change in the future.
      bform->FormSystemMatrix(ess_tdof_list, static_cast<mfem::SparseMatrix&>(*A.Ptr()));
   }
#ifdef MFEM_USE_MPI
   else if (A.Type() == mfem::Operator::Hypre_ParCSR)
   {
      mfem::SparseMatrix &mat = bform->SpMat();
      mfem::ParBilinearForm *pbform = dynamic_cast<mfem::ParBilinearForm*>(bform);

      const bool skip_zeros = false;
      mat.Finalize(skip_zeros);

      // -------- FOR SOME VERY AGGREVATING REASON THIS DOESN'T WORK ---------
      // mfem::ParFiniteElementSpace *pfes = pbform->ParFESpace();
      // OperatorHandle dA(Operator::Hypre_ParCSR);
      // // construct a parallel block-diagonal matrix 'A' based on 'a'
      // dA.MakeSquareBlockDiag(pfes->GetComm(), *engine->MakeLayout(pfes->GlobalTrueVSize()),
      //                        pfes->GetDofOffsets(), &mat);
      // OperatorHandle Ph(pfes->Dof_TrueDof_Matrix());
      // A.MakePtAP(dA, Ph);
      // A.SetOperatorOwner(false);
      // -------- BUT THIS DOES ---------
      pbform->ParallelAssemble(A, &mat);
      A.SetOperatorOwner(false);
      // ---------------------
      mat.Clear();
      mat_e.Clear();
      std::cout << "operator size (FormSystemMatrix): " << A.Ptr()->InLayout()->Size() << " " << A.Ptr()->OutLayout()->Size() << std::endl;

      mat_e.EliminateRowsCols(A, ess_tdof_list);
   }
#endif
   else
   {
      MFEM_ABORT("Operator::Type is not supported, type = " << A.Type());
   }
}

void BilinearForm::FormLinearSystem(const mfem::Array<int> &ess_tdof_list,
                                    mfem::Vector &x, mfem::Vector &b,
                                    mfem::OperatorHandle &A, mfem::Vector &X, mfem::Vector &B,
                                    int copy_interior)
{
   FormSystemMatrix(ess_tdof_list, A);
   std::cout << "operator size (FormLinearSystem 1): " << A.Ptr()->InLayout()->Size() << " " << A.Ptr()->OutLayout()->Size() << std::endl;
   InitRHS(ess_tdof_list, x, b, A, X, B, copy_interior);
}

void BilinearForm::RecoverFEMSolution(const mfem::Vector &X, const mfem::Vector &b,
                                      mfem::Vector &x)
{
   const mfem::Operator *P = GetProlongation();
   if (P)
   {
      // Apply conforming prolongation
      x.Resize(P->OutLayout());
      P->Mult(X, x);
   }
   // Otherwise X and x point to the same data
}

void BilinearForm::Mult(const mfem::Vector &x, mfem::Vector &y) const
{
   trial_fes->ToEVector(x.Get_PVector()->As<Vector>(), x_local);

   y_local.Fill<double>(0.0);
   for (int i = 0; i < tbfi.Size(); i++) tbfi[i]->MultAdd(x_local, y_local);

   test_fes->ToLVector(y_local, y.Get_PVector()->As<Vector>());
}

void BilinearForm::MultTranspose(const mfem::Vector &x, mfem::Vector &y) const
{ mfem_error("mfem::omp::BilinearForm::MultTranspose() is not supported!"); }


ConstrainedOperator::ConstrainedOperator(mfem::Operator *A_,
                                         const mfem::Array<int> &constraint_list_,
                                         bool own_A_)
   : Operator(A_->InLayout()->As<Layout>()),
     A(A_),
     own_A(own_A_),
     // FIXME: @dudouit1 has a general fix for this
     constraint_list(constraint_list_.Get_PArray()->As<Array>()),
     z(OutLayout()->As<Layout>()),
     w(OutLayout()->As<Layout>()),
     mfem_z((z.DontDelete(), z)),
     mfem_w((w.DontDelete(), w)) { }

void ConstrainedOperator::EliminateRHS(const mfem::Vector &mfem_x, mfem::Vector &mfem_b) const
{
   w.Fill<double>(0.0);

   const Vector &x = mfem_x.Get_PVector()->As<Vector>();
   Vector &b = mfem_b.Get_PVector()->As<Vector>();

   const double *x_data = x.GetData<double>();
   double *b_data = b.GetData<double>();
   double *w_data = w.GetData<double>();
   const int* constraint_data = constraint_list.GetData<int>();

   const std::size_t num_constraint = constraint_list.Size();
   const bool use_target = constraint_list.ComputeOnDevice();
   const bool use_parallel = (use_target || num_constraint > 1000);

   if (num_constraint > 0)
   {
#pragma omp target teams distribute parallel for             \
   map(to: w_data, constraint_data, x_data)                  \
   if (target: use_target)                                   \
   if (parallel: use_parallel)
      for (std::size_t i = 0; i < num_constraint; i++)
         w_data[constraint_data[i]] = x_data[constraint_data[i]];
   }

   A->Mult(mfem_w, mfem_z);

   b.Axpby<double>(1.0, b, -1.0, z);

   if (num_constraint > 0)
   {
#pragma omp target teams distribute parallel for        \
   map(to: b_data, constraint_data, x_data)             \
   if (target: use_target)                              \
   if (parallel: use_parallel)
      for (std::size_t i = 0; i < num_constraint; i++)
         b_data[constraint_data[i]] = x_data[constraint_data[i]];
   }
}

void ConstrainedOperator::Mult(const mfem::Vector &mfem_x, mfem::Vector &mfem_y) const
{
   if (constraint_list.Size() == 0)
   {
      A->Mult(mfem_x, mfem_y);
      return;
   }

   const Vector &x = mfem_x.Get_PVector()->As<Vector>();
   Vector &y = mfem_y.Get_PVector()->As<Vector>();

   const double *x_data = x.GetData<double>();
   double *y_data = y.GetData<double>();
   double *z_data = z.GetData<double>();
   const int* constraint_data = constraint_list.GetData<int>();

   const std::size_t num_constraint = constraint_list.Size();
   const bool use_target = constraint_list.ComputeOnDevice();
   const bool use_parallel = (use_target || num_constraint > 1000);

   z.Assign<double>(x); // z = x

   // z[constraint_list] = 0.0
#pragma omp target teams distribute parallel for             \
   map(to: z_data, constraint_data)                          \
   if (target: use_target)                                   \
   if (parallel: use_parallel)
   for (std::size_t i = 0; i < num_constraint; i++)
      z_data[constraint_data[i]] = 0.0;

   // y = A * z
   A->Mult(mfem_z, mfem_y);

   // y[constraint_list] = x[constraint_list]
#pragma omp target teams distribute parallel for             \
   map(to: y_data, constraint_data, x_data)                  \
   if (target: use_target)                                   \
   if (parallel: use_parallel)
   for (std::size_t i = 0; i < num_constraint; i++)
      y_data[constraint_data[i]] = x_data[constraint_data[i]];
}

// Destructor: destroys the unconstrained Operator @a A if @a own_A is true.
ConstrainedOperator::~ConstrainedOperator()
{
   if (own_A) delete A;
}


} // namespace mfem::omp

} // namespace mfem

#endif // defined(MFEM_USE_BACKENDS) && defined(MFEM_USE_OMP)
