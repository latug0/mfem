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

#ifndef MFEM_BACKENDS_OCCA_ARRAY_HPP
#define MFEM_BACKENDS_OCCA_ARRAY_HPP

#include "../../config/config.hpp"
#if defined(MFEM_USE_BACKENDS) && defined(MFEM_USE_OCCA)

#include <occa.hpp>
#include "layout.hpp"
#include "../base/array.hpp"

namespace mfem
{

namespace occa
{

class Array : public virtual PArray
{
protected:
   //
   // Inherited fields
   //
   // DLayout layout;

   // Always true: Size()*item_size == slice.size() <= data.size()
   mutable ::occa::memory data, slice;

   //
   // Virtual interface
   //

   virtual void *DoGetData() const { return GetBuffer(); }

   virtual PArray *DoClone(bool copy_data, void **buffer,
                           std::size_t item_size) const;

   virtual int DoResize(PLayout &new_layout, void **buffer,
                        std::size_t item_size);

   virtual void *DoPullData(void *buffer, std::size_t item_size);

   virtual void DoFill(const void *value_ptr, std::size_t item_size);

   virtual void DoPushData(const void *src_buffer, std::size_t item_size);

   virtual void DoAssign(const PArray &src, std::size_t item_size);

   //
   // Auxiliary methods
   //

   inline void *GetBuffer() const;

   inline int ResizeData(const Layout *lt, std::size_t item_size);

   template <typename T>
   inline void OccaFill(const T *val_ptr)
   { ::occa::linalg::operator_eq<T>(slice, *val_ptr); }

public:
   Array(Layout &lt, std::size_t item_size)
      : PArray(lt),
        data(lt.Alloc(lt.Size()*item_size)),
        slice(data)
   { }

   virtual ~Array() { }

   inline void MakeRef(Array &master);

   Layout &OccaLayout() const
   { return *static_cast<Layout *>(layout.Get()); }

   ::occa::memory &OccaMem() { return slice; }
   const ::occa::memory &OccaMem() const { return slice; }
};


//
// Inline methods
//

inline void *Array::GetBuffer() const
{
   if (!slice.getDevice().hasSeparateMemorySpace())
   {
      return slice.ptr();
   }
   return NULL;
}

inline int Array::ResizeData(const Layout *lt, std::size_t item_size)
{
   const std::size_t new_bytes = lt->Size()*item_size;
   if (data.size() < new_bytes ||
       data.getDHandle() != lt->OccaEngine().GetDevice().getDHandle())
   {
      data = lt->Alloc(new_bytes);
      slice = data;
      // If memory allocation fails - an exception is thrown.
   }
   else if (slice.size() != new_bytes)
   {
      slice = data.slice(0, new_bytes);
   }
   return 0;
}

inline void Array::MakeRef(Array &master)
{
   layout = master.layout;
   data = master.data;
   slice = master.slice;
}

} // namespace mfem::occa

} // namespace mfem

#endif // defined(MFEM_USE_BACKENDS) && defined(MFEM_USE_OCCA)

#endif // MFEM_BACKENDS_OCCA_ARRAY_HPP
