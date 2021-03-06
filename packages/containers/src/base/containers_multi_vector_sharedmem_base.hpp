/*
//@HEADER
// ************************************************************************
//
// containers_multi_vector_sharedmem_base.hpp
//                     		  Pressio
//                             Copyright 2019
//    National Technology & Engineering Solutions of Sandia, LLC (NTESS)
//
// Under the terms of Contract DE-NA0003525 with NTESS, the
// U.S. Government retains certain rights in this software.
//
// Pressio is licensed under BSD-3-Clause terms of use:
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its
// contributors may be used to endorse or promote products derived
// from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Francesco Rizzi (fnrizzi@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef CONTAINERS_BASE_CONTAINERS_MULTI_VECTOR_SHAREDMEM_BASE_HPP_
#define CONTAINERS_BASE_CONTAINERS_MULTI_VECTOR_SHAREDMEM_BASE_HPP_

namespace pressio{ namespace containers{

template<typename derived_type>
class MultiVectorSharedMemBase : public ContainerSharedMemBase<derived_type>
{
  using traits = ::pressio::containers::details::traits<derived_type>;

  static_assert(traits::is_shared_mem==1,
		"OOPS: the concrete vector type inheriting from sharedMem base is not shared mem!");

  using this_t = MultiVectorSharedMemBase<derived_type>;
  using sc_t  = typename traits::scalar_t;
  using ord_t = typename traits::ordinal_t;
  friend derived_type;

public:
  ord_t numVectors() const{
    return static_cast<const derived_type &>(*this).numVectors();
  }
};//end class

}}//end namespace pressio::containers
#endif  // CONTAINERS_BASE_CONTAINERS_MULTI_VECTOR_SHAREDMEM_BASE_HPP_
