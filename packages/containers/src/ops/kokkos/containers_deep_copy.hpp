/*
//@HEADER
// ************************************************************************
//
// containers_deep_copy.hpp
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

#ifdef PRESSIO_ENABLE_TPL_KOKKOS
#ifndef CONTAINERS_CONTAINER_OPS_KOKKOS_DEEP_COPY_HPP_
#define CONTAINERS_CONTAINER_OPS_KOKKOS_DEEP_COPY_HPP_

namespace pressio{ namespace containers{ namespace ops{

//--------------------------------------------------------------------------
// for wrappers, we overload the = operator
//--------------------------------------------------------------------------
template<
  typename T,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_vector_wrapper_kokkos<T>::value
    > * = nullptr
  >
void deep_copy(const T & src, T & dest){
  dest = src;
}

template<
  typename T1, typename T2,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_vector_wrapper_kokkos<T1>::value and
    ::pressio::containers::meta::is_expression<T2>::value
    > * = nullptr
  >
void deep_copy(const T1 & src, T2 & dest){
  Kokkos::deep_copy(dest(), *src.data());
}

template<
  typename T1, typename T2,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_expression<T1>::value and
    ::pressio::containers::meta::is_vector_wrapper_kokkos<T2>::value
    > * = nullptr
  >
void deep_copy(const T1 & src, T2 & dest){
  Kokkos::deep_copy( *dest.data(), src() );
}

template<
  typename T1,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_expression<T1>::value and
    (::pressio::containers::meta::is_vector_wrapper_kokkos<
     typename ::pressio::containers::details::traits<T1>::data_t
      >::value or
     ::pressio::containers::meta::is_matrix_wrapper_kokkos<
     typename ::pressio::containers::details::traits<T1>::data_t
     >::value)
    > * = nullptr
  >
void deep_copy(const T1 & src, T1 & dest){
  Kokkos::deep_copy(dest(), src());
}

}}}//end namespace pressio::containers::ops
#endif
#endif
