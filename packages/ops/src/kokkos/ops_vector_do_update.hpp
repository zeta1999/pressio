/*
//@HEADER
// ************************************************************************
//
// ops_vector_do_update.hpp
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

#ifndef OPS_KOKKOS_OPS_VECTOR_DO_UPDATE_HPP_
#define OPS_KOKKOS_OPS_VECTOR_DO_UPDATE_HPP_

#include <KokkosBlas1_axpby.hpp>

namespace pressio{ namespace ops{

//----------------------------------------------------------------------
// computing:  V = a * V + b * V1
//----------------------------------------------------------------------
template<typename T, typename scalar_t>
::pressio::mpl::enable_if_t<
  ::pressio::containers::predicates::is_vector_wrapper_kokkos<T>::value and
  !::pressio::containers::predicates::is_diag_expression<T>::value
  >
do_update(T & v, const scalar_t & a,
	  const T & v1, const scalar_t & b)
{
  // v = a*v + b * v1
  KokkosBlas::axpby(b, *v1.data(), a, *v.data());
}

template<typename T, typename scalar_t>
::pressio::mpl::enable_if_t<
  ::pressio::containers::predicates::is_vector_wrapper_kokkos<T>::value and
  !::pressio::containers::predicates::is_diag_expression<T>::value
  >
do_update(T & v, const T & v1, const scalar_t & b)
{
  // v = b*v1
  constexpr auto zero = ::pressio::utils::constants<scalar_t>::zero();
  KokkosBlas::axpby(b, *v1.data(), zero, *v.data());
}


// specialize for when we have vector-like diag expressions,
// because for these we currrently do not support the data() method
// because we need to figure out how to create a native diagonal
// expression in kokkos. we only specialize this case and not the ones below
// because currently the diagonal view is only used with such an op.
template<typename T1, typename T2, typename scalar_t>
::pressio::mpl::enable_if_t<
  ::pressio::containers::predicates::is_vector_wrapper_kokkos<T1>::value and
  ::pressio::containers::predicates::is_vector_wrapper_kokkos<T2>::value and
  ::pressio::containers::predicates::is_diag_expression<T1>::value and
  ::pressio::containers::predicates::is_diag_expression<T2>::value
  >
do_update(T1 & v, const scalar_t & a,
	  const T2 & v1, const scalar_t & b)
{
  static_assert
    (std::is_same<
     typename ::pressio::containers::details::traits<T1>::memory_space,
     typename ::pressio::containers::details::traits<T2>::memory_space>::value,
     "Kokkos diag Expressions do not have matching memory space");

  // static_assert
  //   (Kokkos::SpaceAccessibility<
  //    Kokkos::HostSpace,
  //    typename ::pressio::containers::details::traits<T1>::memory_space>::accessible,
  //    "Kokkos diag Expressions must be host accessible.");

  auto kvLhs = *v.getUnderlyingObject().data();
  auto kvRhs = *v1.getUnderlyingObject().data();
  Kokkos::parallel_for
    (v.extent(0),
     KOKKOS_LAMBDA (const int& i) {
      kvLhs(i,i) = a*kvLhs(i,i) + b*kvRhs(i,i);
    });
}



//----------------------------------------------------------------------
//  overloads for computing this: V = a * V + b * V1 + c * V2
//----------------------------------------------------------------------
template<typename T, typename scalar_t>
::pressio::mpl::enable_if_t<
  ::pressio::containers::predicates::is_vector_wrapper_kokkos<T>::value
  >
do_update(T & v, const scalar_t &a,
	  const T & v1, const scalar_t &b,
	  const T & v2, const scalar_t &c)
{
  using view_t = typename ::pressio::containers::details::traits<T>::wrapped_t;
  using fnctr_t = ::pressio::ops::impl::DoUpdateTwoTermsFunctor<view_t, scalar_t>;
  fnctr_t F(*v.data(), *v1.data(), *v2.data(), a, b, c);
  Kokkos::parallel_for(v.extent(0), F);
}

template<typename T, typename scalar_t>
::pressio::mpl::enable_if_t<
  ::pressio::containers::predicates::is_vector_wrapper_kokkos<T>::value
  >
do_update(T & v,
	  const T & v1, const scalar_t &b,
	  const T & v2, const scalar_t &c)
{
  using view_t = typename ::pressio::containers::details::traits<T>::wrapped_t;
  using fnctr_t = ::pressio::ops::impl::DoUpdateTwoTermsFunctor<view_t, scalar_t>;
  fnctr_t F(*v.data(), *v1.data(), *v2.data(), b, c);
  Kokkos::parallel_for(v.extent(0), F);
}


//----------------------------------------------------------------------
//  overloads for computing:
//	V = a * V + b * V1 + c * V2 + d * V3
//----------------------------------------------------------------------
template<typename T, typename scalar_t>
::pressio::mpl::enable_if_t<
  ::pressio::containers::predicates::is_vector_wrapper_kokkos<T>::value
  >
do_update(T & v, const scalar_t &a,
	  const T & v1, const scalar_t &b,
	  const T & v2, const scalar_t &c,
	  const T & v3, const scalar_t &d)
{
  using view_t = typename ::pressio::containers::details::traits<T>::wrapped_t;
  using fnctr_t = ::pressio::ops::impl::DoUpdateThreeTermsFunctor<view_t, scalar_t>;
  fnctr_t F(*v.data(), *v1.data(), *v2.data(), *v3.data(), a, b, c, d);
  Kokkos::parallel_for(v.extent(0), F);
}

template<typename T, typename scalar_t>
::pressio::mpl::enable_if_t<
  ::pressio::containers::predicates::is_vector_wrapper_kokkos<T>::value
  >
do_update(T & v,
	  const T & v1, const scalar_t &b,
	  const T & v2, const scalar_t &c,
	  const T & v3, const scalar_t &d)
{
  using view_t = typename ::pressio::containers::details::traits<T>::wrapped_t;
  using fnctr_t = ::pressio::ops::impl::DoUpdateThreeTermsFunctor<view_t, scalar_t>;
  fnctr_t F(*v.data(), *v1.data(), *v2.data(), *v3.data(), b, c, d);
  Kokkos::parallel_for(v.extent(0), F);
}


//----------------------------------------------------------------------
//  overloads for computing:
//	V = a * V + b * V1 + c * V2 + d * V3 + e * V4
//----------------------------------------------------------------------
template<typename T, typename scalar_t>
::pressio::mpl::enable_if_t<
  ::pressio::containers::predicates::is_vector_wrapper_kokkos<T>::value
  >
do_update(T & v, const scalar_t &a,
	  const T & v1, const scalar_t &b,
	  const T & v2, const scalar_t &c,
	  const T & v3, const scalar_t &d,
	  const T & v4, const scalar_t &e)
{
  using view_t = typename ::pressio::containers::details::traits<T>::wrapped_t;
  using fnctr_t = ::pressio::ops::impl::DoUpdateFourTermsFunctor<view_t, scalar_t>;
  fnctr_t F(*v.data(), *v1.data(), *v2.data(), *v3.data(), *v4.data(), a, b, c, d, e);
  Kokkos::parallel_for(v.extent(0), F);
}

template<typename T, typename scalar_t>
::pressio::mpl::enable_if_t<
  ::pressio::containers::predicates::is_vector_wrapper_kokkos<T>::value
  >
do_update(T & v,
	  const T & v1, const scalar_t &b,
	  const T & v2, const scalar_t &c,
	  const T & v3, const scalar_t &d,
	  const T & v4, const scalar_t &e)
{
  using view_t = typename ::pressio::containers::details::traits<T>::wrapped_t;
  using fnctr_t = ::pressio::ops::impl::DoUpdateFourTermsFunctor<view_t, scalar_t>;
  fnctr_t F(*v.data(), *v1.data(), *v2.data(), *v3.data(), *v4.data(), b, c, d, e);
  Kokkos::parallel_for(v.extent(0), F);
}


}}//end namespace pressio::ops
#endif  // OPS_KOKKOS_OPS_VECTOR_DO_UPDATE_HPP_
