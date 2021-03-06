/*
//@HEADER
// ************************************************************************
//
// ode_implicit_residual_policy.hpp
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

#ifndef ODE_WILL_BE_CONCEPTS_POLICIES_ODE_IMPLICIT_RESIDUAL_POLICY_HPP_
#define ODE_WILL_BE_CONCEPTS_POLICIES_ODE_IMPLICIT_RESIDUAL_POLICY_HPP_

namespace pressio{ namespace ode{ namespace concepts {

template<
  typename T,
  typename tag,
  std::size_t numPrevStates,
  typename state_t,
  typename residual_t,
  typename system_t,
  typename scalar_t,
  typename enable = void
  >
struct implicit_residual_policy : std::false_type{};


template<
  typename T,
  typename tag,
  std::size_t numPrevStates,
  typename state_t,
  typename residual_t,
  typename system_t,
  typename scalar_t
  >
struct implicit_residual_policy<
  T, tag, numPrevStates, state_t, residual_t, system_t, scalar_t,
  ::pressio::mpl::enable_if_t<
    // is callable with two args
    std::is_same<
      residual_t,
      decltype
      (
       std::declval<T const>().create(std::declval<system_t const &>())
       )
      >::value
    and

    // is callable with six
    std::is_void<
      decltype
      (
       std::declval<T const>().template compute
       <tag>(
	     std::declval<state_t const &>(),
	     std::declval<::pressio::ode::AuxStatesManager<state_t, numPrevStates> const &>(),
	     std::declval<system_t const &>(),
	     std::declval<scalar_t const &>(),
	     std::declval<scalar_t const &>(),
	     std::declval<::pressio::ode::types::step_t>(),
	     std::declval<residual_t &>(),
	     ::pressio::Norm::Undefined,
	     std::declval<scalar_t &>()
	     )
       )
      >::value
    >
  > : std::true_type{};
//------------------------------------------------------------------


template<typename T, typename ... args>
using implicit_euler_residual_policy =
  implicit_residual_policy<
  T, ::pressio::ode::implicitmethods::Euler, 1, args...>;

template<typename T, typename ... args>
using implicit_bdf2_residual_policy =
  implicit_residual_policy<
  T, ::pressio::ode::implicitmethods::BDF2, 2, args...>;

}}} // namespace pressio::ode::concepts
#endif  // ODE_WILL_BE_CONCEPTS_POLICIES_ODE_IMPLICIT_RESIDUAL_POLICY_HPP_
