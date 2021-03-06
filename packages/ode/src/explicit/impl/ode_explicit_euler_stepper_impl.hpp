/*
//@HEADER
// ************************************************************************
//
// ode_explicit_euler_stepper_impl.hpp
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

#ifndef ODE_EXPLICIT_IMPL_ODE_EXPLICIT_EULER_STEPPER_IMPL_HPP_
#define ODE_EXPLICIT_IMPL_ODE_EXPLICIT_EULER_STEPPER_IMPL_HPP_

namespace pressio{ namespace ode{ namespace explicitmethods{ namespace impl{

template<
  typename scalar_type,
  typename state_type,
  typename system_type,
  typename velocity_type,
  typename velocity_policy_type,
  typename standard_velocity_policy_type,
  typename ops_t
  >
class ExplicitEulerStepperImpl
{

public:
  static constexpr bool is_implicit = false;
  static constexpr bool is_explicit = true;

  static constexpr types::stepper_order_t order_value = 1;

private:
  using velocity_storage_t = VelocitiesContainer<velocity_type, 1>;
  using system_wrapper_t   = ::pressio::ode::impl::OdeSystemWrapper<system_type>;

  system_wrapper_t systemObj_;

  typename std::conditional<
    std::is_same<velocity_policy_type, standard_velocity_policy_type>::value,
    const standard_velocity_policy_type,
    const velocity_policy_type &
  >::type policy_;

  velocity_storage_t veloAuxStorage_;
  const ops_t * udOps_ = nullptr;

public:
  ExplicitEulerStepperImpl() = delete;
  ~ExplicitEulerStepperImpl() = default;
  ExplicitEulerStepperImpl(const ExplicitEulerStepperImpl &) = delete;
  ExplicitEulerStepperImpl & operator=(const ExplicitEulerStepperImpl &) = delete;
  ExplicitEulerStepperImpl(ExplicitEulerStepperImpl &&) = delete;
  ExplicitEulerStepperImpl & operator=(ExplicitEulerStepperImpl &&) = delete;

  // the following cnstr is enabled if we are using pressio ops
  template <
    typename _ops_t = ops_t,
    mpl::enable_if_t< std::is_void<_ops_t>::value, int > = 0
  >
  ExplicitEulerStepperImpl(const state_type & state,
                           const system_type & systemObj,
                           const velocity_policy_type & policy)
    : systemObj_(systemObj),
      policy_(policy),
      veloAuxStorage_(policy_.create(systemObj))
  {}

  // the following cnstr is enabled if we are using user-defined ops
  template <
    typename _ops_t = ops_t,
    mpl::enable_if_t< !std::is_void<_ops_t>::value, int > = 0
  >
  ExplicitEulerStepperImpl(const state_type & state,
                           const system_type & systemObj,
                           const velocity_policy_type & policy,
                           const _ops_t & udOps)
    : systemObj_(systemObj),
      policy_(policy),
      veloAuxStorage_(policy_.create(systemObj)),
      udOps_(&udOps)
  {}

  // the following cnstr is only enabled if
  // policy is default constructible and we are using pressio ops
  template <
    typename _vel_pol_t = velocity_policy_type,
    typename _ops_t = ops_t,
    mpl::enable_if_t<
      std::is_void<_ops_t>::value and
      std::is_default_constructible<_vel_pol_t>::value,
      int > = 0
  >
  ExplicitEulerStepperImpl(const state_type & state,
                           const system_type & systemObj)
    : systemObj_(systemObj),
      policy_(), // default construct policy
      veloAuxStorage_(policy_.create(systemObj))
  {}

  // the following cnstr is only enabled if
  // policy is default constructible and we are using user-defined ops
  template <
    typename _vel_pol_t = velocity_policy_type,
    typename _ops_t = ops_t,
    mpl::enable_if_t<
      !std::is_void<_ops_t>::value and
      std::is_default_constructible<_vel_pol_t>::value,
      int > = 0
  >
  ExplicitEulerStepperImpl(const state_type & state,
                           const system_type & systemObj,
                           const _ops_t & udOps)
    : systemObj_(systemObj),
      policy_(), // default construct policy
      veloAuxStorage_(policy_.create(systemObj)),
      udOps_(&udOps)
  {}

public:
  types::stepper_order_t order() const
  {
    return order_value;
  }

  /* user does NOT provide custom ops, so we use ops */
  template<typename _ops_t = ops_t>
  mpl::enable_if_t< std::is_void<_ops_t>::value >
  doStep(state_type & stateInOut,
	 const scalar_type & time,
	 const scalar_type & dt,
	 const types::step_t & step)
  {
    auto & auxRhs0 = veloAuxStorage_(0);
    //eval RHS
    policy_.compute(stateInOut, auxRhs0, systemObj_.get(), time);
    // y = y + dt * rhs
    constexpr auto one  = ::pressio::utils::constants<scalar_type>::one();
    ::pressio::ops::do_update(stateInOut, one, auxRhs0, dt);
  }

  /* user provides custom ops */
  template<typename _ops_t = ops_t>
  mpl::enable_if_t<!std::is_void<_ops_t>::value>
  doStep(state_type & stateInOut,
	 const scalar_type & time,
	 const scalar_type & dt,
	 const types::step_t & step)
  {
    auto & auxRhs0 = veloAuxStorage_(0);
    policy_.compute(stateInOut, auxRhs0, systemObj_.get(), time);
    // y = y + dt * rhs
    constexpr auto one  = ::pressio::utils::constants<scalar_type>::one();
    udOps_->do_update(*stateInOut.data(), one, *auxRhs0.data(), dt);
  }
};

}}}}//end namespace pressio::ode::explicitmethods::impl
#endif  // ODE_EXPLICIT_IMPL_ODE_EXPLICIT_EULER_STEPPER_IMPL_HPP_
