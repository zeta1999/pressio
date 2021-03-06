/*
//@HEADER
// ************************************************************************
//
// rom_lspg_unsteady_problem_discrete_time_api.hpp
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

#ifndef ROM_LSPG_IMPL_UNSTEADY_DISCRETE_TIME_API_ROM_LSPG_UNSTEADY_PROBLEM_DISCRETE_TIME_API_HPP_
#define ROM_LSPG_IMPL_UNSTEADY_DISCRETE_TIME_API_ROM_LSPG_UNSTEADY_PROBLEM_DISCRETE_TIME_API_HPP_

#include "./policies/rom_lspg_unsteady_residual_policy_discrete_time_api.hpp"
#include "./policies/rom_lspg_unsteady_jacobian_policy_discrete_time_api.hpp"

#include "./traits/rom_lspg_unsteady_common_traits_discrete_time_api.hpp"
#include "./traits/rom_lspg_unsteady_default_problem_traits_discrete_time_api.hpp"

namespace pressio{ namespace rom{ namespace lspg{ namespace impl{ namespace unsteady{

template <
  template <class, class, class, class ...> class lspg_type,
  typename stepper_tag,
  typename fom_system_type,
  typename lspg_state_type,
  typename ...Args
  >
class ProblemDiscreteTimeApi
{

public:
  // define the type holding types for the problem
  using lspg_problem_t = lspg_type<stepper_tag, fom_system_type, lspg_state_type, Args...>;

  using fom_system_t		= typename lspg_problem_t::fom_system_t;
  using scalar_t		= typename lspg_problem_t::scalar_t;
  using fom_native_state_t	= typename lspg_problem_t::fom_native_state_t;
  using fom_native_residual_t	= typename lspg_problem_t::fom_native_residual_t;
  using fom_state_t		= typename lspg_problem_t::fom_state_t;

  using decoder_t		= typename lspg_problem_t::decoder_t;
  using fom_state_reconstr_t	= typename lspg_problem_t::fom_state_reconstr_t;
  using fom_states_manager_t	= typename lspg_problem_t::fom_states_manager_t;
  using ud_ops_t		= typename lspg_problem_t::ud_ops_t;

  using lspg_state_t		= typename lspg_problem_t::lspg_state_t;
  using lspg_residual_t		= typename lspg_problem_t::lspg_residual_t;
  using lspg_matrix_t		= typename lspg_problem_t::lspg_matrix_t;

  using lspg_residual_policy_t	= typename lspg_problem_t::lspg_residual_policy_t;
  using lspg_jacobian_policy_t	= typename lspg_problem_t::lspg_jacobian_policy_t;

  using lspg_stepper_t		= typename lspg_problem_t::lspg_stepper_t;

private:
  ::pressio::ode::types::step_t step0_;
  scalar_t			t0_;
  scalar_t			dt0_;
  const fom_state_t		fomStateReference_;
  const fom_state_reconstr_t	fomStateReconstructor_;
  fom_states_manager_t		fomStatesMngr_;

  lspg_residual_policy_t	residualPolicy_;
  lspg_jacobian_policy_t	jacobianPolicy_;
  lspg_stepper_t		stepperObj_;

public:
  lspg_stepper_t & getStepperRef(){
    return stepperObj_;
  }

  const fom_native_state_t & viewCurrentFomState() const{
    return *fomStatesMngr_.getCRefToCurrentFomState().data();
  }

  const fom_state_reconstr_t & getFomStateReconstructorCRef() const{
    return fomStateReconstructor_;
  }

public:

  template <
  typename _ud_ops_t = ud_ops_t,
  mpl::enable_if_t< std::is_void<_ud_ops_t>::value, int > = 0
  >
  ProblemDiscreteTimeApi(const fom_system_t & fomSystemObj,
			 const fom_native_state_t & fomStateReferenceNative,
			 decoder_t	 & decoder,
			 lspg_state_t & yROM,
			 scalar_t	t0)
    : step0_{},
      t0_{t0},
      dt0_{},
      fomStateReference_(fomStateReferenceNative),
      fomStateReconstructor_(fomStateReference_, decoder),
      fomStatesMngr_(fomStateReconstructor_, fomStateReference_),
      // construct policies
      residualPolicy_(fomStatesMngr_),
      jacobianPolicy_(fomStatesMngr_, decoder),
      // construct stepper
      stepperObj_(yROM, fomSystemObj, residualPolicy_, jacobianPolicy_)
  {}

  template <
    typename _ud_ops_t = ud_ops_t,
    mpl::enable_if_t< !std::is_void<_ud_ops_t>::value, int > = 0
    >
  ProblemDiscreteTimeApi(const fom_system_t & fomSystemObj,
			 const fom_native_state_t & fomStateReferenceNative,
			 decoder_t & decoder,
			 lspg_state_t & yROM,
			 scalar_t t0,
			 const _ud_ops_t & udOps)
    : step0_{},
      t0_{t0},
      dt0_{},
      fomStateReference_(fomStateReferenceNative),
      fomStateReconstructor_(fomStateReference_, decoder, udOps),
      fomStatesMngr_(fomStateReconstructor_, &udOps, fomStateReference_),
      // construct policies
      residualPolicy_(fomStatesMngr_),
      jacobianPolicy_(fomStatesMngr_, decoder),
      // construct stepper
      stepperObj_(yROM, fomSystemObj, residualPolicy_, jacobianPolicy_)
  {}

};

}}}}}//end namespace pressio::rom::lspg::unsteady::impl
#endif  // ROM_LSPG_IMPL_UNSTEADY_DISCRETE_TIME_API_ROM_LSPG_UNSTEADY_PROBLEM_DISCRETE_TIME_API_HPP_
