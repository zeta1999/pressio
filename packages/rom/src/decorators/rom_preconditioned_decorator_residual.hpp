/*
//@HEADER
// ************************************************************************
//
// rom_preconditioned_decorator_residual.hpp
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

#ifndef ROM_DECORATORS_ROM_PRECONDITIONED_DECORATOR_RESIDUAL_HPP_
#define ROM_DECORATORS_ROM_PRECONDITIONED_DECORATOR_RESIDUAL_HPP_

namespace pressio{ namespace rom{ namespace decorator{

template <typename preconditionable_policy>
class PreconditionedResidualPolicy : public preconditionable_policy
{

  using typename preconditionable_policy::residual_t;
  using preconditionable_policy::fomStatesMngr_;

public:
  PreconditionedResidualPolicy() = delete;

  PreconditionedResidualPolicy(const preconditionable_policy & obj)
    : preconditionable_policy(obj)
  {}

  template <typename ... Args>
  PreconditionedResidualPolicy(Args && ... args)
    : preconditionable_policy(std::forward<Args>(args)...)
  {}

  ~PreconditionedResidualPolicy() = default;

public:
  template <typename fom_system_t>
  residual_t create(const fom_system_t & systemObj) const
  {
    return preconditionable_policy::create(systemObj);
  }

  //-------------------------------
  // unsteady case
  //-------------------------------
  template <
    typename stepper_tag,
    typename prev_states_t,
    typename fom_system_t,
    typename scalar_t,
    typename state_t
    >
  void compute(const state_t & currentState,
	       const prev_states_t & prevStates,
	       const fom_system_t & systemObj,
	       const scalar_t & t,
	       const scalar_t & dt,
	       const ::pressio::ode::types::step_t & step,
	       residual_t & residual,
	       ::pressio::Norm normKind,
	       scalar_t & normValue) const
  {
    preconditionable_policy::template compute<
      stepper_tag>(currentState, prevStates, systemObj, t, dt, step, residual, normKind, normValue);

    const auto & yFom = fomStatesMngr_.getCRefToCurrentFomState();
    systemObj.applyPreconditioner(*yFom.data(), t, *residual.data());
  }

  //-------------------------------
  // steady case
  //-------------------------------
  template <typename lspg_state_t, typename fom_t, typename norm_val_t>
  void compute(const lspg_state_t & currentState,
              residual_t & residual,
              const fom_t & systemObj,
              ::pressio::Norm normKind,
              norm_val_t & normValue) const
  {
    preconditionable_policy::compute(currentState, residual, systemObj, normKind, normValue);

    const auto & yFom = fomStatesMngr_.getCRefToCurrentFomState();
    systemObj.applyPreconditioner(*yFom.data(), *residual.data());
  }

};//end class

}}} //end namespace pressio::rom::decorator
#endif  // ROM_DECORATORS_ROM_PRECONDITIONED_DECORATOR_RESIDUAL_HPP_
