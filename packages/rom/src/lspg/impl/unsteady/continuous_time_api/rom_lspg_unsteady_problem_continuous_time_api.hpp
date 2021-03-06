/*
//@HEADER
// ************************************************************************
//
// rom_lspg_unsteady_problem_continuous_time_api.hpp
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

#ifndef ROM_LSPG_IMPL_UNSTEADY_CONTINUOUS_TIME_API_ROM_LSPG_UNSTEADY_PROBLEM_CONTINUOUS_TIME_API_HPP_
#define ROM_LSPG_IMPL_UNSTEADY_CONTINUOUS_TIME_API_ROM_LSPG_UNSTEADY_PROBLEM_CONTINUOUS_TIME_API_HPP_

#include "./discrete_time_functions/rom_lspg_time_discrete_residual.hpp"
#include "./discrete_time_functions/rom_lspg_time_discrete_jacobian.hpp"

#include "./policies/rom_lspg_unsteady_residual_policy_continuous_time_api.hpp"
#include "./policies/rom_lspg_unsteady_jacobian_policy_continuous_time_api.hpp"

#include "./traits/rom_lspg_unsteady_common_traits_continuous_time_api.hpp"
#include "./traits/rom_lspg_unsteady_default_problem_traits_continuous_time_api.hpp"
#include "./traits/rom_lspg_unsteady_preconditioned_problem_traits_continuous_time_api.hpp"
#include "./traits/rom_lspg_unsteady_masked_problem_traits_continuous_time_api.hpp"


namespace pressio{ namespace rom{ namespace lspg{ namespace impl{ namespace unsteady{

template <
  template <class, class, class, class ...> class lspg_type,
  typename stepper_tag,
  typename fom_system_type,
  typename lspg_state_type,
  typename ...Args
  >
struct ProblemContinuousTimeApi
{
public:
  // define the type holding types for the problem
  using lspg_problem_t = lspg_type<stepper_tag, fom_system_type, lspg_state_type, Args...>;

  using fom_system_t		= typename lspg_problem_t::fom_system_t;
  using scalar_t		= typename lspg_problem_t::scalar_t;
  using fom_native_state_t	= typename lspg_problem_t::fom_native_state_t;
  using fom_state_t		= typename lspg_problem_t::fom_state_t;
  using fom_velocity_t		= typename lspg_problem_t::fom_velocity_t;
  using lspg_state_t		= typename lspg_problem_t::lspg_state_t;
  using decoder_t		= typename lspg_problem_t::decoder_t;
  using fom_state_reconstr_t	= typename lspg_problem_t::fom_state_reconstr_t;
  using fom_states_manager_t	= typename lspg_problem_t::fom_states_manager_t;
  using ud_ops_t		= typename lspg_problem_t::ud_ops_t;
  using lspg_matrix_t		= typename lspg_problem_t::lspg_matrix_t;
  using lspg_residual_policy_t	= typename lspg_problem_t::lspg_residual_policy_t;
  using lspg_jacobian_policy_t	= typename lspg_problem_t::lspg_jacobian_policy_t;
  using aux_stepper_t		= typename lspg_problem_t::aux_stepper_t;
  using lspg_stepper_t		= typename lspg_problem_t::lspg_stepper_t;

private:
  const fom_state_t		fomStateReference_;
  const fom_velocity_t		fomVelocityRef_;
  const fom_state_reconstr_t	fomStateReconstructor_;
  fom_states_manager_t		fomStatesMngr_;
  lspg_matrix_t			jPhiMatrix_;
  lspg_residual_policy_t	residualPolicy_;
  lspg_jacobian_policy_t	jacobianPolicy_;

  /* here we use conditional type if auxiliary stepper is non-void,
   * otherwise we set it to a dummy type and we dont construct it */
  typename std::conditional<
    std::is_void<aux_stepper_t>::value,
    ::pressio::utils::impl::empty, aux_stepper_t
    >::type auxStepperObj_ = {};

  // actual stepper object
  lspg_stepper_t  stepperObj_;

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
  ProblemContinuousTimeApi() = delete;
  ~ProblemContinuousTimeApi() = default;

  /* specialize for:
   * - the fom_system_t is regular c++
   * - aux stepper is NOT needed (e.g. for BDF1)
   * - ud_ops_t == void
   */
  template <
    typename _fom_system_t = fom_system_t,
    typename _aux_stepper_t = aux_stepper_t,
    typename _ud_ops_t = ud_ops_t,
    ::pressio::mpl::enable_if_t<
      !::pressio::ops::predicates::is_object_pybind<_fom_system_t>::value and
      std::is_void<_aux_stepper_t>::value and
      std::is_void<_ud_ops_t>::value,
      int > = 0
  >
  ProblemContinuousTimeApi(const _fom_system_t	& fomSystemObj,
			   const fom_native_state_t & fomStateReferenceNative,
			   const decoder_t & decoder,
			   lspg_state_t	 & yROM,
			   scalar_t	 t0)
    : fomStateReference_(fomStateReferenceNative),
      fomVelocityRef_(fomSystemObj.createVelocity()),
      fomStateReconstructor_(fomStateReference_, decoder),
      fomStatesMngr_(fomStateReconstructor_, fomStateReference_),
      jPhiMatrix_(fomSystemObj.createApplyJacobianResult( *decoder.getReferenceToJacobian().data() )),
      // here we pass a fom velocity object to the residual policy to
      // use it to initialize the residual data
      // since the lspg residual is of same type and size of the fom velocity
      // (this is true w and w/o hyperreduction)
      residualPolicy_(fomVelocityRef_, fomStatesMngr_),
      jacobianPolicy_(fomStatesMngr_, jPhiMatrix_, decoder),
      auxStepperObj_{},
      stepperObj_(yROM, fomSystemObj, residualPolicy_, jacobianPolicy_)
  {}


  /* specialize for:
   * - the fom_system_t is regular c++
   * - aux stepper is NOT needed (e.g. for BDF1)
   * - ud_ops_t == non-void
   */
  template <
    typename _fom_system_t = fom_system_t,
    typename _aux_stepper_t = aux_stepper_t,
    typename _ud_ops_t = ud_ops_t,
    ::pressio::mpl::enable_if_t<
      !::pressio::ops::predicates::is_object_pybind<_fom_system_t>::value and
      std::is_void<_aux_stepper_t>::value and
      !std::is_void<_ud_ops_t>::value,
      int > = 0
  >
  ProblemContinuousTimeApi(const _fom_system_t & fomSystemObj,
			   const fom_native_state_t & fomStateReferenceNative,
			   const decoder_t & decoder,
			   lspg_state_t & yROM,
			   scalar_t t0,
			   const _ud_ops_t & udOps)
    : fomStateReference_(fomStateReferenceNative),
      fomVelocityRef_(fomSystemObj.createVelocity()),
      fomStateReconstructor_(fomStateReference_, decoder, udOps),
      fomStatesMngr_(fomStateReconstructor_, fomStateReference_),
      jPhiMatrix_(fomSystemObj.createApplyJacobianResult(*decoder.getReferenceToJacobian().data())),
      // here we pass a fom velocity object to the residual policy to
      // use it to initialize the residual data
      // since the lspg residual is of same type and size of the fom velocity
      // (this is true w and w/o hyperreduction)
      residualPolicy_(fomVelocityRef_, fomStatesMngr_, udOps),
      jacobianPolicy_(fomStatesMngr_, jPhiMatrix_, decoder, udOps),
      auxStepperObj_{},
      stepperObj_(yROM, fomSystemObj, residualPolicy_, jacobianPolicy_)
  {}


  /* specialize for:
   * - the fom_system_t is regular c++
   * - aux stepper is needed (e.g. for BDF2)
   * - ud_ops_t == void
   */
  template <
    typename _fom_system_t = fom_system_t,
    typename _aux_stepper_t = aux_stepper_t,
    typename _ud_ops_t = ud_ops_t,
    ::pressio::mpl::enable_if_t<
      !::pressio::ops::predicates::is_object_pybind<_fom_system_t>::value and
      !std::is_void<_aux_stepper_t>::value and
      std::is_void<_ud_ops_t>::value,
      int > = 0
    >
  ProblemContinuousTimeApi(const _fom_system_t & fomSystemObj,
			   const fom_native_state_t & fomStateReferenceNative,
			   const decoder_t & decoder,
			   lspg_state_t	& yROM,
			   scalar_t t0)
    : fomStateReference_(fomStateReferenceNative),
      fomVelocityRef_(fomSystemObj.createVelocity()),
      fomStateReconstructor_(fomStateReference_, decoder),
      fomStatesMngr_(fomStateReconstructor_, fomStateReference_),
      jPhiMatrix_(fomSystemObj.createApplyJacobianResult(*decoder.getReferenceToJacobian().data())),
      // here we pass a fom velocity object to the residual policy to
      // use it to initialize the residual data
      // since the lspg residual is of same type and size of the fom velocity
      // (this is true w and w/o hyperreduction)
      residualPolicy_(fomVelocityRef_, fomStatesMngr_),
      jacobianPolicy_(fomStatesMngr_, jPhiMatrix_, decoder),
      auxStepperObj_(yROM, fomSystemObj, residualPolicy_, jacobianPolicy_),
      stepperObj_(yROM, fomSystemObj, residualPolicy_, jacobianPolicy_, auxStepperObj_)
  {}


#ifdef PRESSIO_ENABLE_TPL_PYBIND11
  /*- the fom_system_t is a pybind11::object
   * - aux stepper is NOT needed (e.g. for BDF1)
   * - ud_ops_t == void
   */
  template <
    typename _fom_system_t = fom_system_t,
    typename _lspg_state_t = lspg_state_t,
    typename _aux_stepper_t = aux_stepper_t,
    typename _ud_ops_t = ud_ops_t,
    ::pressio::mpl::enable_if_t<
      ::pressio::ops::predicates::is_object_pybind<_fom_system_t>::value and
      ::pressio::containers::predicates::is_vector_wrapper_pybind<_lspg_state_t>::value and
      std::is_void<_aux_stepper_t>::value and
      std::is_void<_ud_ops_t>::value,
      int > = 0
  >
  ProblemContinuousTimeApi(const _fom_system_t & fomSystemObj,
			   const fom_native_state_t fomStateReferenceIn,
			   const decoder_t & decoder,
			   typename ::pressio::containers::details::traits<_lspg_state_t>::wrapped_t & yROM,
			   scalar_t t0)
    : fomStateReference_(fomStateReferenceIn),
      fomVelocityRef_( fomSystemObj.attr("velocity")(fomStateReferenceIn, t0) ),
      fomStateReconstructor_(fomStateReference_, decoder),
      fomStatesMngr_(fomStateReconstructor_, fomStateReference_),
      jPhiMatrix_( fomSystemObj.attr("applyJacobian")(fomStateReferenceIn, *decoder.getReferenceToJacobian().data(), t0)),
      residualPolicy_(fomVelocityRef_, fomStatesMngr_),
      jacobianPolicy_(fomStatesMngr_, jPhiMatrix_, decoder),
      stepperObj_(_lspg_state_t(yROM), fomSystemObj, residualPolicy_, jacobianPolicy_)
  {}
#endif

};

}}}}}
#endif  // ROM_LSPG_IMPL_UNSTEADY_CONTINUOUS_TIME_API_ROM_LSPG_UNSTEADY_PROBLEM_CONTINUOUS_TIME_API_HPP_
