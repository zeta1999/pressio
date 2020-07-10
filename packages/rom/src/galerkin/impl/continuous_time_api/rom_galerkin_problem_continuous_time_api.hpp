
#ifndef rom_galerkin_problem_continuous_time_api_HPP_
#define rom_galerkin_problem_continuous_time_api_HPP_

#include "rom_galerkin_explicit_velocity_policy.hpp"
#include "./traits/rom_galerkin_common_traits_continuous_time_api.hpp"
#include "./traits/rom_galerkin_default_problem_traits_continuous_time_api.hpp"

namespace pressio{ namespace rom{ namespace galerkin{ namespace impl{

template <
  template <class ...> class galerkin_type,
  typename stepper_tag,
  typename fom_type,
  typename rom_state_type,
  typename ...Args
  >
class ProblemContinuousTimeApi
{

public:
  // define the type holding types for the problem
  using problem_t = galerkin_type<stepper_tag, fom_type, rom_state_type, Args...>;

  using fom_t			= typename problem_t::fom_t;
  using scalar_t		= typename problem_t::scalar_t;
  using fom_native_state_t	= typename problem_t::fom_native_state_t;
  using fom_state_t		= typename problem_t::fom_state_t;
  using fom_velocity_t		= typename problem_t::fom_velocity_t;

  using galerkin_state_t	= typename problem_t::galerkin_state_t;
  using galerkin_native_state_t	= typename problem_t::galerkin_native_state_t;
  using decoder_t		= typename problem_t::decoder_t;
  using fom_state_reconstr_t	= typename problem_t::fom_state_reconstr_t;
  using fom_states_manager_t		= typename problem_t::fom_states_manager_t;
  using ud_ops_t		= typename problem_t::ud_ops_t;

  using residual_policy_t	= typename problem_t::residual_policy_t;
  using stepper_t		= typename problem_t::stepper_t;

private:
  fom_state_t			fomStateReference_;
  fom_state_reconstr_t		fomStateReconstructor_;
  fom_velocity_t		fomVelocityRef_;
  fom_states_manager_t		fomStatesMngr_;
  residual_policy_t		residualPolicy_;
  stepper_t			stepperObj_;

public:
  stepper_t & getStepperRef(){
    return stepperObj_;
  }

  const fom_state_reconstr_t & getFomStateReconstructorCRef() const{
    return fomStateReconstructor_;
  }

public:
  ProblemContinuousTimeApi() = delete;
  ~ProblemContinuousTimeApi() = default;

  /*
   * ud_ops_t = void, C++ types
  */
  template <
    typename _ud_ops_t = ud_ops_t,
    ::pressio::mpl::enable_if_t<
      std::is_void<_ud_ops_t>::value and
      ::pressio::containers::predicates::is_wrapper<galerkin_state_t>::value
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
      and !::pressio::containers::predicates::is_vector_wrapper_pybind<galerkin_state_t>::value
#endif
      , int> = 0
  >
  ProblemContinuousTimeApi(const fom_t   & appObj,
			      const fom_native_state_t & yFomRefNative,
			      const decoder_t	    & decoder,
			      galerkin_state_t	    & yROM,
			      scalar_t		    t0)
    : fomStateReference_(yFomRefNative),
      fomStateReconstructor_(fomStateReference_, decoder),
      fomVelocityRef_(appObj.createVelocity()),
      fomStatesMngr_(fomStateReconstructor_, fomStateReference_),
      residualPolicy_(fomVelocityRef_, fomStatesMngr_, decoder),
      stepperObj_(yROM, appObj, residualPolicy_)
  {}

  /*
   * ud_ops_t != void, C++ types
  */
  template <
    typename _ud_ops_t = ud_ops_t,
    ::pressio::mpl::enable_if_t<
      !std::is_void<_ud_ops_t>::value and
      ::pressio::containers::predicates::is_wrapper<galerkin_state_t>::value
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
      and !::pressio::containers::predicates::is_vector_wrapper_pybind<galerkin_state_t>::value
#endif
      , int> = 0
  >
  ProblemContinuousTimeApi(const fom_t   & appObj,
			      const fom_native_state_t & yFomRefNative,
			      const decoder_t	    & decoder,
			      galerkin_state_t	    & yROM,
			      scalar_t		    t0,
			      const _ud_ops_t & udOps)
    : fomStateReference_(yFomRefNative),
      fomStateReconstructor_(fomStateReference_, decoder, udOps),
      fomVelocityRef_(appObj.createVelocity()),
      fomStatesMngr_(fomStateReconstructor_, &udOps, fomStateReference_),
      residualPolicy_(fomVelocityRef_, fomStatesMngr_, decoder, udOps),
      stepperObj_(yROM, appObj, residualPolicy_)
  {}

// #ifdef PRESSIO_ENABLE_TPL_PYBIND11
//   /*
//    * ud_ops_t == void and state_type is wrapper of pybind11::array
//   */
//   template <
//     typename _ud_ops_t = ud_ops_t,
//     ::pressio::mpl::enable_if_t<
//       std::is_void<_ud_ops_t>::value and
//       ::pressio::containers::predicates::is_vector_wrapper_pybind<galerkin_state_t>::value, 
//       int > = 0
//   >
//   ProblemContinuousTimeApi(const fom_t   & appObj,
// 			      fom_native_state_t	    yFomRefNative,
// 			      const decoder_t	    & decoder,
// 			      galerkin_native_state_t  yROM,
// 			      scalar_t		    t0)
//     : fomStateReference_(yFomRefNative),
//       fomStateReconstructor_(fomStateReference_, decoder),
//       fomVelocityRef_( appObj.attr("velocity")(*fomStateReference_.data(), t0) ),
//       fomStatesMngr_(fomStateReconstructor_, fomStateReference_),
//       residualPolicy_(fomVelocityRef_, fomStatesMngr_, decoder),
//       stepperObj_(galerkin_state_t(yROM), appObj, residualPolicy_)
//   {}
// #endif

};

}}}}//end namespace pressio::rom::galerkin::impl
#endif