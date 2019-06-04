
#ifndef ODE_STEPPERS_IMPLICIT_STEPPERS_BASE_IMPLICIT_STEPPER_BASE_HPP_
#define ODE_STEPPERS_IMPLICIT_STEPPERS_BASE_IMPLICIT_STEPPER_BASE_HPP_

#include "ode_implicit_stepper_traits.hpp"
#include "../policies/meta/ode_is_implicit_jacobian_standard_policy.hpp"
#include "../policies/meta/ode_is_implicit_residual_standard_policy.hpp"
#include "../policies/meta/ode_is_legitimate_implicit_jacobian_policy.hpp"
#include "../policies/meta/ode_is_legitimate_implicit_residual_policy.hpp"
#include "../../ode_storage.hpp"
#include "../../ode_aux_data.hpp"

namespace rompp{ namespace ode{

template<typename concrete_stepper_type, int nAuxStates>
class ImplicitStepperBase
  : private core::details::CrtpBase<ImplicitStepperBase<concrete_stepper_type, nAuxStates>>
{
  using traits		  = typename details::traits<concrete_stepper_type>;
  using sc_t		  = typename traits::scalar_t;
  using state_t		  = typename traits::state_t;
  using residual_t	  = typename traits::residual_t;
  using jacobian_t	  = typename traits::jacobian_t;
  using standard_res_policy_t = typename traits::standard_res_policy_t;
  using standard_jac_policy_t = typename traits::standard_jac_policy_t;
  using residual_pol_t = typename traits::residual_policy_t;
  using jacobian_pol_t = typename traits::jacobian_policy_t;
  using model_t		  = typename traits::model_t;

  //do checking here that things are as supposed
  static_assert( meta::is_legitimate_implicit_state_type<state_t>::value,
       "OOPS: STATE_TYPE IN SELECTED IMPLICIT STEPPER IS NOT VALID");
  static_assert( meta::is_legitimate_implicit_residual_type<residual_t>::value,
       "OOPS: RESIDUAL_TYPE IN SELECTED IMPLICIT STEPPER IS NOT VALID");
  static_assert( meta::is_legitimate_jacobian_type<jacobian_t>::value,
       "OOPS: JACOBIAN_TYPE IN SELECTED IMPLICIT STEPPER IS NOT VALID");

protected:
  impl::OdeStorage<state_t, residual_t, nAuxStates> odeStorage_;
  impl::ImpOdeAuxData<model_t, sc_t> auxData_;

  typename std::conditional<
    mpl::is_same<standard_res_policy_t, residual_pol_t>::value,
    const residual_pol_t,
    const residual_pol_t &
    >::type residual_obj_;

  typename std::conditional<
    mpl::is_same<standard_jac_policy_t, jacobian_pol_t>::value,
    const jacobian_pol_t,
    const jacobian_pol_t &
    >::type jacobian_obj_;

public:
  decltype(traits::order_value) order() const{
    return traits::order_value;
  }

  void residual(const state_t & y,
		residual_t & R) const{
    this->residual_obj_.template operator()<
      traits::enum_id,
      traits::steps
      >(y, R, odeStorage_.auxStates_, auxData_.model_, auxData_.t_, auxData_.dt_);
  }

  void jacobian(const state_t & y,
		jacobian_t & J) const{
    this->jacobian_obj_.template operator()<
      traits::enum_id
      >(y, J, auxData_.model_, auxData_.t_, auxData_.dt_);
  }

  residual_t residual(const state_t & y) const{
    return this->residual_obj_.template operator()<
      traits::enum_id,
      traits::steps
      >(y, odeStorage_.auxStates_, auxData_.model_, auxData_.t_, auxData_.dt_);
  }

  jacobian_t jacobian(const state_t & y) const{
    return this->jacobian_obj_.template operator()<
      traits::enum_id
      >(y, auxData_.model_, auxData_.t_, auxData_.dt_);
  }


private:
  ImplicitStepperBase(const state_t & y0,
		      const model_t & model,
		      const residual_pol_t & resPolicyObj,
		      const jacobian_pol_t & jacPolicyObj)
    : odeStorage_{y0},
      auxData_{model},
      residual_obj_{resPolicyObj},
      jacobian_obj_{jacPolicyObj}{}

  // cstr for standard residual and jacob policies
  template <
    typename T1 = standard_res_policy_t,
    typename T2 = standard_jac_policy_t,
    ::rompp::mpl::enable_if_t<
      mpl::is_same<T1, residual_pol_t>::value and
      mpl::is_same<T2, jacobian_pol_t>::value
      > * = nullptr
    >
  ImplicitStepperBase(const state_t & y0,
  		      const model_t & model)
    : odeStorage_{y0},
      auxData_{model},
      residual_obj_{},
      jacobian_obj_{}{}

  // cstr for standard jacob policies
  template <
    typename T2 = standard_jac_policy_t,
    ::rompp::mpl::enable_if_t<
      mpl::is_same<T2, jacobian_pol_t>::value
      > * = nullptr
    >
  ImplicitStepperBase(const state_t & y0,
  		      const model_t & model,
  		      const residual_pol_t & resPolicyObj)
    : odeStorage_{y0},
      auxData_{model},
      residual_obj_{resPolicyObj},
      jacobian_obj_{}{}

  ImplicitStepperBase() = delete;
  ~ImplicitStepperBase() = default;

  /* workaround for nvcc issue with templates, see https://devtalk.nvidia.com/default/topic/1037721/nvcc-compilation-error-with-template-parameter-as-a-friend-within-a-namespace/ */
  template<typename DummyType> struct dummy{using type = DummyType;};
  friend typename dummy<concrete_stepper_type>::type;

  friend core::details::CrtpBase<ImplicitStepperBase<concrete_stepper_type, nAuxStates>>;

};//end class

}}//end namespace rompp::ode
#endif
