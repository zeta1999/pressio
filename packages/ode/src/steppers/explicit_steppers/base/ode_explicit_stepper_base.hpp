
#ifndef ODE_EXPLICIT_STEPPER_BASE_HPP_
#define ODE_EXPLICIT_STEPPER_BASE_HPP_

#include "ode_ConfigDefs.hpp"
#include "../ode_explicit_stepper_traits.hpp"
#include "../../../meta/ode_meta.hpp"
#include "../../../meta/ode_meta_explicit.hpp"
#include "../../../policies/meta/ode_explicit_policies_meta.hpp"

namespace ode{

template<typename stepper_type>
class explicitStepperBase
{
private:
  using step_traits = ode::details::traits<stepper_type>;

  using state_t = typename step_traits::state_t; 
  using residual_t = typename step_traits::residual_t;
  using model_t = typename step_traits::model_t;
  using time_t = typename step_traits::time_t;
  using residual_policy_t = typename step_traits::residual_policy_t;
  using sizer_t = typename step_traits::sizer_t;

  using order_t = typename step_traits::order_t;
  static constexpr order_t order_value = step_traits::order_value;

  //do checking here that things are as supposed
  static_assert( meta::isLegitimateExplicitStateType<state_t>::value,
		 "OOPS: STATE_TYPE IN SELECTED EXPLICIT STEPPER IS NOT VALID");
  static_assert( meta::isLegitimateExplicitResidualType<residual_t>::value,
		 "OOPS: RESIDUAL_TYPE IN SELECTED EXPLICIT STEPPER IS NOT VALID");
  static_assert( meta::isLegitimateTimeType<time_t>::value,
		 "OOPS: TIME_TYPE IN SELECTED EXPLICIT STEPPER IS NOT VALID");
  static_assert( meta::isLegitimateExplicitResidualPolicy<residual_policy_t>::value,
      "RESIDUAL_POLICY NOT ADMISSIBLE, MAYBE NOT A CHILD OF EXPLICIT POLICY BASE");

public:
  order_t order() const{
    return order_value;
  }
  template<typename step_t>
  void doStep(state_t &inout, time_t t, time_t dt, step_t step){
    this->stepper()->doStepImpl(inout, t, dt, step);
  }

private:    
  explicitStepperBase(model_t & model,
		      residual_policy_t & res_policy_obj)
    : model_(&model),
      residual_obj_(&res_policy_obj)
  {}
  explicitStepperBase() = delete;
  ~explicitStepperBase() = default;

private:
  friend stepper_type;

  stepper_type * stepper( ){
    return static_cast< stepper_type* >( this );
  }
  const stepper_type * stepper( ) const{
    return static_cast< const stepper_type* >( this );
  }

protected:
  model_t * model_;
  residual_policy_t * residual_obj_;
  
};//end class
}//end namespace
#endif