
#ifndef ODE_IMPLICIT_STEPPER_TRAITS_HPP_
#define ODE_IMPLICIT_STEPPER_TRAITS_HPP_

#include "ode_forward_declarations.hpp"

namespace ode{
namespace details{
  
template<typename state_type,
	 typename residual_type,
	 typename jacobian_type,
	 typename scalar_type,
	 typename model_type,
	 typename time_type,
	 typename sizer_type,
	 typename solver_policy_type,
	 typename residual_policy_type,
	 typename jacobian_policy_type>
struct traits< impl::implicitEulerStepperImpl<state_type,
					      residual_type,
					      jacobian_type,
					      scalar_type,
					      model_type,
					      time_type,
					      sizer_type,
					      solver_policy_type,
					      residual_policy_type,
					      jacobian_policy_type>>
{
  using stepper_t =
    impl::implicitEulerStepperImpl<state_type,
				   residual_type,
				   jacobian_type,
				   scalar_type,
				   model_type,
				   time_type,
				   sizer_type,
				   solver_policy_type,
				   residual_policy_type,
				   jacobian_policy_type>;
  using state_t =  state_type;
  using residual_t = residual_type;
  using jacobian_t =  jacobian_type;
  using scalar_t = scalar_type;    
  using model_t = model_type;
  using time_t = time_type;
  using sizer_t = sizer_type;
  using solver_policy_t = solver_policy_type;
  using residual_policy_t = residual_policy_type;
  using jacobian_policy_t = jacobian_policy_type;

  static constexpr bool advanceIncrement = residual_policy_t::advanceIncrement;
  static_assert(residual_policy_t::advanceIncrement ==
		jacobian_policy_t::advanceIncrement,
		"Residual and jacobian policies BOTH need to advance full state or just increment wrt initial condition. In this case they are not");

  using order_t = unsigned int;
  static constexpr order_t order_value = 1;
  static constexpr order_t steps = 1;
};

  
}//end namespace details
}//end namespace ode

#endif
