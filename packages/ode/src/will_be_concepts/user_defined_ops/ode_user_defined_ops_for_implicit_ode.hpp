
#ifndef ODE_IS_VALID_USER_DEFINED_OPS_IMPLICIT_ODE_HPP_
#define ODE_IS_VALID_USER_DEFINED_OPS_IMPLICIT_ODE_HPP_

namespace pressio{ namespace ode{ namespace concepts {

template<
  typename T,
  typename scalar_t,
  typename state_t,
  typename residual_t,
  typename jacobian_t,
  typename enable = void
  >
struct user_defined_ops_for_implicit_ode : std::false_type{};

template<
  typename T,
  typename scalar_t,
  typename state_t,
  typename residual_t,
  typename jacobian_t
  >
struct user_defined_ops_for_implicit_ode<
  T, scalar_t, state_t, residual_t, jacobian_t,
  mpl::enable_if_t<
    is_valid_user_defined_ops_for_implicit_euler<
      T, scalar_t, state_t, residual_t, jacobian_t
      >::value 
    and
    is_valid_user_defined_ops_for_implicit_bdf2<
      T, scalar_t, state_t, residual_t, jacobian_t
      >::value
    >
  > : std::true_type{};

}}} // namespace pressio::ode::concepts
#endif