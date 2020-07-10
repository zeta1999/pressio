
#ifndef ODE_META_HAS_TIME_DISCRETE_jacobian_TYPEDEF_HPP_
#define ODE_META_HAS_TIME_DISCRETE_jacobian_TYPEDEF_HPP_

namespace pressio{ namespace ode{ namespace predicates {

template <typename T, typename enable = void>
struct has_discrete_time_jacobian_typedef : std::false_type{};

template <typename T>
struct has_discrete_time_jacobian_typedef<
  T,
  mpl::enable_if_t<
    !std::is_void<
      typename T::discrete_time_jacobian_type
      >::value
    >
  > : std::true_type{};

}}}//end namespace pressio::ode::predicates
#endif