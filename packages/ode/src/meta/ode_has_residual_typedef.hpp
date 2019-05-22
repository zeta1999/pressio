
#ifndef ODE_META_HAS_RESIDUAL_TYPEDEF_HPP_
#define ODE_META_HAS_RESIDUAL_TYPEDEF_HPP_

#include <type_traits>

namespace rompp{ namespace ode{ namespace meta {

template <typename T, typename enable = void>
struct has_residual_typedef : std::false_type{};

template <typename T>
struct has_residual_typedef<
  T,
  mpl::enable_if_t<
    !std::is_void<
      typename T::residual_type
      >::value
    >
  > : std::true_type{};

}}}//end namespace rompp::ode::meta
#endif