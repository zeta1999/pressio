
#ifdef HAVE_TRILINOS
#ifndef CONTAINERS_IS_VECTOR_WRAPPER_TPETRA_HPP_
#define CONTAINERS_IS_VECTOR_WRAPPER_TPETRA_HPP_

#include "../containers_vector_traits.hpp"

namespace rompp{ namespace containers{ namespace meta {

template <typename T, typename enable = void>
struct is_vector_wrapper_tpetra : std::false_type {};

template <typename T>
struct is_vector_wrapper_tpetra<
  T, ::rompp::mpl::enable_if_t<
       containers::details::traits<T>::is_vector &&
       containers::details::traits<T>::wrapped_vector_identifier==
       containers::details::WrappedVectorIdentifier::Tpetra
       >
  > : std::true_type{};

}}}//end namespace rompp::containers::meta
#endif
#endif