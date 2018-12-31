
#ifdef HAVE_ARMADILLO
#ifndef CORE_NATIVE_ARMADILLO_VECTOR_META_HPP_
#define CORE_NATIVE_ARMADILLO_VECTOR_META_HPP_

#include "core_meta_basic.hpp"
#include <armadillo>

namespace rompp{ namespace core{ namespace meta {

template <typename T, typename enable = void>
struct is_armadillo_column_vector : std::false_type {};

template <typename T>
struct is_armadillo_column_vector<T,
	 core::meta::enable_if_t<
	   std::is_same<T,
     	    arma::Col<typename T::elem_type>
			>::value
	   >
      > : std::true_type{};

template <typename T, typename enable = void>
struct is_armadillo_row_vector : std::false_type {};

template <typename T>
struct is_armadillo_row_vector<T,
	 core::meta::enable_if_t<
	   std::is_same<T,
     	    arma::Row<typename T::elem_type>
			>::value
	   >
      > : std::true_type{};

template <typename T, typename enable = void>
struct is_vector_armadillo : std::false_type {};

template <typename T>
struct is_vector_armadillo<T,
	 core::meta::enable_if_t<
	   is_armadillo_row_vector<T>::value or
	   is_armadillo_column_vector<T>::valu
	   >
      > : std::true_type{};

}}}//end namespace rompp::core::meta
#endif
#endif