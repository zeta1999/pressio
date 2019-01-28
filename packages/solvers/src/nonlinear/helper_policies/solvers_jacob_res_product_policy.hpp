
#ifndef SOLVERS_IMPL_JACOBIAN_RESIDUAL_PRODUCT_POLICY_HPP
#define SOLVERS_IMPL_JACOBIAN_RESIDUAL_PRODUCT_POLICY_HPP

#include "../../solvers_ConfigDefs.hpp"
#include "../../../../CORE_OPS"

namespace rompp{ namespace solvers{ namespace iterative{ namespace impl{

template<typename J_t, typename enable = void>
struct JacobianTranspResProdHelper;


// when J is a matrix wrapper, then J^T*R
// is computed doing regular mat-vec product
template<typename J_t>
struct JacobianTranspResProdHelper<
  J_t,
  core::meta::enable_if_t<
    core::meta::is_core_matrix_wrapper<J_t>::value
    >>{

  template <typename resid_t, typename result_t>
  void operator()(J_t & J, resid_t & R, result_t & result) const
  {
    constexpr bool transposeJ = true;
    ::rompp::core::ops::product<J_t, resid_t, result_t,
				transposeJ>(J, R, result);
  }
};


// when J is multivector wrapper, then J^T*R
// can be computed doing the DOT of J*R
template<typename J_t>
struct JacobianTranspResProdHelper<
  J_t,
  core::meta::enable_if_t<
    core::meta::is_core_multi_vector_wrapper<J_t>::value
    >>{

  template <typename resid_t, typename result_t>
  void operator()(J_t & J, resid_t & R, result_t & result) const{
    ::rompp::core::ops::dot(J, R, result);
  }
};


}}}} //end namespace rompp::solvers::iterative::impl
#endif