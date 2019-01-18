
#ifndef SOLVERS_IMPL_NORM_HELPER_MIXIN_HPP
#define SOLVERS_IMPL_NORM_HELPER_MIXIN_HPP

#include "../../solvers_ConfigDefs.hpp"
#include "../../../../CORE_OPS"

namespace rompp{ namespace solvers{ namespace iterative{ namespace impl{

template <typename convergence_tag>
struct NormSelectorHelper{
  using norm_t = L2Norm;
};

template <typename norm_type>
struct NormSelectorHelper<
  converged_when::absoluteNormCorrectionBelowTol<norm_type>
  >{
  using norm_t = norm_type;
};
//---------------------------------------------------------


template <typename norm_t>
struct ComputeNormHelper;

template <>
struct ComputeNormHelper<::rompp::solvers::L2Norm>{
  template <typename vec_t, typename scalar_t>
  void operator()(const vec_t & vecIn, scalar_t & result) const{
    result = ::rompp::core::ops::norm2(vecIn);
  }
};

template <>
struct ComputeNormHelper<::rompp::solvers::L1Norm>{
  template <typename vec_t, typename scalar_t>
  void operator()(const vec_t & vecIn, scalar_t & result) const{
    result = ::rompp::core::ops::norm1(vecIn);
  }
};
//---------------------------------------------------------


}}}} //end namespace rompp::solvers::iterative::impl
#endif