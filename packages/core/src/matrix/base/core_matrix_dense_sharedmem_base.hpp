
#ifndef CORE_MATRIX_BASE_MATRIX_DENSE_SHAREDMEM_BASE_HPP_
#define CORE_MATRIX_BASE_MATRIX_DENSE_SHAREDMEM_BASE_HPP_

#include "../core_matrix_traits.hpp"
#include "../../shared_base/core_operators_base.hpp"

namespace core{
    
template<typename derived_type>
class MatrixDenseSharedMemBase
  : private core::details::CrtpBase<
     MatrixDenseSharedMemBase<derived_type>>,
  public Subscripting2DOperatorsBase<
    MatrixDenseSharedMemBase<derived_type>, 
    typename details::traits<derived_type>::scalar_t,
    typename details::traits<derived_type>::ordinal_t>
{

  static_assert( details::traits<derived_type>::isSharedMem==1,
  "OOPS: distributed matrix inheriting from dense sharedMem base!");
  
private:
  using this_t = MatrixDenseSharedMemBase<derived_type>;
  friend derived_type;
  friend core::details::CrtpBase<this_t>;

  using sc_t = typename details::traits<derived_type>::scalar_t;
  using ord_t = typename details::traits<derived_type>::ordinal_t;
  friend Subscripting2DOperatorsBase< this_t, sc_t, ord_t>;
  
  MatrixDenseSharedMemBase() = default;
  ~MatrixDenseSharedMemBase() = default;
  
};//end class
  
} // end namespace core
#endif
