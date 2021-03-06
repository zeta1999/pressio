/*
//@HEADER
// ************************************************************************
//
// containers_matrix_subspan_expression.hpp
//                     		  Pressio
//                             Copyright 2019
//    National Technology & Engineering Solutions of Sandia, LLC (NTESS)
//
// Under the terms of Contract DE-NA0003525 with NTESS, the
// U.S. Government retains certain rights in this software.
//
// Pressio is licensed under BSD-3-Clause terms of use:
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its
// contributors may be used to endorse or promote products derived
// from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Francesco Rizzi (fnrizzi@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef CONTAINERS_EXPRESSIONS_SUBSPAN_CONTAINERS_MATRIX_SUBSPAN_EXPRESSION_HPP_
#define CONTAINERS_EXPRESSIONS_SUBSPAN_CONTAINERS_MATRIX_SUBSPAN_EXPRESSION_HPP_

namespace pressio{ namespace containers{ namespace expressions{

template <typename matrix_t>
struct SubspanExpr<
  matrix_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::predicates::is_dense_matrix_wrapper_eigen<matrix_t>::value
    >
  >
  : public MatrixSharedMemBase< SubspanExpr<matrix_t> >
{

  using this_t = SubspanExpr<matrix_t>;
  using mytraits = typename details::traits<this_t>;
  using sc_t = typename mytraits::scalar_t;
  using ord_t = typename mytraits::ordinal_t;
  using size_t = typename mytraits::size_t;

  using ref_t = typename mytraits::reference_t;
  using const_ref_t = typename mytraits::const_reference_t;

  using native_expr_t = typename mytraits::native_expr_t;
  using data_return_t = typename mytraits::data_return_t;
  using const_data_return_t = typename mytraits::const_data_return_t;

  using pair_t = std::pair<std::size_t, std::size_t>;

private:
  matrix_t & matObj_;
  ord_t rowStart_;
  ord_t colStart_;
  ord_t endRow_;
  ord_t endCol_;
  ord_t numRows_ = {};
  ord_t numCols_ = {};
  native_expr_t nativeExprObj_;

public:
  SubspanExpr() = delete;
  ~SubspanExpr() = default;
  SubspanExpr(const SubspanExpr & other) = default;
  SubspanExpr(SubspanExpr && other) = default;
  SubspanExpr & operator=(const SubspanExpr & other) = default;
  SubspanExpr & operator=(SubspanExpr && other) = default;

  SubspanExpr(matrix_t & matObjIn,
	      const pair_t rowRangeIn,
	      const pair_t colRangeIn)
    : matObj_(matObjIn),
    rowStart_(std::get<0>(rowRangeIn)),
    colStart_(std::get<0>(colRangeIn)),
    endRow_(std::get<1>(rowRangeIn)-1),
    endCol_(std::get<1>(colRangeIn)-1),
    numRows_(endRow_ - rowStart_ + 1),
    numCols_(endCol_ - colStart_ + 1),
    nativeExprObj_(matObj_.data()->block(rowStart_, colStart_, numRows_, numCols_))
  {
    assert( rowStart_ >= 0 and rowStart_ < matObjIn.extent(0) );
    assert( (int)std::get<1>(rowRangeIn) <= matObjIn.extent(0) );
    assert( colStart_ >= 0 and colStart_ < matObjIn.extent(1) );
    assert( (int)std::get<1>(colRangeIn) <= matObjIn.extent(1) );

    // here the ranges are exclusive of the last index (like Kokkos and Python)
    // so the indices of the last row and col included are:
    assert(endRow_ >= rowStart_);
    assert(endCol_ >= colStart_);
  }

  size_t extent(size_t i) const{
    return (i==0) ? numRows_ : numCols_;
  }

  const_data_return_t data() const{
    return &nativeExprObj_;
  }

  data_return_t data(){
    return &nativeExprObj_;
  }

  ref_t operator()(const ord_t & i, const ord_t & j)
  {
    assert(i < numRows_);
    assert(j < numCols_);
    return nativeExprObj_(i, j);
  }

  const_ref_t operator()(const ord_t & i, const ord_t & j) const
  {
    assert(i < numRows_);
    assert(j < numCols_);
    return nativeExprObj_(i, j);
  }
};



#ifdef PRESSIO_ENABLE_TPL_KOKKOS
template <typename matrix_t>
struct SubspanExpr<
  matrix_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::predicates::is_matrix_wrapper_kokkos<matrix_t>::value
    >
  >
  : public MatrixSharedMemBase< SubspanExpr<matrix_t> >
{

  using this_t = SubspanExpr<matrix_t>;
  using mytraits = typename details::traits<this_t>;
  using sc_t = typename mytraits::scalar_t;
  using ord_t = typename mytraits::ordinal_t;
  using size_t = typename mytraits::size_t;
  using pair_t = typename mytraits::pair_t;

  using ref_t = typename mytraits::reference_t;
  using const_ref_t = typename mytraits::const_reference_t;

  using native_expr_t = typename mytraits::native_expr_t;
  using data_return_t = typename mytraits::data_return_t;
  using const_data_return_t = typename mytraits::const_data_return_t;

private:
  matrix_t & matObj_;
  std::size_t rowStart_;
  std::size_t colStart_;
  std::size_t endRow_;
  std::size_t endCol_;
  std::size_t numRows_ = {};
  std::size_t numCols_ = {};
  native_expr_t nativeExprObj_;

public:
  SubspanExpr() = delete;
  ~SubspanExpr() = default;
  SubspanExpr(const SubspanExpr & other) = default;
  SubspanExpr(SubspanExpr && other) = default;
  SubspanExpr & operator=(const SubspanExpr & other) = default;
  SubspanExpr & operator=(SubspanExpr && other) = default;

  SubspanExpr(matrix_t & matObjIn,
	      const pair_t rowRangeIn,
	      const pair_t colRangeIn)
    : matObj_(matObjIn),
    rowStart_(std::get<0>(rowRangeIn)),
    colStart_(std::get<0>(colRangeIn)),
    endRow_(std::get<1>(rowRangeIn)-1),
    endCol_(std::get<1>(colRangeIn)-1),
    numRows_(endRow_ - rowStart_ + 1),
    numCols_(endCol_ - colStart_ + 1),
    nativeExprObj_(Kokkos::subview(*matObj_.data(),
                   std::make_pair(rowStart_, rowStart_+numRows_),
                   std::make_pair(colStart_, colStart_+numCols_)))
  {
    assert( rowStart_ >= 0 and rowStart_ < matObjIn.extent(0) );
    assert( std::get<1>(rowRangeIn) <= matObjIn.extent(0) );
    assert( colStart_ >= 0 and colStart_ < matObjIn.extent(1) );
    assert( std::get<1>(colRangeIn) <= matObjIn.extent(1) );

    // here the ranges are exclusive of the last index (like Kokkos and Python)
    // so the indices of the last row and col included are:
    assert(endRow_ >= rowStart_);
    assert(endCol_ >= colStart_);
  }

  size_t extent(size_t i) const{
    return (i==0) ? numRows_ : numCols_;
  }

  ref_t operator()(const size_t & i, const size_t & j){
    assert(i < numRows_);
    assert(j < numCols_);
    return nativeExprObj_(i, j);
  }

  const_ref_t const & operator()(const size_t & i, const size_t & j) const{
    assert(i < numRows_);
    assert(j < numCols_);
    return nativeExprObj_(i, j);
  }

  const_data_return_t data() const{
    return &nativeExprObj_;
  }

  data_return_t data(){
    return &nativeExprObj_;
  }
};
#endif

}}} //end namespace pressio::containers::expressions
#endif  // CONTAINERS_EXPRESSIONS_SUBSPAN_CONTAINERS_MATRIX_SUBSPAN_EXPRESSION_HPP_
