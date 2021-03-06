/*
//@HEADER
// ************************************************************************
//
// qr_epetra_multi_vector_tsqr_impl.hpp
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

#ifndef QR_IMPL_EPETRA_QR_EPETRA_MULTI_VECTOR_TSQR_IMPL_HPP_
#define QR_IMPL_EPETRA_QR_EPETRA_MULTI_VECTOR_TSQR_IMPL_HPP_

#include "Epetra_TsqrAdaptor.hpp"

namespace pressio{ namespace qr{ namespace impl{

template<typename matrix_t, typename R_t, typename MV_t, template<typename...> class Q_type>
class EpetraMVTSQR<matrix_t, R_t, MV_t, Q_type, void>
{

  using int_t	     = int;
  using sc_t	     = typename containers::details::traits<matrix_t>::scalar_t;
  using serden_mat_t = Teuchos::SerialDenseMatrix<int_t, sc_t>;
  using trcp_mat     = Teuchos::RCP<serden_mat_t>;
  using Q_t	     = Q_type<MV_t>;
  using tsqr_adaptor_type = Epetra::TsqrAdaptor;

public:
  EpetraMVTSQR() = default;
  ~EpetraMVTSQR() = default;

  void computeThinOutOfPlace(matrix_t & A) {
    auto nVecs = A.numVectors();
    auto & ArowMap = A.data()->Map();
    createQIfNeeded(ArowMap, nVecs);
    createLocalRIfNeeded(nVecs);
    tsqrAdaptor_.factorExplicit(*A.data(),
				*Qmat_->data(),
				*localR_.get(),
				false);
    //::pressio::utils::io::print_stdout(*localR_.get());

// #ifdef PRESSIO_ENABLE_DEBUG_PRINT
//     int myrank{};
//     MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
//     auto orthoErr = OM_->orthonormError(*Qmat_->data());
//     if (myrank==0)
//       std::cout << "orthoErr = " << orthoErr << std::endl;
// #endif

//     assert(computedRank_ == nVecs);
  }

  // void computeThinInPlace(matrix_t & A) {}

  template <typename vector_t>
  void doLinSolve(const vector_t & rhs, vector_t & y)const {
    qr::impl::solve<vector_t, trcp_mat>(rhs, this->localR_, y);
  }


  template < typename vector_in_t, typename vector_out_t>
  void applyQTranspose(const vector_in_t & vecIn, vector_out_t & vecOut) const
  {
    constexpr auto beta  = ::pressio::utils::constants<sc_t>::zero();
    constexpr auto alpha = ::pressio::utils::constants<sc_t>::one();
    ::pressio::ops::product(::pressio::transpose(), alpha, *this->Qmat_, vecIn, beta, vecOut);
  }

  template < typename vector_in_t, typename vector_out_t>
  void applyRTranspose(const vector_in_t & vecIn, vector_out_t & y) const
  {
    constexpr auto beta  = ::pressio::utils::constants<sc_t>::zero();
    constexpr auto alpha = ::pressio::utils::constants<sc_t>::one();
    ::pressio::ops::product(::pressio::transpose(), alpha, *this->localR_, vecIn, beta, y);
  }

  // if R_type != wrapper of Teuchos::SerialDenseMatrix
  template <typename T = R_t>
  ::pressio::mpl::enable_if_t<
    !containers::predicates::is_dense_matrix_wrapper_teuchos<T>::value and !std::is_void<T>::value,
    const T &
  >
  getCRefRFactor() const {
    this->Rmat_ = std::make_shared<T>(this->localR_->values());
    return *this->Rmat_;
  }

  // if R_type == wrapper of Teuchos::SerialDenseMatrix
  template <typename T = R_t>
  ::pressio::mpl::enable_if_t<
    containers::predicates::is_dense_matrix_wrapper_teuchos<T>::value and !std::is_void<T>::value,
    const T &
  >
  getCRefRFactor() const {
    this->Rmat_ = std::make_shared<T>(*this->localR_, Teuchos::View);
    return *this->Rmat_;
  }

  const Q_t & getCRefQFactor() const {
    return *this->Qmat_;
  }

private:
  void createLocalRIfNeeded(int newsize){
    if (localR_.is_null() or
    	(localR_->numRows()!=newsize and localR_->numCols()!=newsize)){
      localR_ = Teuchos::rcp(new serden_mat_t(newsize, newsize) );
    }
  }

  template <typename map_t>
  void createQIfNeeded(const map_t & map, int cols){
    if (!Qmat_ or !Qmat_->data()->Map().SameAs(map))
      Qmat_ = std::make_shared<Q_t>(map, cols);
  }

private:
  tsqr_adaptor_type tsqrAdaptor_;
  trcp_mat localR_			= {};
  int computedRank_			= {};
  mutable std::shared_ptr<Q_t> Qmat_	= nullptr;
  mutable std::shared_ptr<R_t> Rmat_	= nullptr;
};

}}} // end namespace pressio::qr::impl
#endif  // QR_IMPL_EPETRA_QR_EPETRA_MULTI_VECTOR_TSQR_IMPL_HPP_
