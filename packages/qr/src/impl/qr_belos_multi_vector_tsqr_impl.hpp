
#if defined HAVE_TRILINOS and defined HAVE_BELOS_TSQR
#ifndef QR_BELOS_MULTI_VECTOR_TSQR_IMPL_HPP_
#define QR_BELOS_MULTI_VECTOR_TSQR_IMPL_HPP_

#include "../qr_rfactor_solve_impl.hpp"
#include "BelosConfigDefs.hpp"
#include "BelosEpetraAdapter.hpp"
#include "BelosTpetraAdapter.hpp"
#include <BelosTsqrOrthoManager.hpp>

namespace rompp{ namespace qr{ namespace impl{

template<typename matrix_t, typename R_t, int n, int m,
	 typename MV_t, template<typename...> class Q_type>
class BelosMVTSQR<matrix_t, R_t, n, m, MV_t, Q_type, void>{

  using int_t	     = int;
  using sc_t	     = typename core::details::traits<matrix_t>::scalar_t;
  using ortho_t      = Belos::TsqrOrthoManager<sc_t, MV_t>;
  using serden_mat_t = Teuchos::SerialDenseMatrix<int_t, sc_t>;
  using trcp_mat     = Teuchos::RCP<serden_mat_t>;
  using Q_t	     = Q_type<MV_t>;
  const std::string label_ = "Belos";

public:
  BelosMVTSQR() = default;
  ~BelosMVTSQR() = default;

  void computeThinOutOfPlace(matrix_t & A) {
    auto nVecs = A.globalNumVectors();
    auto & ArowMap = A.getDataMap();
    createQIfNeeded(ArowMap, nVecs);
    createLocalRIfNeeded(nVecs);
    computedRank_ = OM_->normalizeOutOfPlace(*A.data(),
					     *Qmat_->data(),
					     localR_);
#ifdef DEBUG_PRINT
    int myrank{};
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    auto orthoErr = OM_->orthonormError(*Qmat_->data());
    if (myrank==0)
      std::cout << "orthoErr = " << orthoErr << std::endl;
#endif

    assert(computedRank_ == nVecs);
  }

  void computeThinInPlace(matrix_t & A) {
    auto nVecs = A.globalNumVectors();
    createLocalRIfNeeded(nVecs);
    computedRank_ = OM_->normalize(*A.data(), localR_);
    assert(computedRank_ == nVecs);
  }

  template <typename vector_t>
  void doLinSolve(const vector_t & rhs, vector_t & y)const {
      qr::impl::solve<vector_t, trcp_mat, n>(rhs, this->localR_, y);
  }


  template < typename vector_in_t, typename vector_out_t>
  void project(const vector_in_t & vecIn,
  		   vector_out_t & vecOut) const{
    core::ops::dot( *this->Qmat_, vecIn, vecOut );
  }

  // if R_type != wrapper of Teuchos::SerialDenseMatrix
  template <typename T = R_t,
  	    core::meta::enable_if_t<
  	      !core::meta::is_dense_matrix_wrapper_teuchos<T>::value and
	      !std::is_void<T>::value
  	      > * = nullptr>
  const T & getCRefRFactor() const {
    this->Rmat_ = std::make_shared<T>(this->localR_->values());
    return *this->Rmat_;
  }

  // if R_type == wrapper of Teuchos::SerialDenseMatrix
  template <typename T = R_t,
  	    core::meta::enable_if_t<
  	      core::meta::is_dense_matrix_wrapper_teuchos<T>::value and
	      !std::is_void<T>::value
  	      > * = nullptr>
  const T & getCRefRFactor() const {
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
    if (!Qmat_ or !Qmat_->hasRowMapEqualTo(map) )
      Qmat_ = std::make_shared<Q_t>(map, cols);
  }

private:
  std::shared_ptr< ortho_t > OM_	= std::make_shared<ortho_t>(label_);
  trcp_mat localR_			= {};
  int computedRank_			= {};

  // todo: these must be moved somewhere else
  mutable std::shared_ptr<Q_t> Qmat_	= nullptr;
  mutable std::shared_ptr<R_t> Rmat_	= nullptr;
};

}}} // end namespace rompp::qr::impl
#endif
#endif