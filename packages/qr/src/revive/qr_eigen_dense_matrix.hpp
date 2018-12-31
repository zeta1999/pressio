
#ifndef QR_EIGEN_DENSE_MATRIX_HPP_
#define QR_EIGEN_DENSE_MATRIX_HPP_

#include "qr_ConfigDefs.hpp"
#include "qr_forward_declarations.hpp"
#include "qr_solver_base.hpp"
// #include "../../CORE_ALL"
#include <Eigen/OrderingMethods>
#include<Eigen/SparseQR>

namespace rompp{ namespace qr{

// the input data is a wrapper for eigen DENSE matrix
template<typename matrix_type,
	 typename R_type,
	 template <typename...> class Q_type>
class QRSolver<matrix_type,
	       ::rompp::qr::Householder,
	       R_type,
	       Q_type,
	       core::meta::enable_if_t<
		 core::meta::is_eigen_dense_matrix_wrapper<matrix_type>::value and
		 core::meta::is_core_matrix_wrapper<R_type>::value and
		 core::details::traits<R_type>::is_shared_mem and
		 core::details::traits<R_type>::is_dense
		 >
	       >
  : public QRSolverBase<QRSolver<matrix_type, ::rompp::qr::Householder, R_type, Q_type>,
			R_type, Q_type<Eigen::MatrixXd>, matrix_type>{

  using sc_t = typename core::details::traits<matrix_type>::scalar_t;
  using Q_t = Q_type<Eigen::MatrixXd>;
  using this_t = QRSolver<matrix_type, ::rompp::qr::Householder, R_type, Q_type>;
  using base_t = QRSolverBase<this_t, R_type, Q_t, matrix_type>;
  friend base_t;
  using base_t::Qmat_;
  using base_t::Rmat_;

public:
  QRSolver() = default;
  ~QRSolver() = default;

private:

  void computeThinImpl(matrix_type & A){

    auto m = A.rows();
    auto n = A.cols();

    using native_mat_type = typename core::details::traits<matrix_type>::wrapped_t;
    Eigen::HouseholderQR<native_mat_type> eQR(*A.data());

    auto Qm = eQR.householderQ() * Eigen::MatrixXd::Identity(m,n);
    Qmat_ = std::make_shared<Q_t>(Qm);

    auto & Rm = eQR.matrixQR().template triangularView<Eigen::Upper>();
    Rmat_ = std::make_shared<R_type>( Rm );
  }

  const Q_t & cRefQFactorImpl() const {
    return *Qmat_;
  }

  const R_type & cRefRFactorImpl() const {
    return *Rmat_;
  }

};//end class

}} // end namespace rompp::qr
#endif