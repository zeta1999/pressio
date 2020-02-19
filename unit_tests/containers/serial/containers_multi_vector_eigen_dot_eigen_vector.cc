
#include <gtest/gtest.h>
#include "pressio_containers.hpp"


TEST(containers_multi_vector_serial_eigen_dynamic_class,
     dot_WithEigenVector_dynamic){

  using eigdmat_t = Eigen::MatrixXd;
  using myMV_t = pressio::containers::MultiVector<eigdmat_t>;

  //construct multivector
  myMV_t A(6,3);
  A(0,0) = 1.; A(0,1) = 2.; A(0,2) = 3.;
  A(1,0) = 3.; A(1,1) = 2.; A(1,2) = 1.;
  A(2,0) = 0.; A(2,1) = 0.; A(2,2) = 1.;
  A(3,0) = 0.; A(3,1) = 1.; A(3,2) = 0.;
  A(4,0) = 1.; A(4,1) = 0.; A(4,2) = 0.;
  A(5,0) = 0.; A(5,1) = 1.; A(5,2) = 1.;
  
  // eigen vector
  pressio::containers::Vector<Eigen::VectorXd> b(6);
  b(0) = 1.; b(1) = 1.; b(2) = 1.;
  b(3) = 1.; b(4) = 2.; b(5) = 1.;
    
  pressio::containers::Vector<Eigen::VectorXd> c(3);
  constexpr auto beta  = ::pressio::utils::constants::zero<scalar_t>();
  constexpr auto alpha = ::pressio::utils::constants::one<scalar_t>();
  ::pressio::containers::ops::product(::pressio::transpose(), alpha, A, b, beta, c);
  

  for (auto i=0; i<c1.extent(0); i++)
    EXPECT_DOUBLE_EQ( c(i), c1(i));
  
}



TEST(containers_multi_vector_serial_eigen_dynamic_class,
     dot_WithEigenVector_static){

  using eigdmat_t = Eigen::MatrixXd;
  using myMV_t = pressio::containers::MultiVector<eigdmat_t>;

  //construct multivector
  myMV_t A(6,3);
  A(0,0) = 1.; A(0,1) = 2.; A(0,2) = 3.;
  A(1,0) = 3.; A(1,1) = 2.; A(1,2) = 1.;
  A(2,0) = 0.; A(2,1) = 0.; A(2,2) = 1.;
  A(3,0) = 0.; A(3,1) = 1.; A(3,2) = 0.;
  A(4,0) = 1.; A(4,1) = 0.; A(4,2) = 0.;
  A(5,0) = 0.; A(5,1) = 1.; A(5,2) = 1.;
  
  // eigen vector
  pressio::containers::Vector<Eigen::VectorXd> b(6);
  b(0) = 1.; b(1) = 1.; b(2) = 1.;
  b(3) = 1.; b(4) = 2.; b(5) = 1.;
  
  using eig_v_st = Eigen::Matrix<double, 3, 1>;
  eig_v_st a;
  pressio::containers::Vector<eig_v_st> c(a);
  constexpr auto beta  = ::pressio::utils::constants::zero<scalar_t>();
  constexpr auto alpha = ::pressio::utils::constants::one<scalar_t>();
  ::pressio::containers::ops::product(::pressio::transpose(), alpha, A, b, beta, c);

  ASSERT_EQ(c.extent(0), 3);
  EXPECT_DOUBLE_EQ( c(0), 6.);
  EXPECT_DOUBLE_EQ( c(1), 6.);
  EXPECT_DOUBLE_EQ( c(2), 6.);  
}
