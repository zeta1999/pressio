
#include "tpetra_only_fixtures.hpp"

TEST_F(tpetraMultiVectorR9C4VecS9Fixture, MVVecDotProduct){
  using namespace pressio;

  // --------------------------------------------
  // construct and fill multivector wrapper
  // --------------------------------------------
  using nat_mv_t = typename tpetraMultiVectorR9C4VecS9Fixture::mvec_t;
  using mvec_t = containers::MultiVector<nat_mv_t>;
  using mv_device_t = typename containers::details::traits<mvec_t>::device_t;

  mvec_t MV( *mv_ );
  MV.setZero();
  // get trilinos tpetra multivector object
  auto trilD = MV.data();
  trilD->sync<Kokkos::HostSpace>();

  /*--------------------------------------------
   * (1): modify the host view and then sync
  //--------------------------------------------*/

  auto v2d = trilD->getLocalView<Kokkos::HostSpace>();
  //we are going to change the host view
  trilD->modify<Kokkos::HostSpace>();

  if(rank_==0){
    v2d(0,0) = 3.2; v2d(0,1) = 1.2;
    v2d(1,0) = 1.2;
    v2d(2,1) = 4;
  }
  if(rank_==1){
    v2d(2,2) = 3;
  }
  // sync from host to device
  trilD->sync<mv_device_t> ();

  // --------------------------------------------
  // construct and fill vector wrapper
  // --------------------------------------------
  using nat_v_t = typename tpetraMultiVectorR9C4VecS9Fixture::vec_t;
  using vec_t = containers::Vector<nat_v_t>;
  using v_device_t = typename containers::details::traits<vec_t>::device_t;

  //check that device_t of mv and vec match
  ::testing::StaticAssertTypeEq<mv_device_t, v_device_t>();

  vec_t VV( *x_ );
  VV.setZero();
  // get trilinos tpetra vector object
  auto vvtrilD = VV.data();
  vvtrilD->putScalar(static_cast<double>(1));
  // sync from host to device
  vvtrilD->sync<v_device_t>();

  //--------------------------------------------

  //MV dot b = c
  //  auto c = containers::ops::dot(MV, b);
  // EXPECT_EQ((int) c.size(), 4);
  // EXPECT_NEAR(c[0], 4.4, 1e-12);
  // EXPECT_NEAR(c[1], 5.2, 1e-12);
  // EXPECT_NEAR(c[2], 3., 1e-12);
  // EXPECT_NEAR(c[3], 0., 1e-12);

  // //MV dot b = c2
  containers::Vector<Eigen::VectorXd> c2;
  c2.resize(4);
  containers::ops::dot(MV, VV, c2);
  EXPECT_EQ( c2.size(), 4);
  EXPECT_NEAR(c2[0], 4.4, 1e-12);
  EXPECT_NEAR(c2[1], 5.2, 1e-12);
  EXPECT_NEAR(c2[2], 3., 1e-12);
  EXPECT_NEAR(c2[3], 0., 1e-12);

  // store into teuchos serial dense vector wrapper
  //MV dot b = c3
  using natvec_t3 = Teuchos::SerialDenseVector<int, double>;
  containers::Vector<natvec_t3> c3(4);
  containers::ops::dot(MV, VV, c3);
  EXPECT_EQ( c3.size(), 4);
  EXPECT_NEAR(c3[0], 4.4, 1e-12);
  EXPECT_NEAR(c3[1], 5.2, 1e-12);
  EXPECT_NEAR(c3[2], 3., 1e-12);
  EXPECT_NEAR(c3[3], 0., 1e-12);

}