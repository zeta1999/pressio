
#ifndef QR_UTEST_FIXTURES_HPP_
#define QR_UTEST_FIXTURES_HPP_

#include <gtest/gtest.h>
#include "CORE_ALL"
#include "Epetra_MpiComm.h"
#include "Eigen/Dense"
#include "qr_gold.hpp"

struct epetraR9Fixture
  : public ::testing::Test{

public:
  using nat_mv_t = Epetra_MultiVector;
  using mymvec_t = rompp::core::MultiVector<nat_mv_t>;
  using nat_v_t = Epetra_Vector;
  using myvec_t = rompp::core::Vector<nat_v_t>;

  std::shared_ptr<Epetra_MpiComm> comm_;
  std::shared_ptr<Epetra_Map> rowMap_;
  std::shared_ptr<mymvec_t> A_;
  std::shared_ptr<myvec_t> v_;

  int rank_;
  int numProc_;
  int numGlobalEntries_;
  const int numVectors_ = 4;

  virtual void SetUp(){
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    comm_ = std::make_shared<Epetra_MpiComm>(MPI_COMM_WORLD);
    rank_ = comm_->MyPID();
    numProc_ = comm_->NumProc();
    assert(numProc_ == 2);
    numGlobalEntries_ = 9;

    rowMap_ = std::make_shared<Epetra_Map>(numGlobalEntries_, 0, *comm_);
    A_ = std::make_shared<mymvec_t>(*rowMap_, numVectors_);
    v_ = std::make_shared<myvec_t>(*rowMap_);

  }//setup

  void fillMatrix(){
    if(rank_==0){
      (*A_)(0,0) = 3.2; (*A_)(0,1) = 1.2;  (*A_)(0,2) = 1.;
      (*A_)(1,0) = 1.2; (*A_)(1,2) = -2.2;
      (*A_)(2,1) = 4.0; (*A_)(2,3) = -2.;
      (*A_)(3,1) = 4.;
      (*A_)(4,2) = -1.; (*A_)(4,3) = -4.;
    }
    if(rank_==1){
      (*A_)(0,0) = 0.2; (*A_)(0,1) = 5.;     (*A_)(0,2) = 1.;
      (*A_)(1,0) = 1.;  (*A_)(1,1) = 1.1;    (*A_)(1,2) = 1.25; (*A_)(1,3) = -3.;
      (*A_)(2,2) = 1.;  (*A_)(2,1) = 0.1111; (*A_)(2,3) = 6.;
    }
    //A_->data()->Print(std::cout);
  }

  void fillVector(){
    v_->putScalar(1.0);
  }

  virtual void TearDown(){}
};
// --------------------------------------------


struct tpetraR9Fixture
  : public ::testing::Test{

public:
  using tcomm = Teuchos::Comm<int>;
  using map_t = Tpetra::Map<>;
  using nat_mvec_t = Tpetra::MultiVector<>;
  using ST = typename nat_mvec_t::scalar_type;
  using LO = typename nat_mvec_t::local_ordinal_type;
  using GO = typename nat_mvec_t::global_ordinal_type;
  using nat_vec_t = Tpetra::Vector<>;

  using mymvec_t = rompp::core::MultiVector<nat_mvec_t>;
  using mv_device_t = typename rompp::core::details::traits<mymvec_t>::device_t;
  using myvec_t = rompp::core::Vector<nat_vec_t>;

  // gold solution
  rompp::qr::test::qrGoldSol<double> gold;

  int rank_;
  int numProc_;
  const int numVectors_ = 4;
  int numGlobalEntries_;
  Teuchos::RCP<const tcomm> comm_;
  Teuchos::RCP<const map_t> contigMap_;
  std::shared_ptr<mymvec_t> A_;
  std::shared_ptr<myvec_t> v_;

  virtual void SetUp(){
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    comm_ = Teuchos::rcp (new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
    rank_ = comm_->getRank();
    numProc_ = comm_->getSize();

    assert(numProc_==2);

    numGlobalEntries_ = 9;
    contigMap_ = Teuchos::rcp(new map_t(numGlobalEntries_,0,comm_));
    A_ = std::make_shared<mymvec_t>(contigMap_, numVectors_);
    v_ = std::make_shared<myvec_t>(contigMap_);
  }

  void fillMatrix(){
    // get trilinos tpetra multivector object
    auto trilD = A_->data();
    trilD->sync<Kokkos::HostSpace>();

    /*--------------------------------------------
     * (1): modify the host view and then sync
     * most likely, host and device will be same unless we run CUDA
     * so in theory we should not worry about syncing but it
     * does not hurt to do it anyway
     //--------------------------------------------*/
    auto v2d = trilD->getLocalView<Kokkos::HostSpace>();
    auto c0 = Kokkos::subview(v2d, Kokkos::ALL(), 0);
    //we are going to change the host view
    trilD->modify<Kokkos::HostSpace>();

    if(rank_==0){
      v2d(0,0) = 3.2; v2d(0,1) = 1.2;  v2d(0,2) = 1.;
      v2d(1,0) = 1.2; v2d(1,2) = -2.2;
      v2d(2,1) = 4.0; v2d(2,3) = -2.;
      v2d(3,1) = 4.;
      v2d(4,2) = -1.; v2d(4,3) = -4.;
    }
    if(rank_==1){
      v2d(0,0) = 0.2; v2d(0,1) = 5.;     v2d(0,2) = 1.;
      v2d(1,0) = 1.;  v2d(1,1) = 1.1;    v2d(1,2) = 1.25; v2d(1,3) = -3.;
      v2d(2,2) = 1.;  v2d(2,1) = 0.1111; v2d(2,3) = 6.;
    }
    // sync from host to device
    trilD->sync<mv_device_t>();
  }

  void fillVector(){
    v_->putScalar(1.0);
  }

  virtual void TearDown(){}
};


#endif /* QR_FIXTURES_HPP_ */
