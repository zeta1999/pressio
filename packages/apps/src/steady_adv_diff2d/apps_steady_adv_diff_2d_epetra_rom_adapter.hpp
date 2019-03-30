
#ifndef ROMPP_APPS_STEADY_ADV_DIFF_2D_EPETRA_rom_adapter_HPP_
#define ROMPP_APPS_STEADY_ADV_DIFF_2D_EPETRA_rom_adapter_HPP_

#include "apps_steady_adv_diff_2d_epetra.hpp"

namespace rompp{ namespace apps{

class SteadyAdvDiff2dEpetraRomAdapter{
  using mv_t = Epetra_MultiVector;

public:
  /* these types exposed because need to be detected */
  using scalar_type	= double;
  using state_type	= Epetra_Vector;
  using residual_type	= state_type;

public:
  template <typename ... Args>
  SteadyAdvDiff2dEpetraRomAdapter(Args&& ... args)
    : appObj_{std::forward<Args>(args)...}
  {
    appObj_.setup();
  }

public:
  Epetra_Map const & getDataMap()const {
    return appObj_.getDataMap();
  };

  std::shared_ptr<state_type>
  getState() const {
    return appObj_.getState();
  }

  void residual(const state_type & yState,
		residual_type & rhs) const{
    appObj_.assembleMatrix();
    appObj_.fillRhs();
    auto A = appObj_.getMatrix();
    A->Multiply(false, yState, rhs);
    // now, rhs = A*u so we just subtract f to obtain residual
    auto f = appObj_.getForcing();
    rhs.Update(-1., (*f), 1.0);
  }

  residual_type residual(const state_type & yState) const{
    residual_type R( appObj_.getDataMap() );
    residual(yState, R);
    return R;
  };

  // computes: C = Jac B where B is a multivector
  void applyJacobian(const state_type & yState,
		     const mv_t & B,
		     mv_t & C) const
  {
    appObj_.assembleMatrix();
    auto A = appObj_.getMatrix();
    assert( A->NumGlobalCols() == B.GlobalLength() );
    assert( C.GlobalLength() == A->NumGlobalRows() );
    assert( C.NumVectors() == B.NumVectors() );
    A->Multiply(false, B, C);
  }

  // computes: A = Jac B where B is a multivector
  mv_t applyJacobian(const state_type & yState,
		     const mv_t & B) const{
    mv_t C( appObj_.getDataMap(), B.NumVectors() );
    applyJacobian(yState, B, C);
    return C;
  };

private:
  SteadyAdvDiff2dEpetra appObj_;

};

}} //namespace rompp::apps
#endif
