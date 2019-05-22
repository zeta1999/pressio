
#ifndef ROMPP_APPS_NONLIN_ADV_DIFF_REACTION_FLAME_2D_EIGEN_HPP_
#define ROMPP_APPS_NONLIN_ADV_DIFF_REACTION_FLAME_2D_EIGEN_HPP_

#include "../../../CORE_ALL"
#include "Eigen/Dense"
#include "Eigen/SparseCore"
#include <cmath>

namespace rompp{ namespace apps{

class UnsteadyNonLinAdvDiffReacFlame2dEigen{
protected:
  using this_t		= UnsteadyNonLinAdvDiffReacFlame2dEigen;
  using nativeVec	= Eigen::VectorXd;
  using mv_t		= Eigen::MatrixXd;


public:
  using scalar_type	= double;
  using state_type	= nativeVec;
  using residual_type	= state_type;
  using jacobian_type	= Eigen::SparseMatrix<scalar_type, Eigen::RowMajor, int>;

  typedef Eigen::Triplet<scalar_type> Tr;
  using mat4_t		= Eigen::Matrix<scalar_type, 4, 4>;
  static constexpr auto zero = ::rompp::core::constants::zero<scalar_type>();
  static constexpr auto one = ::rompp::core::constants::one<scalar_type>();
  static constexpr auto two = ::rompp::core::constants::two<scalar_type>();

public:
  UnsteadyNonLinAdvDiffReacFlame2dEigen
  (int Nx, int Ny,
   scalar_type K	= static_cast<scalar_type>(2), // cm^2/s
   scalar_type preExp	= static_cast<scalar_type>(5.5*1e12),
   scalar_type E	= static_cast<scalar_type>(8.0*1000.),
   std::array<scalar_type,3> W = {{2.016, 31.9, 18}})
    : K_{K},
      preExp_{preExp},
      negE_{-E},
      W_(W),
      rhoOvWH2{rho_/W_[0]},
      rhoOvWO2{rho_/W_[1]},
      WH2ovRho{W_[0]/rho_},
      WO2ovRho{W_[1]/rho_},
      WH2OovRho{W_[2]/rho_},
      Nx_{Nx},
      Ny_{Ny},
      dx_{Lx_/(Nx)},
      dy_{Ly_/(Ny)},
      dxSqInv_{one/(dx_*dx_)},
      dySqInv_{one/(dy_*dy_)},
      dx2Inv_{one/(two*dx_)},
      dy2Inv_{one/(two*dy_)}
  {}

public:
  state_type const & getInitialState() const{
    return state_;
  };

  void setupPhysicalGrid();
  void setupFields();
  void setup();

  int getUnknownCount() const{ return this_t::numSpecies_*Nx_*Ny_; }
  scalar_type getDx() const { return dx_; };
  scalar_type getDy() const { return dy_; };
  nativeVec getX() const { return x_; }
  nativeVec getY() const { return y_; }
  nativeVec getU() const { return u_; }
  nativeVec getV() const { return v_; }

public:
  void residual(const state_type & yState,
		residual_type & rhs,
		scalar_type t) const{
    residual_impl(yState, rhs);
  }

  residual_type residual(const state_type & yState,
			 scalar_type t) const{
    residual_type R(numDof_);
    residual_impl(yState, R);
    return R;
  };

  void jacobian(const state_type & u,
  		jacobian_type & jac,
  		const scalar_type /*t*/) const{
    this->jacobian_impl(u, jac);
  }

  jacobian_type jacobian(const state_type & u,
  			 const scalar_type t) const{
    jacobian_type JJ(u.size(), u.size());
    this->jacobian_impl(u, JJ);
    return JJ;
  }

  // // computes: C = Jac B where B is a multivector
  // void applyJacobian(const state_type & yState,
  // 		     const mv_t & B,
  // 		     mv_t & C,
  // 		     scalar_type t) const{
  //   applyJacobian_impl(yState, B, C);
  // }

  // mv_t applyJacobian(const state_type & yState,
  // 		     const mv_t & B,
  // 		     scalar_type t) const{
  //   mv_t A( yState.size(), B.cols() );
  //   applyJacobian_impl(yState, B, A);
  //   return A;
  // };

private:
  void globalIDToGiGj(int ID, int & gi, int & gj) const{
    gj = ID/Nx_;
    gi = ID % Nx_;
  }

  void compute_dsdw(scalar_type wT, scalar_type wH2,
		    scalar_type wO2, scalar_type wH2O)const;
  void compute_sources(scalar_type wT, scalar_type wH2,
		       scalar_type wO2, scalar_type wH2O)const;

  void residual_impl(const state_type & yState,
		     residual_type & R) const;

  void jacobian_impl(const state_type & yState,
  		     jacobian_type & J) const;

  // void applyJacobian_impl(const state_type & yState,
  // 			  const mv_t & B,
  // 			  mv_t & C) const{
  //   jacobian_type JJ(yState.size(), yState.size());
  //   JJ.setZero();
  //   this->jacobian_impl(yState, JJ);
  //   C = JJ * B;
  // }

protected:
  // T, H2, O2, H2O
  static constexpr int numSpecies_{4};
  const scalar_type Lx_{1.8};
  const scalar_type Ly_{0.9};
  const std::array<scalar_type,2> xAxis_{{0., Lx_}};
  const std::array<scalar_type,2> yAxis_{{0., Ly_}};

  const scalar_type initTemp_ = {300};
  const scalar_type rho_{0.00139};
  const scalar_type gasR_{8.314};
  const scalar_type Q_{9800};
  const scalar_type K_{};
  const scalar_type preExp_ = {};
  const scalar_type negE_ = {};
  const std::array<scalar_type,3> W_ = {};
  const scalar_type rhoOvWH2 = {};
  const scalar_type rhoOvWO2 = {};
  const scalar_type WH2ovRho = {};
  const scalar_type WO2ovRho = {};
  const scalar_type WH2OovRho = {};

  const std::array<scalar_type,numSpecies_> bcLeftGamma13_{{300, 0, 0, 0}};
  const std::array<scalar_type,numSpecies_> bcLeftGamma2_{{950, 0.0282, 0.2259, 0}};

  // grid points (cell centered)
  const int Nx_{};
  const int Ny_{};

  const scalar_type dx_{};
  const scalar_type dy_{};
  const scalar_type dxSqInv_{};
  const scalar_type dySqInv_{};
  const scalar_type dx2Inv_{};
  const scalar_type dy2Inv_{};

  // note that dof refers to the degress of freedom,
  // which is NOT same as grid points. for this problem,
  // the dof = numSpecies_ * number_of_unknown_grid_points
  int numGpt_;
  int numDof_;

  nativeVec x_;
  nativeVec y_;
  nativeVec u_;
  nativeVec v_;
  mutable std::vector<Tr> tripletList;
  mutable jacobian_type J_;
  mutable nativeVec state_;

  mutable std::array<scalar_type, numSpecies_> s_;
  mutable mat4_t dsdw_;
  // to label a grid point based on whether it belongs
  //to lower, mid or upper region of the flow, see diagram in paper
  mutable nativeVec regionLabel_;
};

}} //namespace rompp::apps
#endif
