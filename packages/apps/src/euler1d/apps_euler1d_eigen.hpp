/*
//@HEADER
// ************************************************************************
//
// apps_euler1d_eigen.hpp
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

#ifndef PRESSIOAPPS_Euler1D_EIGEN_HPP_
#define PRESSIOAPPS_Euler1D_EIGEN_HPP_

#include "Eigen/Dense"
#include "Eigen/SparseCore"
#include <iostream>

namespace pressio{ namespace apps{

class Euler1dEigen{
  using eigVec = Eigen::VectorXd;

  using ui_t = unsigned int;

public:
  using scalar_type	= double;
  using state_type	= eigVec;
  using velocity_type	= eigVec;
  using dense_matrix_type = Eigen::MatrixXd;

  using eig_sp_mat = Eigen::SparseMatrix<scalar_type, Eigen::RowMajor, int>;
  using jacobian_type	= eig_sp_mat;

  typedef Eigen::Triplet<scalar_type> Tr;

public:
  explicit Euler1dEigen(eigVec params, ui_t Ncell=1001)
    : mu_(params), Ncell_(Ncell){
    this->setup();
  }

  Euler1dEigen() = delete;
  ~Euler1dEigen() = default;

public:
  state_type const & getInitialState() const {
    return U0_;
  };

  void velocity(const state_type & u,
  		const scalar_type /* t */,
    velocity_type & rhs) const{

	// Compute primitive variables from state vector
    this->computePrimitiveVariables(u);

    // Compute velocity
    this->computeFluxContributions(rhs);
    this->computeSourceContributions(rhs);

  }

  velocity_type velocity(const state_type & u,
  			 const scalar_type t) const{
    velocity_type RR(Ncell_);
    this->velocity(u, t, RR);
    return RR;
  }

  void
  computePrimitiveVariables(const state_type & U) const
  {
	for (ui_t i=0; i<Ncell_; ++i)
	{
	  r_(i) = U(Nvars_ * i + 0) / S_(i);
	  u_(i) = U(Nvars_ * i + 1) / U(Nvars_ * i + 0);
	  E_(i) = U(Nvars_ * i + 2) / U(Nvars_ * i + 0);
	  c_(i) = std::sqrt(gamma_*gm1_*(E_(i) - 0.5 * u_(i) * u_(i)));
	  p_(i) = gm1_*(U(Nvars_ * i + 2)  -.5*U(Nvars_ * i + 1)*U(Nvars_ * i + 1)/U(Nvars_ * i + 0))/S_(i);
	}
  }

  void
  computeFluxContributions(velocity_type & rhs) const
  {
	scalar_type Lm, Lp, Lmax;

	scalar_type L[3];
	scalar_type DSP[3];
	scalar_type Dm[3];
	scalar_type Dp[3];
	scalar_type Dpp[3];
	scalar_type fl[3], fr[3];

	// compute flux at each interior cell interface
	for (ui_t i=1; i<Ncell_+1; ++i)
	{
		// Maximum eigenvalue of flux Jacobian at cell interface
		Lm 	= c_(i)    + std::fabs(u_(i));
		Lp 	= c_(i+1)  + std::fabs(u_(i+1));
		Lmax 	= std::max(Lm,Lp);


        if (spatial_order_ == 1)
        {
        	// 1st order
        	for (ui_t j=0; j< Nvars_; j++)
				DSP[j] 	= .5*Lmax*.5*(S_(i)+S_(i+1))*
							(U_(Nvars_*i +j)/S_(i)  -U_(Nvars_*(i+1) +j)/S_(i+1));
        }
        else if (spatial_order_ == 2)
        {
        	// 2nd order
        	std::cout << "Still need to implement this... " << std::endl;
        	exit(1);
        }
        else
        {
        	std::cout << "Must specify 1 or 2 for spatial order. " << std::endl;
        	exit(1);
        }

        // Compute cell interface flux
        this->inviscidFlux(i,fl);
        this->inviscidFlux(i,fr);
        for (ui_t j=0; j< Nvars_; j++)
          F_(Nvars_ * i + j) = 0.5 * (fl[j] + fr[j]) + DSP[j];
	}

	// compute flux at boundaries

	// TODO add other BC options?

	// Left boundary
	this->inlet(true);

	// Right boundary
    this->outlet(false);

    // Add flux contributions to residual
    for (ui_t i=1; i<Ncell_; ++i)
      for (ui_t j=0; j< Nvars_; j++)
    	rhs(Nvars_ * i + j) = (F_[Nvars_ * (i+1)+ j] - F_[Nvars_ * i + j]) / dx_;

  }

  void
  computeSourceContributions(velocity_type & rhs) const
  {
	 // Add source contributions to residual
	 for (ui_t i=1; i<Ncell_; ++i)
	   rhs(Nvars_ * i + 1) += p_(i) * dSdx_(i);

  }

  // Compute invisid flux fi at the center of cell i
  void
  inviscidFlux(ui_t i, scalar_type * fi) const
  {
	  fi[0] = S_[i]*r_[i]*u_[i];
	  fi[1] = S_[i]*(r_[i]*u_[i]*u_[i]  +p_[i]);
	  fi[2] = S_[i]*(r_[i]*E_[i]  +p_[i])*u_[i];
  }

  // Compute flux fi for an inlet boundary condition given which side of the domain its on
  void
  inlet(bool is_left) const
  {
	  ui_t ind = 0;
	  if (!is_left)
	    ind = Ncell_-1;

		/* Compute inlet conditions */
	  scalar_type ptL = pt_in_; //3*101325;
	  scalar_type TtL = Tt_in_; //289740/1004.5;

	  scalar_type JiL = -u_(ind) + 5*std::pow(gamma_*p_[ind]/r_[ind],0.5);
	  scalar_type tmp1 = gamma_*R_air_*TtL - 0.5*gm1_*JiL*JiL;
	  scalar_type tmp2 = -4*gogm1_*R_air_*TtL;
	  scalar_type tmp3 = 4*gogm12_*R_air_*TtL - JiL*JiL;
	  scalar_type Mi1 = (-tmp2 + std::pow(tmp2*tmp2 - 4*tmp1*tmp3,0.5))/(2*tmp1);
	  scalar_type Mi2 = (-tmp2 - std::pow(tmp2*tmp2 - 4*tmp1*tmp3,0.5))/(2*tmp1);
	  Mi1 = std::max(Mi1,Mi2);
	//	printf("Mi: %f\n",Mi1);
	  scalar_type pL = ptL*std::pow(1 + 0.5 * gm1_ * Mi1*Mi1,-gogm1_);
	  scalar_type TL = TtL/(1 + 0.5*gm1_*Mi1*Mi1);
	  scalar_type rL = pL/(R_air_*TL);
	  scalar_type cL = std::pow(gamma_*pL/rL,0.5);
	  scalar_type uL = Mi1*cL;

	  F_(Nvars_ * ind + 0) = S_[ind]*rL*uL;
	  F_(Nvars_ * ind + 1)   = S_[ind]*rL*uL*uL + S_[ind]*pL;
	  F_(Nvars_ * ind + 2)   = S_[ind]*(gogm1_*pL + 0.5*rL*uL*uL)*uL;
  }

  // Compute flux fi for an outlet boundary condition given which side of the domain its on
  void
  outlet(bool is_left) const
  {
	  ui_t ind = 0;
	  if (!is_left)
	    ind = Ncell_-1;

	  scalar_type pR = p_ex_; //2*101325;

	  // interior entropy
	  scalar_type SiR = p_[ind]/ std::pow(r_[ind],gamma_);
	  scalar_type JiR = u_[ind] + 5.*std::pow(1.4*p_[ind]/r_[ind],0.5);
	  scalar_type rR = std::pow(pR/SiR,invg_);
	  scalar_type cR = std::sqrt(gamma_*pR/rR);
	  scalar_type uR = JiR - 5.*cR;

	  F_(Nvars_ * ind + 0) = S_[ind]*rR*uR;
	  F_(Nvars_ * ind + 1) = S_[ind]*rR*uR*uR + S_[ind]*pR;
	  F_(Nvars_ * ind + 2) = S_[ind]*(gogm1_*pR + 0.5*rR*uR*uR)*uR;

  }

  // TODO source function


  // TODO post-processing: outputing conserved variables, primitives?

//  // computes: A = Jac B where B is a Eigen::MatrixXd
//  void applyJacobian(const state_type & y,
//		     const dense_matrix_type & B,
//		     scalar_type t,
//		     dense_matrix_type & A) const{
//    auto JJ = jacobian(y, t);
//    // std::cout << "ApplyJacobian" << std::endl;
//    // std::cout << JJ << std::endl;
//    // multiply
//    A = JJ * B;
//  }
//
//  // computes: A = Jac B where B is a Eigen::MatrixXd
//  dense_matrix_type applyJacobian(const state_type & y,
//				  const dense_matrix_type & B,
//				  scalar_type t) const{
//    dense_matrix_type A( y.size(), B.cols() );
//    applyJacobian(y, B, t, A);
//    return A;
//  }
//
//  void jacobian(const state_type & u,
//		const scalar_type /*t*/,
//    jacobian_type & jac) const
//  {
//    //evaluate jacobian
//    if (jac.rows() == 0 || jac.cols()==0 ){
//      jac.resize(u.size(), u.size());
//    }
//    tripletList.clear();
//    tripletList.push_back( Tr( 0, 0, -dxInv_*u(0)) );
//    for (ui_t i=1; i<Ncell_; ++i){
//      tripletList.push_back( Tr( i, i-1, dxInv_ * u(i-1) ) );
//      tripletList.push_back( Tr( i, i, -dxInv_ * u(i) ) );
//    }
//    jac.setFromTriplets(tripletList.begin(), tripletList.end());
//  }
//
//  jacobian_type jacobian(const state_type & u,
//			 const scalar_type t) const{
//
//    jacobian_type JJ(u.size(), u.size());
//    this->jacobian(u, t, JJ);
//    return JJ;
//  }




private:
  void setup(){
    dx_ = (xR_ - xL_)/static_cast<scalar_type>(Ncell_);
    dxInv_ = 1.0/dx_;

    // TODO set physical constants

    // grid

    // TODO read from file, determine dSdx

    // Cell centers
    x_.resize(Ncell_);
    for (ui_t i=0; i<Ncell_; ++i)
      x_(i) = xL_ + dx_*i + 0.5 * dx_;

    // Cell cross-sectional areas
    S_.resize(Ncell_);
    dSdx_.resize(Ncell_);

    scalar_type dS = 1.0 / (xR_ - xL_);
    for (ui_t i=0; i<Ncell_; ++i)
    {
      S_(i) = 1.0 + (x_(i) - xL_) * dS;
      dSdx_(i) = dS;
    }

    // init condition

    // TODO read from restart file?

    scalar_type rL = pt_in_/(R_air_*Tt_in_);
    scalar_type uL = 0.0;
    scalar_type EL = R_air_*Tt_in_/(gamma_-1);

    scalar_type rR = rL;
    scalar_type uR = 0.0;
    scalar_type ER = EL;

    U_.resize(Ncell_ * Nvars_);
    for (ui_t i=0; i<Ncell_; ++i)
    {
      if (x_(i) < 0.0)
      {
    	  U_(Nvars_ * i + 0) = rL * S_[i];
    	  U_(Nvars_ * i + 1) = rL * uL * S_[i];
    	  U_(Nvars_ * i + 2) = rL * EL * S_[i];
      }
      else
      {
    	  U_(Nvars_ * i + 0) = rR * S_[i];
    	  U_(Nvars_ * i + 1) = rR * uR * S_[i];
    	  U_(Nvars_ * i + 2) = rR * ER * S_[i];
      }

    }

    U0_ = U_;

    // Initialize primitive variables
    r_.resize(Ncell_);
    u_.resize(Ncell_);
    p_.resize(Ncell_);

    // Initialize interior fluxes
    F_.resize(Ncell_+1);
    for (ui_t i=0; i<Ncell_+1; ++i)
      F_(i) = 0.0;

  };

private:

  /* Inputs */
  eigVec mu_; // parameters

  /* Solver Parameters */
  ui_t spatial_order_ = 2;  // 1 or 2..
  scalar_type CFL_ = 0.8;
  scalar_type Lt_ = 3.;
  scalar_type dt_ = 2.e-5;
  ui_t 	 Nt_ = Lt_/dt_;
  scalar_type kappa_;
  std::string limtype_ = "VanLeer";

  /* Physical Parameters */
  ui_t     Nvars_ = 3;
  scalar_type  gamma_ = 1.4;
  scalar_type  R_air_ = 287;


  scalar_type  gm1_ = gamma_ - 1.0;
  scalar_type  gogm1_ = gamma_ / gm1_;
  scalar_type  gogm12_ = gogm1_ / gm1_;
  scalar_type  invg_ = 1.0 / gamma_;



  /* State and fluxes */
  mutable state_type U_; // state vector
  mutable state_type U0_; // initial state vector
  mutable state_type r_; // density field
  mutable state_type u_; // velocity field
  mutable state_type E_; // energy field
  mutable state_type p_; // pressure field
  mutable state_type c_; // local speed of sound field
  mutable state_type F_; // cell interface fluxes

  /* Boundary Conditions */
  scalar_type pt_in_ = 3*101325.;
  scalar_type Tt_in_ = 289740/1004.5;
  scalar_type p_ex_ = 1.64*101325;


  /* Geometry variables */
  const scalar_type xL_ = 0.0; //left side of domain
  const scalar_type xR_ = 10.0; // right side of domain
  ui_t Ncell_; // # of cells
  scalar_type dx_; // cell size
  scalar_type dxInv_; // inv of cell size
  eigVec xCell_; // mesh points coordinates
  mutable state_type x_; // Cell centers
  mutable state_type S_; // Cell cross-sectional area
  mutable state_type dSdx_; // Cross-sectional area derivative
  bool 	constArea_ = true; // TODO turnoff when mesh reader is enabled...

};//end class

}} //namespace pressio::apps
#endif
