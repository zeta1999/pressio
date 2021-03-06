/*
//@HEADER
// ************************************************************************
//
// solvers_nonlinear_compose.hpp
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

#ifndef SOLVERS_NONLINEAR_IMPL_SOLVERS_NONLINEAR_COMPOSE_HPP_
#define SOLVERS_NONLINEAR_IMPL_SOLVERS_NONLINEAR_COMPOSE_HPP_

#include "./operators/solvers_hessian_gradient_operators.hpp"
#include "./operators/solvers_residual_jacobian_operators.hpp"
#include "./correction_mixins/solvers_hessian_gradient_corrector.hpp"
#include "./correction_mixins/solvers_qr_corrector.hpp"
#include "./correction_mixins/solvers_rj_corrector.hpp"
#include "solver.hpp"

namespace pressio{ namespace solvers{ namespace nonlinear{ namespace impl{

struct NewtonRaphson{};
struct GaussNewton{};
struct GaussNewtonQR{};
struct LevenbergMarquardt{};
using LM = LevenbergMarquardt;

#ifdef PRESSIO_ENABLE_TPL_PYBIND11
struct GaussNewtonPy{};
#endif


// ----------------------------------------------------------------------------
// *** COMPOSE CORRECTOR ***
// ----------------------------------------------------------------------------
template<typename tag, typename ... Args>
struct composeCorrector;


// *** GAUSS-NEWTON residual-jacobian API ***
template<
  typename system_t, typename state_t, typename h_t, typename g_t, typename lin_solver_t
  >
struct composeCorrector<
  GaussNewton,
  mpl::enable_if_t<
    pressio::solvers::concepts::system_residual_jacobian<system_t>::value or
    pressio::solvers::concepts::system_fused_residual_jacobian<system_t>::value>,
  system_t, state_t, h_t, g_t, lin_solver_t
  >
{
  using r_t = typename system_t::residual_type;
  using j_t = typename system_t::jacobian_type;
  static constexpr auto norm = pressio::Norm::L2;
  using operators_t = HessianGradientOperatorsRJApi<h_t, g_t, r_t, j_t>;
  using type	    = HessianGradientCorrector<operators_t, state_t, lin_solver_t, norm>;
};


// *** GAUSS-NEWTON residual-jacobian API with ud_ops ***
template<
  typename system_t, typename state_t, typename h_t, typename g_t, typename lin_solver_t,
  typename ud_ops_t
  >
struct composeCorrector<
  GaussNewton,
  mpl::enable_if_t<
    pressio::solvers::concepts::system_residual_jacobian<system_t>::value or
    pressio::solvers::concepts::system_fused_residual_jacobian<system_t>::value>,
  system_t, state_t, h_t, g_t, lin_solver_t, ud_ops_t
  >
{
  using r_t = typename system_t::residual_type;
  using j_t = typename system_t::jacobian_type;
  static constexpr auto norm = pressio::Norm::L2;
  using operators_t = HessianGradientOperatorsRJApi<h_t, g_t, r_t, j_t, ud_ops_t>;
  using type	    = HessianGradientCorrector<operators_t, state_t, lin_solver_t, norm>;
};

// *** GAUSS-NEWTON hessian-gradient API ***
template<
  typename system_t, typename state_t, typename h_t, typename g_t, typename lin_solver_t
  >
struct composeCorrector<
  GaussNewton,
  mpl::enable_if_t<
    pressio::solvers::concepts::system_hessian_gradient<system_t>::value or
    pressio::solvers::concepts::system_fused_hessian_gradient<system_t>::value>,
  system_t, state_t, h_t, g_t, lin_solver_t
  >
{
  static constexpr auto norm = pressio::Norm::L2;
  using operators_t = HessianGradientOperatorsHGApi<h_t, g_t>;
  using type	    = HessianGradientCorrector<operators_t, state_t, lin_solver_t, norm>;
};

// *** GAUSS-NEWTON with QR solver ***
template<typename system_t, typename state_t, typename qr_solver_t>
struct composeCorrector<
  GaussNewtonQR,
  mpl::enable_if_t<
    pressio::solvers::concepts::system_residual_jacobian<system_t>::value or
    pressio::solvers::concepts::system_fused_residual_jacobian<system_t>::value>,
  system_t, state_t, qr_solver_t
  >
{
  using r_t = typename system_t::residual_type;
  using j_t = typename system_t::jacobian_type;
  static constexpr auto norm = pressio::Norm::L2;
  using operators_t = ResidualJacobianOperators<r_t, j_t>;
  using type	    = QRCorrector<operators_t, state_t, qr_solver_t, norm>;
};

// *** LEVENBERG-MARQUARDT Residual-jacobian API ***
template<
  typename system_t, typename state_t, typename h_t, typename g_t, typename lin_solver_t
  >
struct composeCorrector<
  LM,
  mpl::enable_if_t<
    pressio::solvers::concepts::system_residual_jacobian<system_t>::value or
    pressio::solvers::concepts::system_fused_residual_jacobian<system_t>::value>,
  system_t, state_t, h_t, g_t, lin_solver_t
  >
{
  using r_t = typename system_t::residual_type;
  using j_t = typename system_t::jacobian_type;
  static constexpr auto norm = pressio::Norm::L2;
  using operators_t = LMHessianGradientOperatorsRJApi<h_t, g_t, r_t, j_t>;
  using type	    = HessianGradientCorrector<operators_t, state_t, lin_solver_t, norm>;
};

// *** LEVENBERG-MARQUARDT hessian-gradient API ***
template<
  typename system_t, typename state_t, typename h_t, typename g_t, typename lin_solver_t
  >
struct composeCorrector<
  LM,
  mpl::enable_if_t<
    pressio::solvers::concepts::system_hessian_gradient<system_t>::value or
    pressio::solvers::concepts::system_fused_hessian_gradient<system_t>::value>,
  system_t, state_t, h_t, g_t, lin_solver_t
  >
{
  static constexpr auto norm = pressio::Norm::L2;
  using operators_t = LMHessianGradientOperatorsHGApi<h_t, g_t>;
  using type	    = HessianGradientCorrector<operators_t, state_t, lin_solver_t, norm>;
};

// *** Newton-Raphson API ***
template<typename system_t, typename state_t, typename lin_solver_t>
struct composeCorrector<
  NewtonRaphson,
  mpl::enable_if_t<
    pressio::solvers::concepts::system_residual_jacobian<system_t>::value or
    pressio::solvers::concepts::system_fused_residual_jacobian<system_t>::value>,
  system_t, state_t, lin_solver_t
  >
{
  using r_t = typename system_t::residual_type;
  using j_t = typename system_t::jacobian_type;
  static constexpr auto norm = pressio::Norm::L2;
  using operators_t = ResidualJacobianOperators<r_t, j_t>;
  using type = RJCorrector<operators_t, state_t, lin_solver_t, norm>;
};



// ----------------------------------------------------------------------------
// *** COMPOSE SOLVER TYPE ***
// ----------------------------------------------------------------------------
template<
  typename system_t,
  typename tag,
  template< typename...> class update,
  typename ... Args>
struct compose
{
  using type = void;
};

template<
  typename system_t,
  typename tag,
  template<typename...> class update_t,
  typename linear_solver_t
  >
struct compose<
  system_t, tag, update_t,
  mpl::enable_if_t<
    std::is_same<tag, GaussNewton>::value or std::is_same<tag, LM>::value
    >, linear_solver_t
  >
{
  // // GN or LM with neq need a valid linear solver for hess/grad system
  // using ic2 = ::pressio::mpl::variadic::find_if_unary_pred_t<
  //   ::pressio::solvers::concepts::linear_solver_for_least_squares_solver, Args...>;
  // using linear_solver_t = ::pressio::mpl::variadic::at_or_t<void, ic2::value, Args...>;
  // static_assert(!std::is_void<linear_solver_t>::value and ic2::value < sizeof... (Args),
  // 		"A valid linear solver type must be passed to GN with normal equations");
  static_assert(::pressio::solvers::concepts::linear_solver_for_least_squares_solver<linear_solver_t>::value,
  		"A valid linear solver type must be passed to GN with normal equations");

  using scalar_t = typename system_t::scalar_type;
  using state_t = typename system_t::state_type;
  // gradient is same as state_t
  using grad_t = state_t;
  // hessian_t is extracted from linear solver
  using hess_t = typename linear_solver_t::matrix_type;

  using corr_mixin = typename composeCorrector<
    tag, void, system_t, state_t, hess_t, grad_t, linear_solver_t>::type;

  // TODO: assert that the update is admissible for the tag
  using update_mixin  = update_t<scalar_t, state_t, corr_mixin>;
  using type = Solver<update_mixin, scalar_t>;
};


#ifdef PRESSIO_ENABLE_TPL_PYBIND11
template<
  typename system_t,
  template<typename...> class update_t,
  typename hessian_t
  >
struct compose<
  system_t, GaussNewton, update_t,
  mpl::enable_if_t<
    // when dealing with pressio4py and the GN solver is created for solving steady or unsteady LSPG,
    // the system_t is NOT a pybind object.
    // The system class is a python object only if one is trying to use this GaussNewton
    // to solve a nonlinear system that is writte in python. For that case, I would say they are
    // better off using other Python libraries for solvers, so disable that scenario for now.
    !::pressio::ops::predicates::is_object_pybind<system_t>::value and
    ::pressio::containers::predicates::is_matrix_wrapper_pybind<hessian_t>::value>,
  pybind11::object, hessian_t
  >
{
  using scalar_t = typename system_t::scalar_type;
  using state_t = typename system_t::state_type;
  using grad_t = state_t;
  using linear_solver_t = pybind11::object;

  using corr_mixin = typename composeCorrector<
    GaussNewton, void, system_t, state_t, hessian_t, grad_t, linear_solver_t>::type;
  using update_mixin  = update_t<scalar_t, state_t, corr_mixin>;
  using type = Solver<update_mixin, scalar_t>;
};
#endif


template<
  typename system_t,
  typename tag,
  template<typename...> class update_t,
  typename linear_solver_t,
  typename ud_ops_t
  >
struct compose<
  system_t, tag, update_t,
  mpl::enable_if_t<
    (std::is_same<tag, GaussNewton>::value or std::is_same<tag, LM>::value) and
    !::pressio::containers::predicates::is_wrapper<ud_ops_t>::value
    >, linear_solver_t, ud_ops_t
  >
{
  static_assert(::pressio::solvers::concepts::linear_solver_for_least_squares_solver<linear_solver_t>::value,
  		"A valid linear solver type must be passed to GN with normal equations");

  // todo: check that ops type is admissible

  using scalar_t = typename system_t::scalar_type;
  using state_t = typename system_t::state_type;
  // gradient is same as state_t
  using grad_t = state_t;
  // hessian_t is extracted from linear solver
  using hess_t = typename linear_solver_t::matrix_type;

  using corr_mixin = typename composeCorrector<
    tag, void, system_t, state_t, hess_t, grad_t, linear_solver_t, ud_ops_t>::type;
  // TODO: assert that the update is admissible for the tag
  using update_mixin  = update_t<scalar_t, state_t, corr_mixin>;
  using type = Solver<update_mixin, scalar_t>;
};

template<
  typename system_t,
  template<typename...> class update_t,
  typename ... Args
  >
struct compose<system_t, GaussNewtonQR, update_t, Args...>
{
  // verify the sequence contains a valid QR solver type
  using ic2 = ::pressio::mpl::variadic::find_if_unary_pred_t<
    ::pressio::solvers::concepts::qr_solver_for_gn_qr, Args...>;
  using qr_solver_t = ::pressio::mpl::variadic::at_or_t<void, ic2::value, Args...>;
  static_assert(!std::is_void<qr_solver_t>::value and
		ic2::value < sizeof... (Args),
  		"A valid QR solver type must be passed to compose a QR-based GN solver");
  using qr_solver_matrix_t = typename ::pressio::qr::details::traits<qr_solver_t>::matrix_t;

  using scalar_t = typename system_t::scalar_type;
  using state_t = typename system_t::state_type;

  using corr_mixin = typename composeCorrector<GaussNewtonQR, void, system_t, state_t, qr_solver_t>::type;
  // TODO: assert that the update is admissible for the tag
  using update_mixin  = update_t<scalar_t, state_t, corr_mixin>;
  using type = Solver<update_mixin, scalar_t>;
};

template<
  typename system_t,
  template<typename...> class update_t,
  typename linear_solver_t,
  typename ... Args
  >
struct compose<system_t, NewtonRaphson, update_t, linear_solver_t, Args...>
{
  using scalar_t = typename system_t::scalar_type;
  using state_t  = typename system_t::state_type;
  using corr_mixin = typename composeCorrector<NewtonRaphson, void, system_t, state_t, linear_solver_t>::type;
  // TODO: assert that the update is admissible for the tag
  using update_mixin  = update_t<scalar_t, state_t, corr_mixin>;
  using type = Solver<update_mixin, scalar_t>;
};

}}}}
#endif  // SOLVERS_NONLINEAR_IMPL_SOLVERS_NONLINEAR_COMPOSE_HPP_
