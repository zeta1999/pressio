/*
//@HEADER
// ************************************************************************
//
// ode_discrete_time_jacobian_impl.hpp
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

#ifndef ODE_IMPLICIT_IMPL_ODE_DISCRETE_TIME_JACOBIAN_IMPL_HPP_
#define ODE_IMPLICIT_IMPL_ODE_DISCRETE_TIME_JACOBIAN_IMPL_HPP_

namespace pressio{ namespace ode{ namespace impl{

template <typename jacobian_type, typename scalar_type>
::pressio::mpl::enable_if_t<
  // (std::is_same<stepper_tag, ::pressio::ode::implicitmethods::Euler>::value) and
  (containers::predicates::is_sparse_matrix_wrapper_eigen<jacobian_type>::value or
#ifdef PRESSIO_ENABLE_TPL_TRILINOS
   containers::predicates::is_sparse_matrix_wrapper_epetra<jacobian_type>::value or
#endif
  containers::predicates::is_dense_matrix_wrapper_eigen<jacobian_type>::value)
>
discrete_time_jacobian(jacobian_type & jac, const scalar_type & dt, ::pressio::ode::implicitmethods::Euler)
{
  constexpr auto cn   = ::pressio::ode::constants::bdf1<scalar_type>::c_n_;
  const auto cf	  = ::pressio::ode::constants::bdf1<scalar_type>::c_f_ * dt;
  // jac.scale(cf);
  // jac.addToDiagonal(cn);
  ::pressio::ops::scale(jac, cf);
  ::pressio::ops::addToDiagonal(jac, cn);
}


#ifdef PRESSIO_ENABLE_TPL_PYBIND11
template <typename jacobian_type, typename scalar_type>
::pressio::mpl::enable_if_t<
  // (std::is_same<stepper_tag, ::pressio::ode::implicitmethods::Euler>::value) and
  containers::predicates::is_array_pybind11<jacobian_type>::value
>
discrete_time_jacobian(jacobian_type & jac, const scalar_type & dt, ::pressio::ode::implicitmethods::Euler)
{
  using namespace ::pressio::ode::constants;

  constexpr auto cn   = ::pressio::ode::constants::bdf1<scalar_type>::c_n_;
  const auto cf	  = ::pressio::ode::constants::bdf1<scalar_type>::c_f_ * dt;

  if (jac.ndim() != 2)
    throw std::runtime_error("Tensors with dim>2 not supported");

  double *ptr = jac.mutable_data();
  const size_t rows = jac.shape()[0];
  const size_t cols = jac.shape()[1];

  for (size_t irow = 0; irow < rows; irow++){
    for (size_t icol = 0; icol < cols; icol++)
    {
      ptr[irow*cols + icol] *= cf;
      if (irow == icol and irow == 2)
	ptr[irow*cols + icol] += cn;
    }
  }
}
#endif


template <typename jacobian_type, typename scalar_type>
::pressio::mpl::enable_if_t<
  // (std::is_same<stepper_tag, ::pressio::ode::implicitmethods::BDF2>::value) and
  containers::predicates::is_sparse_matrix_wrapper_eigen<jacobian_type>::value
>
discrete_time_jacobian(jacobian_type & jac, const scalar_type & dt, ::pressio::ode::implicitmethods::BDF2)
{
  constexpr auto cn   = ::pressio::ode::constants::bdf2<scalar_type>::c_n_;
  const auto cf	  = ::pressio::ode::constants::bdf2<scalar_type>::c_f_ * dt;

  using namespace ::pressio::ode::constants;
  // jac.scale(cf);
  // jac.addToDiagonal(cn);
  ::pressio::ops::scale(jac, cf);
  ::pressio::ops::addToDiagonal(jac, cn);
}


}}}//end namespace pressio::ode::impl
#endif  // ODE_IMPLICIT_IMPL_ODE_DISCRETE_TIME_JACOBIAN_IMPL_HPP_
