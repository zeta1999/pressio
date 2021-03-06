/*
//@HEADER
// ************************************************************************
//
// rom_lspg_steady_common_traits.hpp
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

#ifndef ROM_LSPG_IMPL_STEADY_TRAITS_ROM_LSPG_STEADY_COMMON_TRAITS_HPP_
#define ROM_LSPG_IMPL_STEADY_TRAITS_ROM_LSPG_STEADY_COMMON_TRAITS_HPP_

namespace pressio{ namespace rom{ namespace lspg{ namespace impl{ namespace steady{

template <typename fom_system_t, typename lspg_state_t, typename enable = void>
struct ExtractNativeHelper;

template <typename fom_system_t, typename lspg_state_t>
struct ExtractNativeHelper<
  fom_system_t, lspg_state_t
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
  ,mpl::enable_if_t<
     !::pressio::containers::predicates::is_vector_wrapper_pybind<lspg_state_t>::value
    >
#endif
  >
{
  using fom_native_state_t    = typename fom_system_t::state_type;
  using fom_native_residual_t = typename fom_system_t::residual_type;
};

#ifdef PRESSIO_ENABLE_TPL_PYBIND11
template <typename fom_system_t, typename lspg_state_t>
struct ExtractNativeHelper<
  fom_system_t, lspg_state_t,
  mpl::enable_if_t<
    ::pressio::containers::predicates::is_vector_wrapper_pybind<lspg_state_t>::value
    >
  >
{
  using fom_native_state_t     = typename ::pressio::containers::details::traits<lspg_state_t>::wrapped_t;
  using fom_native_residual_t = fom_native_state_t;
};
#endif
// ------------------------------------------------------------------------------------


template <
  typename fom_system_type,
  typename decoder_type,
  typename lspg_state_type,
  typename ...Args
  >
struct CommonTraits
{
  // the scalar type
  using scalar_t = typename ::pressio::containers::details::traits<lspg_state_type>::scalar_t;

  using fom_system_t		= fom_system_type;
  using fom_native_state_t	= typename ExtractNativeHelper<fom_system_t, lspg_state_type>::fom_native_state_t;
  using fom_native_residual_t	= typename ExtractNativeHelper<fom_system_t, lspg_state_type>::fom_native_residual_t;

  // fom wrapper types
  using fom_state_t	= ::pressio::containers::Vector<fom_native_state_t>;
  using fom_residual_t	= ::pressio::containers::Vector<fom_native_residual_t>;

  // rom state type (passed in)
  using lspg_state_t		= lspg_state_type;

  // for LSPG, the rom residual type = containers::wrapper of application rhs
  // i.e. the wrapped fom rhs type
  using lspg_residual_t		= fom_residual_t;

  // decoder types (passed in)
  using decoder_t		= decoder_type;
  using decoder_jac_t = typename decoder_t::jacobian_type;

  // fro now, later on this needs to be detected from args
  using ud_ops_t = void;

  // fom state reconstructor type
  using fom_state_reconstr_t	= FomStateReconstructor<scalar_t, fom_state_t, decoder_t>;

  // class type holding fom states data: we only need to store one FOM state
  using fom_states_manager_t = ::pressio::rom::ManagerFomStatesStatic<fom_state_t, 1, fom_state_reconstr_t, void>;
};

}}}}}
#endif  // ROM_LSPG_IMPL_STEADY_TRAITS_ROM_LSPG_STEADY_COMMON_TRAITS_HPP_
