/*
//@HEADER
// ************************************************************************
//
// rom_compose_galerkin_impl.hpp
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

#ifndef ROM_GALERKIN_IMPL_ROM_COMPOSE_GALERKIN_IMPL_HPP_
#define ROM_GALERKIN_IMPL_ROM_COMPOSE_GALERKIN_IMPL_HPP_

#include "./continuous_time_api/rom_galerkin_problem_continuous_time_api.hpp"
#include "./discrete_time_api/rom_galerkin_problem_discrete_time_api.hpp"

namespace pressio{ namespace rom{ namespace galerkin{ namespace impl{

struct Default{};

template<typename tag, typename ...Args>
struct compose{
  using type = void;
};

// default continuous time API
template<typename stepper_tag, typename fom_system_type, typename galerkin_state_t, typename ...Args>
struct compose<
::pressio::rom::galerkin::impl::Default,
mpl::enable_if_t<
::pressio::rom::concepts::continuous_time_explicit_system<fom_system_type>::value and
(std::is_same< stepper_tag, ::pressio::ode::explicitmethods::Euler>::value or
std::is_same< stepper_tag, ::pressio::ode::explicitmethods::RungeKutta4>::value)
>,
stepper_tag, fom_system_type, galerkin_state_t, Args...>
{
  using type = ::pressio::rom::galerkin::impl::ProblemContinuousTimeApi<
            ::pressio::rom::galerkin::impl::DefaultProblemTraitsContinuousTimeApi,
            stepper_tag, fom_system_type, galerkin_state_t, Args...>;
};

// default discrete time api
template<typename stepper_tag, typename fom_system_type, typename galerkin_state_t, typename ...Args>
struct compose<
::pressio::rom::galerkin::impl::Default,
mpl::enable_if_t<
::pressio::rom::concepts::discrete_time_system<fom_system_type>::value and
std::is_same< stepper_tag, ::pressio::ode::implicitmethods::Arbitrary>::value
>,
stepper_tag, fom_system_type, galerkin_state_t, Args...>
{
  using type = ::pressio::rom::galerkin::impl::ProblemDiscreteTimeApi<
            ::pressio::rom::galerkin::impl::DefaultProblemTraitsDiscreteTimeApi,
            stepper_tag, fom_system_type, galerkin_state_t, Args...>;
};

}}}}
#endif  // ROM_GALERKIN_IMPL_ROM_COMPOSE_GALERKIN_IMPL_HPP_
