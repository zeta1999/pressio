/*
//@HEADER
// ************************************************************************
//
// rom_query_fom_velocity.hpp
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

#ifndef ROM_FOM_QUERY_ROM_QUERY_FOM_VELOCITY_HPP_
#define ROM_FOM_QUERY_ROM_QUERY_FOM_VELOCITY_HPP_

namespace pressio{ namespace rom{

//------------------------------------------
// for native c++
//------------------------------------------
template <typename fom_system_t, typename state_t, typename rhs_t, typename time_t>
mpl::enable_if_t<
  !::pressio::ops::predicates::is_object_pybind<fom_system_t>::value
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
  and !::pressio::containers::predicates::is_vector_wrapper_pybind<state_t>::value
  and !::pressio::containers::predicates::is_vector_wrapper_pybind<rhs_t>::value
#endif
  >
queryFomVelocity(const fom_system_t & fomObj,
		 const state_t & fomState,
		 rhs_t & rhs,
		 const time_t & time)
{
  fomObj.velocity(*fomState.data(), time, *rhs.data());
}

//------------------------------------------
// for python
//------------------------------------------
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
template <typename state_t, typename rhs_t, typename time_t>
mpl::enable_if_t<
  ::pressio::containers::predicates::is_vector_wrapper_pybind<state_t>::value and
  ::pressio::containers::predicates::is_vector_wrapper_pybind<rhs_t>::value
>
queryFomVelocity(const pybind11::object & fomObj,
		 const state_t & fomState,
		 rhs_t & rhs,
		 const time_t & time)
{
  *rhs.data() = fomObj.attr("velocity")(*fomState.data(), time);
}
#endif

}} //end namespace pressio::rom
#endif  // ROM_FOM_QUERY_ROM_QUERY_FOM_VELOCITY_HPP_
