/*
//@HEADER
// ************************************************************************
//
//                          pop_front.hpp
//                         dharma_new
//              Copyright (C) 2017 NTESS, LLC
//
// Under the terms of Contract DE-NA-0003525 with NTESS, LLC,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact somebody@sandia.gov
//
// ************************************************************************
//@HEADER
*/

#ifndef FRONTEND_INCLUDE_TINYMPL_POP_FRONT_HPP_
#define FRONTEND_INCLUDE_TINYMPL_POP_FRONT_HPP_

#include <tinympl/sequence.hpp>
#include <tinympl/as_sequence.hpp>
#include <tinympl/variadic/pop_front.hpp>

namespace tinympl {

template<class Sequence,
  template <class...> class Out = as_sequence<Sequence>::template rebind
>
struct pop_front : pop_front<as_sequence_t<Sequence>, Out> {};

template<template <class...> class Out, class... Args>
struct pop_front<sequence<Args...>, Out> : variadic::pop_front<Out, Args...> {};

namespace types_only {

template<class Sequence>
struct pop_front : tinympl::pop_front<Sequence, as_sequence<Sequence>::template rebind> { };

} // end namespace types_only

} // namespace tinympl



#endif /* FRONTEND_INCLUDE_TINYMPL_POP_FRONT_HPP_ */