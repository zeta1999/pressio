/*
//@HEADER
// ************************************************************************
//
// ode_aux_states_container.hpp
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

#ifndef ODE_AUX_STATES_CONTAINER_HPP_
#define ODE_AUX_STATES_CONTAINER_HPP_

#include "ode_fwd.hpp"
#include "../../containers/src/collection/containers_static_collection.hpp"

namespace pressio{ namespace ode{

template<bool is_explicit, typename T, std::size_t n>
class AuxStatesContainer;


// partially specialize for explicit scheme
template<typename T, std::size_t n>
class AuxStatesContainer<true, T, n>{
public:
  using data_type = ::pressio::containers::StaticCollection<T, n>;

  template <typename ... Args>
  AuxStatesContainer(Args && ... args)
    : data_( std::forward<Args>(args)... ){}

  ~AuxStatesContainer() = default;

public:
  static constexpr std::size_t size() {
    return data_type::size();
  }

  T & operator()(std::size_t i){
    assert( i<n );
    return data_(i);
  }

  T const & operator()(std::size_t i) const{
    assert( i<n );
    return data_(i);
  }

private:
  data_type data_;
};


// partially specialize for implicit scheme
template<typename T, std::size_t n>
class AuxStatesContainer<false, T, n>{
public:
  using data_type = ::pressio::containers::StaticCollection<T, n>;

  template <typename ... Args>
  AuxStatesContainer(Args && ... args)
    : data_( std::forward<Args>(args)... ){}

  ~AuxStatesContainer() = default;

public:
  static constexpr std::size_t size() {
    return data_type::size();
  }

  // n-1
  template <typename T1>
  ::pressio::mpl::enable_if_t< std::is_same<T1, ode::nMinusOne>::value, T &>
  get(){ return data_(0); }

  template <typename T1>
  ::pressio::mpl::enable_if_t< std::is_same<T1, ode::nMinusOne>::value, T const &>
  get() const{ return data_(0); }

  // n-2
  template <typename T1>
  ::pressio::mpl::enable_if_t< std::is_same<T1, ode::nMinusTwo>::value, T &>
  get(){ return data_(1); }

  template <typename T1>
  ::pressio::mpl::enable_if_t< std::is_same<T1, ode::nMinusTwo>::value, T const &>
  get() const{ return data_(1); }

  // n-3
  template <typename T1>
  ::pressio::mpl::enable_if_t< std::is_same<T1, ode::nMinusThree>::value, T &>
  get(){ return data_(1); }

  template <typename T1>
  ::pressio::mpl::enable_if_t< std::is_same<T1, ode::nMinusThree>::value, T const &>
  get() const{ return data_(1); }

  // n-4
  template <typename T1>
  ::pressio::mpl::enable_if_t< std::is_same<T1, ode::nMinusFour>::value, T &>
  get(){ return data_(1); }

  template <typename T1>
  ::pressio::mpl::enable_if_t< std::is_same<T1, ode::nMinusFour>::value, T const &>
  get() const{ return data_(1); }

private:
  data_type data_;
};

}}//end namespace pressio::ode
#endif