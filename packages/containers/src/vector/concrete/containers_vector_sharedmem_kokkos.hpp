/*
//@HEADER
// ************************************************************************
//
// containers_vector_sharedmem_kokkos.hpp
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

#ifndef CONTAINERS_VECTOR_CONCRETE_CONTAINERS_VECTOR_SHAREDMEM_KOKKOS_HPP_
#define CONTAINERS_VECTOR_CONCRETE_CONTAINERS_VECTOR_SHAREDMEM_KOKKOS_HPP_

namespace pressio{ namespace containers{

template <typename wrapped_type>
class Vector<wrapped_type,
	     ::pressio::mpl::enable_if_t<
	       containers::predicates::is_vector_kokkos<wrapped_type>::value
	       >
	     >
  : public VectorSharedMemBase< Vector<wrapped_type> >
{

  using this_t = Vector<wrapped_type>;
  using mytraits = typename details::traits<this_t>;
  using sc_t = typename mytraits::scalar_t;
  using ord_t = typename  mytraits::ordinal_t;
  using ref_t = typename mytraits::reference_t;
  using const_ref_t = typename mytraits::const_reference_t;

public:
  // default cnstr
  Vector() = default;

  // Views have "view semantics." copy constructor and
  // operator= only do shallow copies.
  // Here, for the time being, we construct wrapper
  // of a view WITHOUT doing shallow copy.
  // We create a new object and deep_copy original.

  explicit Vector(const wrapped_type src)
    : data_{src.label(), src.extent(0)}{
    Kokkos::deep_copy(data_, src);
  }

  Vector(const std::string & label, ord_t e1) : data_{label, e1}{}

  Vector(const ord_t e1) : data_{"empty", e1}{}

  // copy constructor implements copy semantics (for time being)
  Vector(const Vector & other)
    : data_{other.data_.label(), other.data_.extent(0)}{
    Kokkos::deep_copy(data_, other.data_);
  }

  // copy assign implments copy semantics not view (for time being)
  Vector & operator=(const Vector & other){
    if (&other != this){
      assert(this->extent(0) == other.extent(0));
      Kokkos::deep_copy(data_, *other.data());
    }
    return *this;
  }

  // move cnstr and assign
  Vector(Vector && other)
    : data_{other.data_.label(), other.data_.extent(0)}{
    Kokkos::deep_copy(data_, other.data_);
  }

  Vector & operator=(Vector && other){
    assert(this->extent(0) == other.extent(0));
    Kokkos::deep_copy(data_, *other.data());
    return *this;
  }

  ~Vector() = default;

  template<typename _wrapped_type = wrapped_type>
  mpl::enable_if_t<
    // todo: this is not entirely correct because this would work also
    // for UMV space, needs to be fixed
    std::is_same<typename mytraits::memory_space, Kokkos::HostSpace>::value,
    sc_t &>
  operator [] (ord_t i){
    return data_(i);
  };

  template<typename _wrapped_type = wrapped_type>
  mpl::enable_if_t<
    // todo: this is not entirely correct because this would work also
    // for UMV space, needs to be fixed
    std::is_same<typename mytraits::memory_space, Kokkos::HostSpace>::value,
    sc_t const &>
  operator [] (ord_t i) const{
    return data_(i);
  };

  template<typename _wrapped_type = wrapped_type>
  mpl::enable_if_t<
    // todo: this is not entirely correct because this would work also
    // for UMV space, needs to be fixed
    std::is_same<typename mytraits::memory_space, Kokkos::HostSpace>::value,
    sc_t &>
  operator () (ord_t i){
    return data_(i);
  };

  template<typename _wrapped_type = wrapped_type>
    mpl::enable_if_t<
      // todo: this is not entirely correct because this would work also
      // for UMV space, needs to be fixed
    std::is_same<typename mytraits::memory_space, Kokkos::HostSpace>::value,
    sc_t const &>
  operator () (ord_t i) const{
    return data_(i);
  };

  wrapped_type const * data() const{
    return &data_;
  }
  wrapped_type * data(){
    return &data_;
  }

  bool empty() const{
    return data_.extent(0)==0 ? true : false;
  }

  ord_t extent(ord_t i) const {
    assert( i == 0 );
    return data_.extent(i);
  }

private:
  friend VectorSharedMemBase< this_t >;
  wrapped_type data_ = {};

};//end class

}}//end namespace pressio::containers
#endif  // CONTAINERS_VECTOR_CONCRETE_CONTAINERS_VECTOR_SHAREDMEM_KOKKOS_HPP_
