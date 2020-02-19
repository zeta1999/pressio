/*
//@HEADER
// ************************************************************************
//
// containers_mvec_prod_mvec.hpp
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

#ifdef PRESSIO_ENABLE_TPL_TRILINOS
#ifndef CONTAINERS_SRC_OPS_EPETRA_MULTI_VECTOR_PROD_MULTI_VECTOR_HPP_
#define CONTAINERS_SRC_OPS_EPETRA_MULTI_VECTOR_PROD_MULTI_VECTOR_HPP_

namespace pressio{ namespace containers{ namespace ops{

/*
 * for epetra:
 *
 * C = beta * C + alpha*op(A)*op(B)
 *
*/


template <
  typename A_type, typename B_type, typename scalar_type, typename C_type
  >
::pressio::mpl::enable_if_t<
  ::pressio::containers::meta::is_multi_vector_wrapper_epetra<A_type>::value and
  ::pressio::containers::meta::is_multi_vector_wrapper_epetra<B_type>::value and
  ::pressio::containers::meta::is_dense_matrix_wrapper_eigen<C_type>::value
  >
product(::pressio::transpose modeA,
	::pressio::nontranspose modeB,
	const scalar_type alpha,
	const A_type & A,
	const B_type & B,
	const scalar_type beta,
	C_type & C)
{
  static_assert(containers::meta::are_scalar_compatible<A_type, B_type, C_type>::value,
		"Types are not scalar compatible");

  // how many vectors are in A and B
  const auto numVecsA = A.numVectors();
  const auto numVecsB = B.numVectors();
  assert(A.extent(0) == B.extent(0));
  assert(C.extent(0) == numVecsA);
  assert(C.extent(1) == numVecsB);

  auto const & Adata = *A.data();
  auto const & Bdata = *B.data();
  auto tmp = ::pressio::utils::constants::zero<scalar_type>();
  // compute dot between every column of A with every col of B
  for (size_t i=0; i<(size_t)numVecsA; i++){
    for (size_t j=0; j<(size_t)numVecsB; j++){
      C(i,j) = beta*C(i,j);
      Adata(i)->Dot( *(Bdata(j)), &tmp );
      C(i,j) += alpha*tmp;
    }
  }
}

template <
  typename C_type, typename A_type, typename B_type, typename scalar_type
  >
::pressio::mpl::enable_if_t<
  ::pressio::containers::meta::is_multi_vector_wrapper_epetra<A_type>::value and
  ::pressio::containers::meta::is_multi_vector_wrapper_epetra<B_type>::value and
  ::pressio::containers::meta::is_dynamic_dense_matrix_wrapper_eigen<C_type>::value,
  C_type
  >
product(::pressio::transpose modeA,
	::pressio::nontranspose modeB,
	const scalar_type alpha,
	const A_type & A,
	const B_type & B)
{
  static_assert(containers::meta::are_scalar_compatible<A_type, B_type, C_type>::value,
		"Types are not scalar compatible");
  constexpr auto zero = ::pressio::utils::constants::zero<scalar_type>();

  C_type C(A.numVectors(), B.numVectors());
  product(modeA, modeB, alpha, A, B, zero, C);
  return C;
}



/***********************************
 * special case A==B
**********************************/

template <
  typename A_type, typename scalar_type, typename C_type
  >
::pressio::mpl::enable_if_t<
  ::pressio::containers::meta::is_multi_vector_wrapper_epetra<A_type>::value and
  ::pressio::containers::meta::is_dense_matrix_wrapper_eigen<C_type>::value
  >
product(::pressio::transpose modeA,
	::pressio::nontranspose modeB,
	const scalar_type alpha,
	const A_type & A,
	const scalar_type beta,
	C_type & C)
{
  static_assert(containers::meta::are_scalar_compatible<A_type, C_type>::value,
		"Types are not scalar compatible");

  // how many vectors are in A and B
  const auto numVecsA = A.numVectors();
  assert(C.extent(0) == numVecsA);
  assert(C.extent(1) == numVecsA);
  auto const & Adata = *A.data();

  auto tmp = ::pressio::utils::constants::zero<scalar_type>();

  // A dot A = A^T*A, which yields a symmetric matrix
  // only need to compute half and fill remaining entries accordingly
  for (size_t i=0; i<(size_t)numVecsA; i++)
  {
    C(i,i) = beta*C(i,i);
    Adata(i)->Dot( *(Adata(i)), &tmp );
    C(i,i) += alpha*tmp;

    for (size_t j=i+1; j<(size_t)numVecsA; j++)
    {
      C(i,j) = beta*C(i,j);
      C(j,i) = beta*C(j,i);

      Adata(i)->Dot( *(Adata(j)), &tmp );
      C(i,j) += alpha*tmp;
      C(j,i) += alpha*tmp;
    }
  }
}

template <
  typename C_type, typename A_type, typename scalar_type
  >
::pressio::mpl::enable_if_t<
  ::pressio::containers::meta::is_multi_vector_wrapper_epetra<A_type>::value and
  ::pressio::containers::meta::is_dynamic_dense_matrix_wrapper_eigen<C_type>::value,
  C_type
  >
product(::pressio::transpose modeA,
	::pressio::nontranspose modeB,
	const scalar_type alpha,
	const A_type & A)
{
  static_assert(containers::meta::are_scalar_compatible<A_type, C_type>::value,
		"Types are not scalar compatible");

  constexpr auto zero = ::pressio::utils::constants::zero<scalar_type>();
  C_type C(A.numVectors(), A.numVectors());
  product(modeA, modeB, alpha, A, zero, C);
  return C;
}

}}}//end namespace pressio::containers::ops
#endif
#endif
