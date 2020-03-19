/*
//@HEADER
// ************************************************************************
//
// ops_mvec_prod_mvec.hpp
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
#ifndef OPS_SRC_OPS_TPETRA_BLOCK_MULTI_VECTOR_PROD_MULTI_VECTOR_HPP_
#define OPS_SRC_OPS_TPETRA_BLOCK_MULTI_VECTOR_PROD_MULTI_VECTOR_HPP_

namespace pressio{ namespace ops{


/*
 * for tpetra:
 *
 * C = beta * C + alpha*op(A)*op(B)
 *
*/

/* -------------------------------------------------------------------
 * specialize for op(A) = A^T and op(B) = B
 *-------------------------------------------------------------------*/
template <
  typename A_type, typename B_type, typename scalar_type, typename C_type
  >
::pressio::mpl::enable_if_t<
  ::pressio::containers::meta::is_multi_vector_wrapper_tpetra_block<A_type>::value and
  ::pressio::containers::meta::is_multi_vector_wrapper_tpetra_block<B_type>::value
  >
product(::pressio::transpose modeA,
	::pressio::nontranspose modeB,
	const scalar_type alpha,
	const A_type & A,
	const B_type & B,
	const scalar_type beta,
	::pressio::containers::MatrixSharedMemBase<C_type> & C)
{
  throw std::runtime_error("Error, C = beta*C + alpha*A^T*B for tpetra block not yet supported");

  // static_assert(containers::meta::are_scalar_compatible<A_type, B_type, C_type>::value,
  // 		"Types are not scalar compatible");
//   // how many vectors are in mvA and mvB
//   const auto numVecsA = mvA.globalNumVectors();
//   const auto numVecsB = mvB.globalNumVectors();
//   auto mvA_v = mvA.data()->getMultiVectorView();
//   auto mvB_v = mvB.data()->getMultiVectorView();
//   // compute dot between every column of A with every col of B
//   for (std::size_t i=0; i<(std::size_t)numVecsA; i++)
//   {
//     // colI is a Teuchos::RCP<Vector<...>>
//     const auto colI = mvA_v.getVector(i);
//     for (std::size_t j=0; j<(std::size_t)numVecsB; j++)
//     {
//       const auto colJ = mvB_v.getVector(j);
//       C(i,j) = colI->dot(*colJ);
//     }
//   }
}

template <
  typename C_type, typename A_type, typename B_type, typename scalar_type
  >
::pressio::mpl::enable_if_t<
  ::pressio::containers::meta::is_multi_vector_wrapper_tpetra_block<A_type>::value and
  ::pressio::containers::meta::is_multi_vector_wrapper_tpetra_block<B_type>::value and
  (::pressio::containers::meta::is_dense_matrix_wrapper_eigen<C_type>::value or
   ::pressio::containers::meta::is_dense_matrix_wrapper_kokkos<C_type>::value),
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

  const auto numVecsA = A.numVectors();
  const auto numVecsB = B.numVectors();
  C_type C(numVecsA, numVecsB);
  product(modeA, modeB, alpha, A, B, zero, C);
  return C;
}



// /***********************************
//  * special case A==B
// **********************************/
template <typename A_type, typename scalar_type, typename C_type>
::pressio::mpl::enable_if_t<
  ::pressio::containers::meta::is_multi_vector_wrapper_tpetra_block<A_type>::value and
  ::pressio::containers::meta::is_dense_matrix_wrapper_eigen<C_type>::value
  >
product(::pressio::transpose modeA,
	::pressio::nontranspose modeB,
	const scalar_type alpha,
	const A_type & A,
	const scalar_type beta,
	C_type & C)
{
  throw std::runtime_error("Error, C = beta*C + alpha*A^T*A for tpetra block not yet supported");

  // static_assert(containers::meta::are_scalar_compatible<A_type, C_type>::value,
  // 		"Types are not scalar compatible");

  // // get a tpetra multivector that views the data
  // const auto mvView = mvA.data()->getMultiVectorView();

  // // how many vectors are in mvA and mvB
  // const auto numVecsA = mvA.globalNumVectors();

  // // A dot A = A^T*A, which yields a symmetric matrix
  // // only need to compute half and fill remaining entries accordingly
  // for (std::size_t i=0; i<(std::size_t)numVecsA; i++)
  // {
  //   // colI is a Teuchos::RCP<Vector<...>>
  //   const auto colI = mvView.getVector(i);
  //   for (std::size_t j=i; j<(std::size_t)numVecsA; j++)
  //   {
  //     const auto colJ = mvView.getVector(j);
  //     C(i,j) = colI->dot(*colJ);
  //     C(j,i) = C(i,j);
  //   }
  // }
}


template <typename A_type, typename scalar_type, typename C_type>
::pressio::mpl::enable_if_t<
  ::pressio::containers::meta::is_multi_vector_wrapper_tpetra_block<A_type>::value and
  ::pressio::containers::meta::is_dense_matrix_wrapper_kokkos<C_type>::value
  >
product(::pressio::transpose modeA,
	::pressio::nontranspose modeB,
	const scalar_type alpha,
	const A_type & A,
	const scalar_type beta,
	C_type & C)
{
  throw std::runtime_error("Error, C = beta*C + alpha*A^T*A for tpetra block not yet supported");

  // static_assert(containers::meta::are_scalar_compatible<A_type, C_type>::value,
  // 		"Types are not scalar compatible");

  // // check traits of the block mv
  // using scalar_t     = typename ::pressio::containers::details::traits<mvec_t>::scalar_t;
  // using tpetra_mvb_t = typename ::pressio::containers::details::traits<mvec_t>::wrapped_t;

  // // from the mvb type we can get the underlying regular tpetra::mv
  // using tpetra_mv_t  = typename tpetra_mvb_t::mv_type;
  // using map_t	     = typename tpetra_mv_t::map_type;

  // // get a tpetra multivector that views the tpetra block data
  // const auto mvView = A.data()->getMultiVectorView();

  // const auto indexBase = mvView.getMap()->getIndexBase();
  // const auto comm = mvView.getMap()->getComm();
  // // C should be symmetric
  // assert( C.rows() == C.cols() );
  // const auto n = C.rows();
  // Teuchos::RCP<const map_t> replMap(new map_t(n, indexBase, comm, Tpetra::LocallyReplicated));
  // // create multivector that views the Kokkos matrix result
  // tpetra_mv_t Cmv(replMap, *C.data());

  // constexpr auto beta  = ::pressio::utils::constants::zero<scalar_t>();
  // constexpr auto alpha = ::pressio::utils::constants::one<scalar_t>();
  // // do the operation C = A^T A
  // Cmv.multiply(Teuchos::ETransp::TRANS, Teuchos::ETransp::NO_TRANS, alpha, mvView, mvView, beta);
}

template <typename C_type, typename A_type, typename scalar_type>
::pressio::mpl::enable_if_t<
  ::pressio::containers::meta::is_multi_vector_wrapper_tpetra_block<A_type>::value and
  (::pressio::containers::meta::is_dynamic_dense_matrix_wrapper_eigen<C_type>::value or
   ::pressio::containers::meta::is_dense_matrix_wrapper_kokkos<C_type>::value),
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

}}//end namespace pressio::ops
#endif
#endif