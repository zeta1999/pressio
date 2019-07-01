
#ifdef HAVE_TRILINOS
#ifndef CONTAINERS_CONTAINER_OPS_MVEC_VEC_PROD_TPETRA_MULTI_VECTOR_PRODUCT_EIGEN_VECTOR_HPP_
#define CONTAINERS_CONTAINER_OPS_MVEC_VEC_PROD_TPETRA_MULTI_VECTOR_PRODUCT_EIGEN_VECTOR_HPP_

#include "../containers_ops_meta.hpp"
#include "../../vector/containers_vector_meta.hpp"
#include "../../multi_vector/containers_multi_vector_meta.hpp"
#include "../../vector/concrete/containers_vector_sharedmem_eigen_dynamic.hpp"

namespace rompp{ namespace containers{ namespace ops{

//--------------------------------------------------------
//Tpetra multivector with eigen vector, result is passed
template <typename mvec_type,
	  typename vec_type,
	  typename res_type,
  ::rompp::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper_tpetra<mvec_type>::value and
    containers::meta::wrapper_pair_have_same_scalar<mvec_type, vec_type>::value and
    containers::meta::is_vector_wrapper_eigen<vec_type>::value and
    containers::meta::is_vector_wrapper_tpetra<res_type>::value
    > * = nullptr
  >
void product(const mvec_type & mvA,
	     const vec_type & vecB,
	     res_type & C){

  //zero out result
  C.setZero();
  // how many vectors are in mvA
  auto numVecs = mvA.globalNumVectors();
  // size of vecB
  size_t vecBLen = vecB.size();
  if (vecBLen != size_t(numVecs))
    assert(size_t(numVecs) == vecBLen);

  // my number of rows
  auto myNrows = mvA.localLength();

  // get the wrapped trilinos tpetra multivector
  auto trilD = mvA.data();
  //  trilD->template sync<Kokkos::HostSpace>();
  auto mv2d = trilD->template getLocalView<Kokkos::HostSpace>();

  // get wrapped data for the result too
  auto C1 = C.data()->template getLocalView<Kokkos::HostSpace>();
  auto C2 = Kokkos::subview(C1, Kokkos::ALL(), 0);
  C.data()->template modify<Kokkos::HostSpace>();

  // loop
  for (decltype(myNrows) i=0; i<myNrows; i++){
    for (decltype(numVecs) j=0; j<numVecs; j++){
     C2[i] += mv2d(i,j) * vecB[j];
    }
  }
  using device_t = typename details::traits<res_type>::device_t;
  C.data()->template sync<device_t>();

}
//-------------------------------------------------------

// result is returned
template <typename mvec_type,
	  typename vec_type,
  ::rompp::mpl::enable_if_t<
   containers::meta::is_multi_vector_wrapper_tpetra<mvec_type>::value and
   containers::meta::wrapper_pair_have_same_scalar<mvec_type, vec_type>::value and
    (containers::meta::is_vector_wrapper_eigen<vec_type>::value)
  > * = nullptr
 >
auto product(const mvec_type & mvA, const vec_type & vecB)
  -> containers::Vector<
  Tpetra::Vector<typename details::traits<mvec_type>::scalar_t,
                 typename details::traits<mvec_type>::local_ordinal_t,
                 typename details::traits<mvec_type>::global_ordinal_t,
                 typename details::traits<mvec_type>::node_t>
                 >
{

  // here, mvA is distrubted, but vecB is NOT.
  // we interpret this as a linear combination of vectors

  // the data map of the multivector
  auto rcpMap = mvA.getRCPDataMap();

  using mvec_traits = typename details::traits<mvec_type>;
  using sc_t = typename mvec_traits::scalar_t;
  using LO_t = typename mvec_traits::local_ordinal_t;
  using GO_t = typename mvec_traits::global_ordinal_t;
  using NO_t = typename mvec_traits::node_t;

  // result is an Tpetra Vector with same distribution of mvA
  using res_nat_t = Tpetra::Vector<sc_t, LO_t, GO_t, NO_t>;
  //  res_nat_t tmp(rcpMap);
  using res_t = containers::Vector<res_nat_t>;
  res_t c(rcpMap);
  product(mvA, vecB, c);
  return c;
}


}}}//end namespace rompp::containers::ops
#endif
#endif //HAVE_TRILINOS