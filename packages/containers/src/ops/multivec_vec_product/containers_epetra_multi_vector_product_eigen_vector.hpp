
#ifdef HAVE_TRILINOS
#ifndef CONTAINERS_CONTAINER_OPS_MVEC_VEC_PROD_EPETRA_MULTI_VECTOR_PRODUCT_EIGEN_VECTOR_HPP_
#define CONTAINERS_CONTAINER_OPS_MVEC_VEC_PROD_EPETRA_MULTI_VECTOR_PRODUCT_EIGEN_VECTOR_HPP_

#include "../containers_ops_meta.hpp"
#include "../../vector/containers_vector_meta.hpp"
#include "../../multi_vector/containers_multi_vector_meta.hpp"
#include "../../vector/concrete/containers_vector_sharedmem_eigen_dynamic.hpp"

namespace rompp{ namespace containers{ namespace ops{

//  Epetra multivector with eigen or armadillo vector
// we pass the result object
template <typename mvec_type,
	  typename vec_type,
  ::rompp::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper_epetra<mvec_type>::value and
    containers::meta::wrapper_pair_have_same_scalar<mvec_type, vec_type>::value and
    containers::meta::is_vector_wrapper_eigen<vec_type>::value
    > * = nullptr
  >
void product(const mvec_type & mvA,
	     const vec_type & vecB,
	     containers::Vector<Epetra_Vector> & C){

  //zero out result
  C.setZero();
  // how many vectors are in mvA
  auto numVecs = mvA.globalNumVectors();
  // size of vecB
  assert(size_t(numVecs) == size_t(vecB.size()));
  // the data map of the multivector
  auto mvMap = mvA.getDataMap();
  // my number of rows
  auto myNrows = mvMap.NumMyElements();

  // loop
  for (decltype(myNrows) i=0; i<myNrows; i++){
    for (decltype(numVecs) j=0; j<numVecs; j++){
      C[i] += mvA(i,j) * vecB[j];
    }
  }
}
//-------------------------------------------------------

// result is returned
template <typename mvec_type,
	  typename vec_type,
  ::rompp::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper_epetra<mvec_type>::value and
    containers::meta::wrapper_pair_have_same_scalar<mvec_type, vec_type>::value and
    containers::meta::is_vector_wrapper_eigen<vec_type>::value
  > * = nullptr
 >
containers::Vector<Epetra_Vector>
product(const mvec_type & mvA, const vec_type & vecB) {

  // here, mvA is distrubted, but vecB is NOT.
  // we interpret this as a linear combination of vectors

  // the data map of the multivector
  auto mvMap = mvA.getDataMap();
  // result is an Epetra Vector with same distribution of mvA
  using res_t = containers::Vector<Epetra_Vector>;
  res_t c(mvMap);
  product(mvA, vecB, c);
  return c;
}

}}}//end namespace rompp::containers::ops
#endif
#endif //HAVE_TRILINOS