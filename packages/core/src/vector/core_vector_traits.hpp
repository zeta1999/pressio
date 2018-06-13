
#ifndef CORE_VECTOR_TRAITS_HPP_
#define CORE_VECTOR_TRAITS_HPP_

#include "core_forward_declarations.hpp"
#include "core_vector_meta.hpp"
#include <vector>
#include "Epetra_Vector.h"
#include "Eigen/Dense"

namespace core{
namespace details{

  //*******************************
  // Eigen vector 
  //******************************* 
  template <typename wrapped_type>
  struct traits<vector<wrapped_type,
  		       typename
		       std::enable_if<
			 core::meta::is_vectorEigen<wrapped_type
						    >::value
			 >::type
  		       >
  		>     		       
  {
    using scalar_t = typename wrapped_type::Scalar;
    using ordinal_t = int;
    using wrapped_t = wrapped_type;
    using derived_t = vector<wrapped_t>;
    static constexpr int isVector = 1;
    static constexpr int isEigen = 1;
    static constexpr int isSerial = 1;
    static constexpr int isSTDVector = 0;
    static constexpr int isDistributed = 0;
    // make these void just to be clear they are not usable
    using local_ordinal_t = void;
    using global_ordinal_t = void;
  };
  
     
  //*******************************
  // for a std vector 
  //******************************* 
  template <typename wrapped_type>
  struct traits<vector<wrapped_type,
		       typename
		       std::enable_if<
			 core::meta::is_vectorStdLib<wrapped_type
						     >::value
			 >::type
		       >
		>
  {
    using scalar_t = typename wrapped_type::value_type;
    using ordinal_t = core::defaultTypes::local_ordinal_t;
    using wrapped_t = wrapped_type;
    using derived_t = vector<wrapped_t>;
    static constexpr int isVector = 1;
    static constexpr int isSTDVector = 1;
    static constexpr int isSerial = 1;
    static constexpr int isDistributed = 0;
    static constexpr int isEigen = 0;
    // make these void just to be clear they are not usable
    using local_ordinal_t = void;
    using global_ordinal_t = void;
  };

  
  //*******************************
  // user-defined serial vector 
  //******************************* 
  template <typename wrapped_type>
  struct traits<vector<wrapped_type,
  		       typename
  		       std::enable_if<
			 core::meta::has_scalarTypedef<wrapped_type
						       >::value &&
			 core::meta::has_ordinalTypedef<wrapped_type
							>::value &&
			 !core::meta::is_vectorStdLib<wrapped_type
						      >::value &&
			 !core::meta::is_vectorEigen<wrapped_type
						     >::value
			 >::type
		       >
  		>
  {
    using scalar_t = typename wrapped_type::scalar_type;
    using ordinal_t = typename wrapped_type::ordinal_type;
    using wrapped_t = wrapped_type;
    using derived_t = vector<wrapped_type>;
    static constexpr int isVector = 1;
    static constexpr int isSerial = 1;
    static constexpr int isDistributed = 0;
    static constexpr int isEigen = 0;
    static constexpr int isSTDVector = 0;
    // make these void just to be clear they are not usable
    using local_ordinal_t = void;
    using global_ordinal_t = void;
  };


  //*******************************
  // for epetra vector 
  //******************************* 
  template<typename wrapped_type>
  struct traits<vector<wrapped_type,
		       typename std::enable_if<
			 std::is_same<wrapped_type,
				      Epetra_Vector
				      >::value
			 >::type
		       >
		>
  {
    using scalar_t = defaultTypes::epetra_scalar_t;
    using local_ordinal_t = core::defaultTypes::epetra_lo_t;
    using global_ordinal_t = core::defaultTypes::epetra_go_t1;
    using wrapped_t = Epetra_Vector;
    using data_map_t = Epetra_BlockMap;
    using communicator_t = Epetra_Comm;
    using derived_t = vector<wrapped_t>;
    static constexpr int isVector = 1;
    static constexpr int isDistributed = 1;
    static constexpr int isEpetra = 1;
    static constexpr int isSTDVector = 0;
    static constexpr int isSerial = 0;
    static constexpr int isEigen = 0;
    // make these void just to be clear they are not usable
    using ordinal_t = void;
  };

  
}//end namespace details  
}//end namespace core
#endif







  // //*******************************
  // // for a general distributed vector 
  // //******************************* 
  // template <typename wrapped_type>
  // struct traits<vector<wrapped_type,
  // 		       typename std::enable_if<!std::is_same<wrapped_type,Epetra_Vector>::value &&
  // 					       core::meta::has_scalarTypedef<wrapped_type>::value &&
  // 					       core::meta::has_localOrdinalTypedef<wrapped_type>::value &&
  // 					       core::meta::has_globalOrdinalTypedef<wrapped_type>::value &&
  // 					       core::meta::has_mapTypedef<wrapped_type>::value && 
  // 					       core::meta::has_commTypedef<wrapped_type>::value
  // 					       >::type
  // 		       >
  // 		>
  // {
  //   using scalar_t = typename wrapped_type::scalar_type;
  //   using local_ordinal_t = typename wrapped_type::local_ordinal_type;
  //   using global_ordinal_t = typename wrapped_type::global_ordinal_type;
  //   using wrapped_t = wrapped_type;
  //   using derived_t = vector<wrapped_t>;
  //   using map_t = typename wrapped_type::map_type;
  //   using comm_t = typename wrapped_type::comm_type;
  //   enum {
  //     isVector = 1,
  //     isSTDVector = 0,
  //     isSerial = 0,
  //     isDistributed = 1,
  //     isEigen = 0
  //   };
  // };