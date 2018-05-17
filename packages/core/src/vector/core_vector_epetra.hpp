
#ifndef CORE_VECTOR_EPETRA_HPP
#define CORE_VECTOR_EPETRA_HPP

#include "core_vectorBaseImpl.hpp"
#include "core_vectorMathImpl.hpp"
#include "core_vectorDistImpl.hpp"

#include "Epetra_Vector.h"


//*******************************
// epetra vector wrapper
//*******************************  


namespace core{
  
template <typename wrapped_type>
class vector<wrapped_type,
	     typename std::enable_if<std::is_same<wrapped_type,Epetra_Vector>::value>::type
	     >
  : public vectorBaseImpl<vector<Epetra_Vector> >,
    public vectorDistImpl<vector<Epetra_Vector> >,
    public vectorMathImpl<vector<Epetra_Vector> >
{
public:    
  using derived_t = vector<Epetra_Vector>;
  using sc_t = typename details::traits<derived_t>::scalar_t;
  using LO_t = typename details::traits<derived_t>::local_ordinal_t;
  using GO_t = typename details::traits<derived_t>::global_ordinal_t;
  using der_t = typename details::traits<derived_t>::derived_t;
  using wrap_t = typename details::traits<derived_t>::wrapped_t;
  using map_t = typename details::traits<derived_t>::map_t;
  using comm_t = typename details::traits<derived_t>::comm_t;

private:
  wrap_t data_;

public:
  vector(){};
  vector(const wrap_t & obj) : data_(obj){};
  ~vector(){};


  sc_t & operator [] (LO_t i){
    return data_[i];
  };
  sc_t const & operator [] (LO_t i) const{
    return data_[i];
  };  
  //-----------------------------------

  wrap_t const * viewImpl() const{
    return &data_;
  };
  wrap_t & getNonConstRefToDataImpl(){
    return data_;
  };
  //-----------------------------------

  sc_t dotImpl(const der_t & b) const{
    // what is this?
    // dot product: <this,b>
    sc_t res = 0.0;
    data_.Dot( *b.view(), &res );
    return res;    
  };

  template <typename op_t>
  void applyOpImpl(op_t op, sc_t a1,
  		   sc_t a2, const der_t & vin)
  {
    // static_assert(std::is_same<op_t,std::plus>::value,
    // 		  "This should be a +");
    // what is this?: this = a1*this op a2*vin;
    data_.Update(a2, *(vin.view()), a1);
  }
  
  sc_t norm2Impl() const{
    sc_t result = 0;
    data_.Norm2(&result);
    return result;
  };
  size_t localSizeImpl() const {
    return data_.MyLength();
  };
  size_t globalSizeImpl() const {
    return data_.GlobalLength();
  };
  map_t const & getMapImpl() const{
    return data_.Map();
  }

};

    
}//end namespace core
#endif

