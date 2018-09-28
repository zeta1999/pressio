
#ifndef ODE_INCREMENTBASED_RESIDUAL_HPP_
#define ODE_INCREMENTBASED_RESIDUAL_HPP_

#include "ode_ConfigDefs.hpp"
#include "../base/ode_residual_policy_base.hpp"
#include "../base/ode_advance_increment_policy_base.hpp"

namespace rompp{
namespace ode{
namespace policy{

template<typename state_type,
	 typename residual_type,
	 typename model_type,
	 typename time_type,
	 typename sizer_type>
class incrementBasedResidual
  : public residualPolicyBase<incrementBasedResidual,
			      state_type,
			      residual_type,
			      model_type,
			      time_type,
			      sizer_type>,
  public advanceIncrementPolicyBase<incrementBasedResidual,
				    state_type, residual_type,
				    model_type, time_type,
				    sizer_type>
{
private:
  using baseIncr_t = advanceIncrementPolicyBase<incrementBasedResidual,
						state_type, residual_type,
						model_type, time_type, sizer_type>;
public:
  incrementBasedResidual(const state_type & y0)
    : baseIncr_t(y0){}
  ~incrementBasedResidual() = default;  

private:
  using baseIncr_t::yFull_;
  using baseIncr_t::y0ptr_;

private:
  // enable if using types from core package
  template <typename U = state_type, typename T = residual_type,
	    typename std::enable_if<
	      core::meta::is_coreVector<U>::value==true &&
	      core::meta::is_coreVector<T>::value==true
	    >::type * = nullptr>
  void computeImpl(const U & y, 
        T & R,
        model_type & model, 
        time_type t)
  { 
    // reconstruct the solution
    yFull_ = *y0ptr_ + y;

    // eval RHS from target model
    R.setZero();
    model.residual(*yFull_.data(), *R.data(), t);
  }  

  void weightTimeDiscreteResidualImpl(const state_type & y,
				      residual_type & R,
				      model_type & model, 
				      time_type t)
  {
    //no op
  }
  
private:
  friend residualPolicyBase<incrementBasedResidual,
			    state_type, residual_type,
			    model_type, time_type, sizer_type>;

  friend baseIncr_t;

};//end class

}//end namespace polices
}//end namespace ode  
}//end namespace rompp
#endif 