
#ifndef ODE_IMPLICIT_EULER_RESIDUAL_STANDARD_POLICY_HPP_
#define ODE_IMPLICIT_EULER_RESIDUAL_STANDARD_POLICY_HPP_

#include "ode_ConfigDefs.hpp"
#include "../base/ode_implicit_residual_policy_base.hpp"
#include "../../ode_residual_impl.hpp"

namespace ode{
namespace policy{

template<typename state_type,
	 typename residual_type,
	 typename model_type,
	 typename time_type,
	 typename sizer_type>
class implicitEulerResidualStandardPolicy
  : public implicitResidualPolicyBase<1,0,
  implicitEulerResidualStandardPolicy<state_type, residual_type,
				      model_type, time_type,
				      sizer_type>
  >
{
public:
  implicitEulerResidualStandardPolicy() = default;
  ~implicitEulerResidualStandardPolicy() = default;  

private:

  // enable when using types from core package
  template <typename U = state_type,
	    typename T = residual_type,
	    typename
	    std::enable_if<
	      core::meta::is_coreVector<U>::value==true &&
	      core::meta::is_coreVector<T>::value==true
	      >::type * = nullptr
	    >
  void computeImpl(const U & y,
		   T & R,
		   const std::array<U, 1> & oldYs,
		   model_type & model,
		   time_type t,
		   time_type dt)
  {
    if (R.empty())
      sizer_type::matchSize(y, R);

    R.setZero();
    model.residual(*y.data(), *R.data(), t);

    // do time discrete residual
    ode::impl::implicit_euler_time_discrete_residual(y, oldYs[0], R, dt);
  }
  //----------------------------------------------------------------

private:
  friend implicitResidualPolicyBase<1,0,
				    implicitEulerResidualStandardPolicy<
				      state_type, residual_type,
				      model_type, time_type,
				      sizer_type> >;
  
};//end class

}//end namespace polices
}//end namespace ode  
#endif 