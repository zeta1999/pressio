
#ifndef ODE_POLICIES_STANDARD_IMPLICIT_EULER_RESIDUAL_STANDARD_POLICY_HPP_
#define ODE_POLICIES_STANDARD_IMPLICIT_EULER_RESIDUAL_STANDARD_POLICY_HPP_

#include "ode_ConfigDefs.hpp"
#include "../base/ode_implicit_residual_policy_base.hpp"
#include "../../ode_residual_impl.hpp"

namespace ode{
namespace policy{

template<typename state_type,
	 typename residual_type,
	 typename model_type>
class ImplicitEulerResidualStandardPolicy
  : public ImplicitResidualPolicyBase<
  ImplicitEulerResidualStandardPolicy<state_type, residual_type,
					  model_type>, 1, 0 >
{
public:
  ImplicitEulerResidualStandardPolicy() = default;
  ~ImplicitEulerResidualStandardPolicy() = default;  

private:
  using scalar_type = typename core::details::traits<state_type>::scalar_t;
  
private:

  template <typename U = state_type,
	    typename T = residual_type,
	    typename
	    std::enable_if<
	      core::meta::is_core_vector_wrapper<U>::value==true &&
	      core::meta::is_core_vector_wrapper<T>::value==true
	      >::type * = nullptr
	    >
  void computeImpl(const U & y,
		   T & R,
		   const std::array<U, 1> & oldYs,
		   model_type & model,
		   scalar_type t,
		   scalar_type dt){
    if (R.empty())
      R.matchLayoutWith(y);

    R.setZero();
    model.residual(*y.data(), *R.data(), t);

    // do time discrete residual
    ode::impl::implicit_euler_time_discrete_residual(y, oldYs[0], R, dt);
  }
  //----------------------------------------------------------------

private:
  friend ImplicitResidualPolicyBase<
				    ImplicitEulerResidualStandardPolicy<
				      state_type, residual_type,
				      model_type>, 1,0>;
};//end class

}//end namespace polices
}//end namespace ode  
#endif 