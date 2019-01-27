
#ifndef ROM_LSPG_RESIDUAL_POLICY_HPP_
#define ROM_LSPG_RESIDUAL_POLICY_HPP_

#include "../rom_forward_declarations.hpp"
#include "../../../ode/src/implicit/ode_residual_impl.hpp"
#include "../../../ode/src/implicit/policies/base/ode_implicit_residual_policy_base.hpp"
#include "../rom_data_fom_rhs.hpp"
#include "../rom_data_fom_states.hpp"

namespace rompp{ namespace rom{

template <typename fom_states_data,
	  typename fom_rhs_data,
	  typename fom_eval_rhs_policy>
class LSPGResidualPolicy
  : public ode::policy::ImplicitResidualPolicyBase<
      LSPGResidualPolicy<fom_states_data,
			 fom_rhs_data,
			 fom_eval_rhs_policy>>,
    protected fom_states_data,
    protected fom_rhs_data,
    protected fom_eval_rhs_policy{

protected:
  using this_t = LSPGResidualPolicy<fom_states_data,
				       fom_rhs_data,
				       fom_eval_rhs_policy>;
  friend ode::policy::ImplicitResidualPolicyBase<this_t>;

  using fom_states_data::yFom_;
  using fom_states_data::yFomOld_;
  using fom_states_data::maxNstates_;
  using fom_rhs_data::fomRhs_;

public:
  static constexpr bool isResidualPolicy_ = true;
  using typename fom_rhs_data::fom_rhs_w_t;

public:
  LSPGResidualPolicy() = delete;
  ~LSPGResidualPolicy() = default;
  LSPGResidualPolicy(const fom_states_data & fomStates,
		     const fom_rhs_data & fomResids,
		     const fom_eval_rhs_policy & fomEvalRhsFunctor)
    : fom_states_data(fomStates),
      fom_rhs_data(fomResids),
      fom_eval_rhs_policy(fomEvalRhsFunctor){}

public:
  template <ode::ImplicitEnum odeMethod,
	    int n,
	    typename ode_state_t,
	    typename ode_residual_t,
	    typename app_t,
	    typename scalar_t>
  void operator()(const ode_state_t		  & odeY,
		  ode_residual_t		  & odeR,
  		  const std::array<ode_state_t,n> & oldYs,
  		  const app_t			  & app,
		  scalar_t			  t,
		  scalar_t			  dt) const
  {
    fom_states_data::template reconstructCurrentFomState(odeY);
    fom_states_data::template reconstructFomOldStates<n>(oldYs);
    fom_eval_rhs_policy::evaluate(app, yFom_, odeR, t);
    ode::impl::time_discrete_residual<
      odeMethod, maxNstates_
      >(yFom_, yFomOld_, odeR, dt);
  }

  template <ode::ImplicitEnum odeMethod,
	    int n,
	    typename ode_state_t,
	    typename app_t,
	    typename scalar_t>
  fom_rhs_w_t operator()(const ode_state_t		   & odeY,
			 const std::array<ode_state_t,n>   & oldYs,
			 const app_t			   & app,
			 scalar_t			   t,
			 scalar_t			   dt) const
  {
    (*this).template operator()<odeMethod, n>(odeY, fomRhs_,
					      oldYs, app, t, dt);
    return fomRhs_;
  }

};//end class

}}//end namespace rompp::rom
#endif
