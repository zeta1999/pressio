
#ifndef ROM_LSPG_TYPE_GENERATOR_PRECONDITIONED_HPP_
#define ROM_LSPG_TYPE_GENERATOR_PRECONDITIONED_HPP_

#include "rom_lspg_type_generator_common.hpp"

namespace rompp{ namespace rom{

template <typename fom_type,
	  ode::ImplicitEnum odeName,
	  typename decoder_type,
	  typename lspg_state_type>
struct PreconditionedLSPGTypeGenerator
  : LSPGCommonTypes<
  fom_type, odeName, decoder_type, lspg_state_type
  >{

  using base_t = LSPGCommonTypes<
    fom_type, odeName, decoder_type, lspg_state_type>;

  using typename base_t::fom_t;
  using typename base_t::scalar_t;
  using typename base_t::fom_state_t;
  using typename base_t::fom_state_w_t;
  using typename base_t::fom_rhs_w_t;
  using typename base_t::decoder_t;
  using typename base_t::decoder_jac_t;
  using typename base_t::lspg_state_t;
  using typename base_t::lspg_residual_t;
  using typename base_t::fom_states_data;
  using typename base_t::fom_rhs_data;

  /* lspg_matrix_t is type of J*decoder_jac_t (in the most basic case) where
   * * J is the jacobian of the fom rhs
   * * decoder_jac_t is the type of the decoder jacobian
   * In more complex cases, we might have (something)*J*decoder_jac_t,
   * where (something) is product of few matrices.
   * For now, set lspg_matrix_t to be of same type as decoder_jac_t
   * if phi is MV<>, then lspg_matrix_t = core::MV<>
   * if phi is Matrix<>, then we have core::Matrix<>
   * not a bad assumption since all matrices are left-applied to decoder_jac_t
   */
  using lspg_matrix_t		= decoder_jac_t;

  // policy for evaluating the rhs of the fom object
  using fom_eval_rhs_policy_t	= rom::policy::EvaluateFomRhsDefault;

  // policy for left multiplying the fom jacobian with decoder_jac_t
  // possibly involving other stuff like explained above
  using fom_apply_jac_policy_t	= rom::policy::ApplyFomJacobianDefault;

  // policy defining how to compute the LSPG time-discrete residual
  using lspg_residual_policy_t =
    rom::decorator::Preconditioned<
    rom::LSPGResidualPolicy<
      fom_states_data, fom_rhs_data, fom_eval_rhs_policy_t
      >
    >;

  // policy defining how to compute the LSPG time-discrete jacobian
  using lspg_jacobian_policy_t	=
    rom::decorator::Preconditioned<
    rom::LSPGJacobianPolicy<
      fom_states_data, lspg_matrix_t, fom_apply_jac_policy_t
      >
    >;

  // auxiliary stepper
  using aux_stepper_t = typename auxStepperHelper<
    odeName, lspg_state_type,
    lspg_residual_t, lspg_matrix_t,
    fom_type, lspg_residual_policy_t,
    lspg_jacobian_policy_t>::type;

  // primary stepper type
  using rom_stepper_t		= ode::ImplicitStepper<
    odeName, lspg_state_type,
    lspg_residual_t, lspg_matrix_t,
    fom_type, aux_stepper_t,
    lspg_residual_policy_t, lspg_jacobian_policy_t>;

};//end class

}}//end  namespace rompp::rom
#endif