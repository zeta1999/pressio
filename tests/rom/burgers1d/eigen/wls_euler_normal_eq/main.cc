#include <iostream>
#include <string>
#include "CONTAINERS_ALL"
#include "ODE_ALL"
#include "SOLVERS_NONLINEAR"
#include "ROM_LSPG_UNSTEADY"
#include "APPS_UNSTEADYBURGERS1D"
#include "utils_eigen.hpp"
#include "wls_apis.hpp"
//#include "/Users/ejparis/pressio_repos/pressio/packages/solvers/src/meta/solvers_basic_meta.hpp"
//#include "/Users/ejparis/pressio_repos/pressio/packages/solvers/src/meta/solvers_system_has_all_needed_jacobian_methods.hpp"
//#include "/Users/ejparis/pressio_repos/pressio/packages/solvers/src/meta/solvers_system_has_all_needed_residual_methods.hpp"
//#include "/Users/ejparis/pressio_repos/pressio/packages/rom/src/experimental/wls/rom_wls_sequential_residual_policy_for_residual_api.hpp"
//#include "/Users/ejparis/pressio_repos/pressio/packages/rom/src/experimental/wls/rom_wls_draft.hpp"

// n is the window size

namespace pressio{ namespace rom{ namespace wls{


template<typename api_tag, typename fom_type, typename wls_state_type,typename ... rest>
struct Specializer{
  using type = void;
};

//struct DefaultApi{};

template<
  typename fom_type, typename wls_state_type,
  typename residual_pol_type, typename jac_pol_type
>
struct Specializer<pressio::rom::wls::WlsSystemDefaultApi<fom_type,wls_state_type,residual_pol_type,jac_pol_type>, fom_type, wls_state_type, residual_pol_type, jac_pol_type>
{ 
//  static_assert( ::pressio::rom::meta::is_legitimate_residual_policy_for_wls<residual_pol_t>,
//     "You are trying to use Wls with default api but the residual policy passed \
is not admissible for this: maybe you have the wrong api? blas blas");
  using type = pressio::rom::wls::WlsSystemDefaultApi<fom_type, wls_state_type, residual_pol_type, jac_pol_type>;
};
template<
  typename fom_type, typename wls_state_type,
  typename jtj_jtr_pol_type
>
struct Specializer<pressio::rom::wls::WlsSystemJTJApi<fom_type,wls_state_type,jtj_jtr_pol_type>, fom_type, wls_state_type, jtj_jtr_pol_type>
{ 
//  static_assert( ::pressio::rom::meta::is_legitimate_residual_policy_for_wls<residual_pol_t>,
//     "You are trying to use Wls with default api but the residual policy passed \
is not admissible for this: maybe you have the wrong api? blas blas");
  using type = pressio::rom::wls::WlsSystemJTJApi<fom_type, wls_state_type,jtj_jtr_pol_type>;
};

//struct DefaultApi{};
//struct JTJRApi{};
template<typename api_type,typename fom_type,typename wls_state_type,typename ... rest>
using WlsSystem = typename pressio::rom::wls::Specializer<api_type, fom_type, wls_state_type, rest...>::type;
}}} // end namespace rom::wls



/*

namespace pressio{ namespace rom{ namespace experimental{

template <std::size_t numAuxStates, typename residual_type, typename fom_states_container_type>
class WlsSequentialResidualPolicyForResidualApi
{
public:
  using this_t = WlsSequentialResidualPolicyForResidualApi<numAuxStates, residual_type, fom_states_container_type>;
  using residual_t = residual_type;

public:
  WlsSequentialResidualPolicyForResidualApi() = delete;
  ~WlsSequentialResidualPolicyForResidualApi() = default;

  WlsSequentialResidualPolicyForResidualApi(fom_states_container_type & fomStatesIn)
    : fomStates_(fomStatesIn){}

public:
  template <typename wls_state_t, typename fom_t>
  void operator()(const wls_state_t                     & wlsState,
                  const fom_t                           & fomObj,
                  residual_t                            & wlsR) const
  {
    this->compute_impl(wlsState, fomObj, wlsR);
  }

  template <typename wls_state_t, typename fom_t>
  residual_t operator()(const wls_state_t                   & wlsState,
                        const fom_t                         & fomObj) const
  {
    fomStates_.template reconstructFomStateAt(wlsState, 0);
    residual_t R( fomObj.createTimeDiscreteResidualObject( *fomStates_[0].data() ));

    return R;
  }
private:
  template <typename wls_state_t, typename fom_t, typename scalar_t>
  void compute_impl(const wls_state_t                   & wlsState,
                    const fom_t                         & fomObj,
                    residual_t                          & wlsR, 
                    const scalar_t                       & dt) const
  {
    for (std::size_t iStep=0; iStep < fomStates_.k_; ++iStep){
      fomStates_.reconstructFomStateAt(wlsState, iStep);
    }
    //wlsR[0] = ...;
    for (std::size_t iStep=1; iStep < fomStates_.n_; ++iStep)
    {
      const auto & currentFomState = fomStates_[iStep];
      auto & currentFomResidual = wlsR[iStep];
      if (numAuxStates == 1){
        fomObj.template timeDiscreteResidual(iStep, time, dt,
                                             *currentFomResidual.data(),
                                             *currentFomState.data(),
                                             *fomStates_[iStep-1].data());
      }
    }
  }

protected:
  fom_states_container_type & fomStates_;
};
}}}



namespace pressio{ namespace rom{ namespace experimental{
// container for the fom states for WLS
template <std::size_t n, typename fom_state_type, typename reconstuctor_type>
class WlsFomStatesContainer{

  // put here usual things and overload operator [] so we can access the fom state at a given index
  // where the index typically is the step number
public:

  static constexpr std::size_t n_ = n;

  WlsFomStatesContainer() = delete;
  ~WlsFomStatesContainer() = default;

  template <
    typename _fom_state_type = fom_state_type,
    ::pressio::mpl::enable_if_t<
      ::pressio::containers::meta::is_vector_wrapper<_fom_state_type>::value
      > * = nullptr
    >
  WlsFomStatesContainer(const _fom_state_type & fomStateIn,
                        const reconstuctor_type & fomStateReconstr)
    : fomStates_{fomStateIn}, // does FOM states now hold all fomStateIn?
      fomStateReconstrObj_(fomStateReconstr)
  {
    this->resetContainersToZero();
  }

public:
  fom_state_type & operator[](std::size_t index) {
    return fomStates_[index];
  }

  fom_state_type const & operator[](std::size_t index) const{
    return fomStates_[index];
  }

  template <typename rom_state_t>
  void reconstructFomStateAt(const rom_state_t & romY, std::size_t index) const
  {
    fomStateReconstrObj_(romY[index], fomStates_[index]);
  }

private:
  // set all entries to zero for all members 
  void resetContainersToZero(){
    for (auto i=0; i<n; i++)
      ::pressio::containers::ops::set_zero(fomStates_[i]);
  }

private:
  mutable std::array<fom_state_type, n> fomStates_;
  const reconstuctor_type & fomStateReconstrObj_  = {};
};

}}}


*/


/*


namespace pressio{ namespace rom{ namespace experimental{
template <std::size_t n, typename fom_type, typename wls_state_type, typename ... Args>
struct DefaultWlsTypeGeneratorResidualApi{
  // native types of the full-order model (fom)
  using fom_t                   = fom_type;
  using scalar_t                = typename fom_t::scalar_type;
  using fom_native_state_t      = typename fom_t::state_type;
  using fom_native_residual_t   = typename fom_t::residual_type;
  // fom_state_t: this is a wrapper of the fom native state, and using a vector wrapper
  // seems like a good idea since typically the fom uses a vector for the state
  using fom_state_t             = ::pressio::containers::Vector<fom_native_state_t>;
  // fom_residual_t: this is a wrapper of the fom native residual, and using a vector wrapper
  // seems like a good idea since typically the fom uses a vector for the residual
  using fom_residual_t          = ::pressio::containers::Vector<fom_native_residual_t>;
  // wls_state_t is the type holding the wls state, so maybe this
  // is a multi-vector or some custom type that we make. On object of this
  // is supposed to own the generalized coordinates at the n steps of the target window.
  // the wls_state_type is passed by the user so it is defined in the main file.
  using wls_state_t             = wls_state_type;
  // wls_residual_t is the type holding the wls residual, so maybe this should be made
  // a multi-vector of fom_residual_t or some custom type that we make.
  // For this DefaultWlsTypeGenerator, we can make this a multi-vector.
  // If we need something else, we can create a different WlsTypeGenerator.
  using wls_residual_t          = ::pressio::containers::MultiVector<fom_residual_t>;
  // wls_matrix_t: an object of this type should hold somehow the wls matrix, so basically the large
  // J*phi that stems from the wls formulation. Since the WLS matrix has a block structure with dense blocks,
  // a basic version could be one where each block is a multivector (as it is done now for lspg) and
  // we use a std::list of multi-vectors to hold the full matrix.
//  using wls_matrix_t          = ;
  // decoder types are easy, see the LspgUnsteady classes
  using ic1 = ::pressio::mpl::variadic::find_if_unary_pred_t<
    ::pressio::rom::meta::is_legitimate_decoder_type, Args...>;
  using decoder_t = ::pressio::mpl::variadic::at_or_t<void, ic1::value, Args...>;
  static_assert(!std::is_void<decoder_t>::value and ic1::value < sizeof... (Args),
                "A valid decoder type must be passed to define a WLS problem");
  //using decoder_t             =  get it from args ;
  using decoder_jac_t           = typename decoder_t::jacobian_t;
  // fom state reconstructor type
  using fom_state_reconstr_t    = FomStateReconstructor<fom_state_t, decoder_t>;

  // fom_states_data: is used to store the fom states, use a custom class for this (see top of this file)
  using fom_states_container_t  = WlsFomStatesContainer<n, fom_state_t, fom_state_reconstr_t>;

  // here we are relying on the residual API, so we expect the user to tell us some details
  // on the stepper scheme used by the fom
  // --- find the total number of states needed ---
  using ic2 = ::pressio::mpl::variadic::find_if_unary_pred_t<
    ::pressio::ode::meta::impl::IsStepperTotalNumStatesSetter, Args...>;
  using tot_n_setter = ::pressio::mpl::variadic::at_or_t<void, ic2::value, Args...>;
  static_assert( !std::is_void<tot_n_setter>::value, "...");
  // numAuxStates is the number of auxiliary states needed by the fom object to
  // compute time discrete operators. This can be extracted from the args... like we do for lspg_unsteady
  static constexpr std::size_t auxStates = tot_n_setter::value - 1;
  // policy to compute the WLS residual: here we use the most basic one that computes things sequentially.
  // For other types of policies, we can create other type generators.
  using wls_residual_policy_t   = ::pressio::rom::experimental::WlsSequentialResidualPolicyForResidualApi<
    auxStates, wls_residual_t, fom_states_container_t>;
  // // policy to compute the WLS jacobian, do something similar for the residual one
  //using wls_jacobian_policy_t =  ;

};//end class
}}}


namespace pressio{ namespace rom{ namespace experimental{
template<
  template <std::size_t n, typename, typename, typename ...> class wls_type,
  //template <std::size_t n,  typename fom_type, typename wls_state_type, typename ... Args> class wls_type,
  std::size_t n,
  typename fom_type,
  typename wls_state_type,
  typename ... Args>
class WlsProblemGeneratorResidualApi
{

public:
  using wls_problem_t = wls_type<n,fom_type, wls_state_type, Args...>;

  using fom_t                   = typename wls_problem_t::fom_t;
  using scalar_t                = typename wls_problem_t::scalar_t;
  using fom_native_state_t      = typename wls_problem_t::fom_native_state_t;
  using fom_native_residual_t   = typename wls_problem_t::fom_native_residual_t;

  using fom_state_t             = typename wls_problem_t::fom_state_t;
  using fom_residual_t          = typename wls_problem_t::fom_residual_t;

  using wls_state_t             = typename wls_problem_t::wls_state_t;
  using wls_residual_t          = typename wls_problem_t::wls_residual_t;
  using wls_matrix_t            = typename wls_problem_t::wls_matrix_t;

  using decoder_t               = typename wls_problem_t::decoder_t;
  using fom_state_reconstr_t    = typename wls_problem_t::fom_state_reconstr_t;
  using fom_states_container_t  = typename wls_problem_t::fom_states_container_t;

  using wls_residual_policy_t   = typename wls_problem_t::wls_residual_policy_t;
  using wls_jacobian_policy_t   = typename wls_problem_t::wls_jacobian_policy_t;
  using wls_system_t            = typename wls_problem_t::wls_system_t;

  wls_state 			wls_state_type;
  wls_residual_policy_t         residualPolicy_;
  wls_jacobian_policy_t         jacobianPolicy_;
  wls_system_t                  systemObj_;
  fom_states_container_t        fomStates_;

  WlsProblemGeneratorResidualApi(const fom_t     & fomObj,
                                 const fom_native_state_t & fomNativeStateReference,
                                 decoder_t       & decoder,
                                 wls_state_t     & wlsInitialState,
                                 scalar_t       t0)
    :   fomStates_(fomNativeStateReference, decoder), 
	residualPolicy_(fomStates_){}  
//      t0_{t0},
//      dt0_{} {}
//      residualQuerier_{},
//      applyJacobQuerier_{},
//      fomStateReference_(fomStateReferenceNative),
//      fomStateReconstructor_(fomStateReference_, decoder),
//      fomStates_(fomStateReconstructor_, fomStateReference_),
//      residualPolicy_(fomStates_, residualQuerier_){}
};//end class
}}}


*/
int main(int argc, char *argv[]){
  using fom_t		= pressio::apps::Burgers1dEigen;
  //using fom_resid_t	= pressio::apps::Burgers1dEigenResidualApi::residual_type;
  //using fom_resid_t2	= fom_t::velocity_type;

  using scalar_t	= typename fom_t::scalar_type;

  using eig_dyn_vec	= Eigen::Matrix<scalar_t, -1, 1>;
  using rom_state_t	= pressio::containers::Vector<eig_dyn_vec>;

  using eig_dyn_mat	= Eigen::Matrix<scalar_t, -1, -1>;
  using decoder_jac_t	= pressio::containers::MultiVector<eig_dyn_mat>;
  using decoder_t	= pressio::rom::LinearDecoder<decoder_jac_t>;

  std::string checkStr {"PASSED"};

  //-------------------------------
  // app object
  int fomSize = 20;
  Eigen::Vector3d mu(5.0, 0.02, 0.02);
  fom_t appObj( mu, fomSize);
  appObj.setup();
  auto t0 = static_cast<scalar_t>(0);
  scalar_t dt = 0.01;

  // read from file the jacobian of the decoder
  int romSize = 11;
  // store modes computed before from file
  decoder_jac_t phi =
    pressio::rom::test::eigen::readBasis("basis.txt", romSize, fomSize);
  int numBasis = phi.numVectors();
  if( numBasis != romSize ) return 0;

  // create decoder Obj
  decoder_t decoderObj(phi);

  // for this problem, my reference state = initial state
  auto & yRef = appObj.getInitialState();




  constexpr std::size_t numStepsInWindow = 5;
  using stepper_stencil = ::pressio::ode::types::StepperTotalNumberOfStates<2>; 
  //using wls_problem      = pressio::rom::experimental:::WlsProblemGeneratorResidualApi<
//    DefaultWlsTypeGeneratorResidualApi, numStepsInWindow,
//    fom_t, wls_state_t, decoder_t, stepper_stencil, scalar_t>;



  // define ROM state
  rom_state_t yROM(romSize);
  using jtr_t = rom_state_t;
  using jtj_t = pressio::containers::MultiVector<eig_dyn_mat>;
  jtj_t jtj(romSize,romSize);
  jtr_t jtr(romSize);
  //using wls_t = pressio::rom::experimental::WlsSystem<fom_t,rom_state_t,JTJ_type,JTR_type,pressio::rom::wls::impl::JTJ_JTR_policy_smart>;
  using jtjr_policy_t = pressio::rom::wls::impl::JTJ_JTR_policy_smart<rom_state_t,fom_t,jtj_t,jtr_t>;
  using wls_api_t   = pressio::rom::wls::WlsSystemJTJApi<fom_t,rom_state_t,jtjr_policy_t>;
  using wls_system_t  = pressio::rom::wls::WlsSystem<wls_api_t, fom_t, rom_state_t,jtjr_policy_t>;
  // construct objects and test
  jtjr_policy_t jtjr_policy;
  wls_system_t wls(jtjr_policy);
  wls.JTJ_JTR(yROM,jtj,jtr); 


  using resid_t = rom_state_t;
  using jacobian_t = pressio::containers::MultiVector<eig_dyn_mat>;
  resid_t residual(romSize);
  jacobian_t jacobian(romSize,romSize);
  using residual_policy_t = pressio::rom::wls::impl::residual_policy_naive<rom_state_t,fom_t,resid_t>;
  using jacobian_policy_t = pressio::rom::wls::impl::jacobian_policy_naive<rom_state_t,fom_t,jacobian_t>;
  using wls_naive_api_t   = pressio::rom::wls::WlsSystemDefaultApi<fom_t,rom_state_t,residual_policy_t,jacobian_policy_t>;
  using wls_naive_system_t  = pressio::rom::wls::WlsSystem<wls_naive_api_t, fom_t, rom_state_t,residual_policy_t,jacobian_policy_t>;
  //construct objects and test
  residual_policy_t residual_policy;
  jacobian_policy_t jacobian_policy;
  wls_naive_system_t wls_naive(residual_policy,jacobian_policy);
  wls_naive.residual(yROM,residual); 


//  wls_system wls(residual_policy_t,jacobian_policy);

  /*
  //pressio::rom::experimental::WlsSystem<lspg_state_t,fom_t> wls(appobj); 
  // initialize to zero (this has to be done)
  yROM.putScalar(0.0);

  // define LSPG type
  //constexpr auto ode_case  = pressio::ode::ImplicitEnum::Euler;
  //using lspg_problem =  pressio::rom::LSPGUnsteadyProblem<pressio::rom::DefaultLSPGUnsteady, ode_case, fom_t, lspg_state_t, decoder_t>;
  using lspg_problem_t =  pressio::rom::LSPGUnsteadyProblem<pressio::rom::DefaultLSPGUnsteady,pressio::ode::ImplicitEnum::Euler, fom_t, lspg_state_t, decoder_t>;
  //using wls_problem_t = wlsProblemGenerator<fom_t, lspg_state_t, decoder_t>;
  using lspg_stepper_t = typename lspg_problem_t::lspg_stepper_t;
  //using wls_stepper_t = typename wls_problem_t::wls_stepper_t; 



  lspg_problem_t lspgProblem(appobj, yRef, decoderObj, yROM, t0);

  using eig_dyn_mat  = Eigen::Matrix<scalar_t, -1, -1>;
  using hessian_t  = pressio::containers::Matrix<eig_dyn_mat>;
  using solver_tag   = pressio::solvers::linear::iterative::LSCG;
  using linear_solver_t = pressio::solvers::iterative::EigenIterative<solver_tag, hessian_t>;
  linear_solver_t linSolverObj;

  // GaussNewton solver
  // hessian comes up in GN solver, it is (J phi)^T (J phi)
  // rom is solved using eigen, hessian is wrapper of eigen matrix
  using eig_dyn_mat  = Eigen::Matrix<scalar_t, -1, -1>;
  using gnsolver_t   =     pressio::solvers::iterative::GaussNewton<lspg_stepper_t, linear_solver_t>;

  gnsolver_t     solver   (lspgProblem.getStepperRef(), yROM, linSolverObj);

  //wlsProblem.residual(yROM);
  //wls_gnsolver_t solver2  (wlsProblem                  ,yROM, linSolverObj);
  solver.setTolerance(1e-13);
  // I know this should converge in few iters every step
  solver.setMaxIterations(2);

  // integrate in time
  pressio::ode::integrateNSteps(lspgProblem.getStepperRef(), yROM, 0.0, dt, 10, solver);

  // compute the fom corresponding to our rom final state
  auto yFomFinal = lspgProblem.getFomStateReconstructorCRef()(yROM);

  // this is a reproducing ROM test, so the final reconstructed state
  // has to match the FOM solution obtained with euler, same time-step, for 10 steps
  // const auto trueY = pressio::apps::test::Burg1DtrueImpEulerN20t010;
  const auto trueY = pressio::apps::test::Burgers1dImpGoldStatesBDF1::get(fomSize, dt, 0.10);
  for (auto i=0; i<yFomFinal.size(); i++){
    if (std::abs(yFomFinal[i] - trueY[i]) > 1e-10)
      checkStr = "FAILED";
  }

  std::cout << std::setprecision(14) << *yFomFinal.data() << std::endl;

  std::cout << checkStr <<  std::endl;
  */
  return 0;
}
