
#include "pressio_rom.hpp"
#include "pressio_apps.hpp"
#include "utils_eigen.hpp"

template<typename fom_state_t, typename scalar_t>
struct MyTimeStepSetter
{
  const fom_state_t & fomState_;
  const scalar_t dt_;

  MyTimeStepSetter(const fom_state_t & fomState, scalar_t dt)
    : fomState_(fomState), dt_(dt)
  {}

  void operator()(const pressio::ode::types::step_t & step,
		  const scalar_t & time,
		  scalar_t & dt) const
  {
    std::cout << step << std::endl;
    dt = dt_;
    // fomState_ holds a reference to the current NATIVE fom state,
    // so you can use it here to do something with it
  }
};

int main(int argc, char *argv[]){
  using fom_t		= pressio::apps::Burgers1dEigen;
  using scalar_t	= typename fom_t::scalar_type;
  using native_state_t  = typename fom_t::state_type;
  using fom_state_t  = pressio::containers::Vector<native_state_t>;

  using eig_dyn_vec	= Eigen::Matrix<scalar_t, -1, 1>;
  using lspg_state_t	= pressio::containers::Vector<eig_dyn_vec>;

  using eig_dyn_mat	= Eigen::Matrix<scalar_t, -1, -1>;
  using decoder_jac_t	= pressio::containers::MultiVector<eig_dyn_mat>;
  using decoder_t	= pressio::rom::LinearDecoder<decoder_jac_t, fom_state_t>;

  std::string checkStr {"PASSED"};

  //-------------------------------
  // app object
  constexpr int numCell = 20;
  Eigen::Vector3d mu(5.0, 0.02, 0.02);
  fom_t appobj( mu, numCell);
  auto t0 = static_cast<scalar_t>(0);
  scalar_t dt = 0.01;

  // read from file the jacobian of the decoder
  constexpr int romSize = 11;
  // store modes computed before from file
  decoder_jac_t phi =
    pressio::rom::test::eigen::readBasis("basis.txt", romSize, numCell);
  const int numBasis = phi.numVectors();
  if( numBasis != romSize ) return 0;

  // create decoder obj
  decoder_t decoderObj(phi);

  // for this problem, my reference state = initial state
  auto & yRef = appobj.getInitialState();

  // define ROM state
  lspg_state_t yROM(romSize);
  // initialize to zero (this has to be done)
  pressio::ops::fill(yROM, 0.0);

  // define LSPG type
  using ode_tag  = pressio::ode::implicitmethods::Euler;
  using lspg_problem = typename pressio::rom::lspg::composeDefaultProblem<
    ode_tag, fom_t, lspg_state_t, decoder_t>::type;
  using lspg_stepper_t = typename lspg_problem::lspg_stepper_t;
  lspg_problem lspgProblem(appobj, yRef, decoderObj, yROM, t0);

  // linear solver
  using eig_dyn_mat  = Eigen::Matrix<scalar_t, -1, -1>;
  using hessian_t  = pressio::containers::Matrix<eig_dyn_mat>;
  using solver_tag   = pressio::solvers::linear::iterative::LSCG;
  using linear_solver_t = pressio::solvers::linear::Solver<solver_tag, hessian_t>;
  linear_solver_t linSolverObj;

  // GaussNewton solver with normal equations
  using nls_t = pressio::solvers::nonlinear::composeGaussNewton_t<
    lspg_stepper_t,
    pressio::solvers::nonlinear::DefaultUpdate,
    linear_solver_t>;
  nls_t solver(lspgProblem.getStepperRef(), yROM, linSolverObj);
  solver.setTolerance(1e-13);
  // I know this should converge in few iters every step
  solver.setMaxIterations(4);

  // my time step setter
  using setter_t = MyTimeStepSetter<native_state_t, scalar_t>;
  static_assert
    (pressio::ode::concepts::time_step_size_manager<
     setter_t, pressio::ode::types::step_t, scalar_t>::value,"");
  setter_t setter(lspgProblem.viewCurrentFomState(), dt);

  // integrate in time
  pressio::ode::advanceNSteps(lspgProblem.getStepperRef(), yROM, 0.0, 10, solver, setter);

  // compute the fom corresponding to our rom final state
  auto yFomFinal = lspgProblem.getFomStateReconstructorCRef()(yROM);

  // this is a reproducing ROM test, so the final reconstructed state
  // has to match the FOM solution obtained with euler, same time-step, for 10 steps
  // const auto trueY = pressio::apps::test::Burg1DtrueImpEulerN20t010;
  const auto trueY = pressio::apps::test::Burgers1dImpGoldStatesBDF1::get(numCell, dt, 0.10);
  for (auto i=0; i<yFomFinal.extent(0); i++){
    if (std::abs(yFomFinal[i] - trueY[i]) > 1e-10)
      checkStr = "FAILED";
  }

  std::cout << std::setprecision(14) << *yFomFinal.data() << std::endl;
  std::cout << checkStr <<  std::endl;
  return 0;
}
