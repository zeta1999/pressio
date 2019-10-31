
#include <gtest/gtest.h>
#include "CONTAINERS_ALL"
#include "ODE_ALL"
#include "SOLVERS_NONLINEAR"

struct MyApp
{

  /*
    dy
    -- = -10*y
    dt

    y(0) = 1;
    y(1) = 2;
    y(2) = 3;
   */
  using scalar_type = double;
  using state_type    = Eigen::VectorXd;
  using velocity_type = state_type;
  using residual_type = state_type;
  using jacobian_type = Eigen::SparseMatrix<double>;

public:
  void velocity(const state_type & yIn,
		const scalar_type & t,
		velocity_type & f) const{
    f = -10. * yIn;
  };

  velocity_type velocity(const state_type & yIn,
			 const scalar_type & t) const{
    velocity_type f(3);
    this->velocity(yIn, t, f);
    return f;
  };

  void jacobian(const state_type & yIn,
  		const scalar_type & t,
		jacobian_type & JJ) const
  {
    typedef Eigen::Triplet<scalar_type> Tr;
    std::vector<Tr> tripletList;
    tripletList.push_back( Tr( 0, 0, -10.) );
    tripletList.push_back( Tr( 1, 1, -10.) );
    tripletList.push_back( Tr( 2, 2, -10.) );
    JJ.setFromTriplets(tripletList.begin(), tripletList.end());
  };

  jacobian_type jacobian(const state_type & yIn,
  			 const scalar_type & t) const{
    jacobian_type JJ(3,3);
    this->jacobian(yIn, t, JJ);
    return JJ;
  };
  //--------------------------------------------

public:
  template <typename step_t, typename ... Args>
  void timeDiscreteResidual(const step_t & step,
                            const scalar_type & time,
                            const scalar_type & dt,
                            residual_type & R,
                            Args && ... states) const
  {
    this->timeDiscreteResidualImpl( step, time, dt, R, std::forward<Args>(states)... );
  }

  template <typename step_t, typename ... Args>
  void timeDiscreteJacobian(const step_t & step,
                            const scalar_type & time,
                            const scalar_type & dt,
                            jacobian_type & J,
                            Args && ... states) const
  {
    this->timeDiscreteJacobianImpl(step, time, dt, J, std::forward<Args>(states)... );
  }

  residual_type createTimeDiscreteResidualObject(const state_type & state) const
  {
    residual_type R(3);
    R.setConstant(0);
    return R;
  }

  jacobian_type createTimeDiscreteJacobianObject(const state_type & state) const
  {
    jacobian_type J(3,3);
    typedef Eigen::Triplet<scalar_type> Tr;
    std::vector<Tr> tripletList;
    tripletList.push_back( Tr( 0, 0, 0.) );
    tripletList.push_back( Tr( 1, 1, 0.) );
    tripletList.push_back( Tr( 2, 2, 0.) );
    J.setFromTriplets(tripletList.begin(), tripletList.end());
    return J;
  }

private:
  template <typename step_t, typename state_type>
  void timeDiscreteResidualImpl(const step_t & step,
				const scalar_type & time,
				const scalar_type & dt,
				residual_type & R,
				const state_type & yn,
				const state_type & ynm1) const
  {
    const auto f =  this->velocity(yn, time);
    R = yn - ynm1 - dt * f;
  }

  template <typename step_t, typename state_t>
  void timeDiscreteJacobianImpl(const step_t & step,
				const scalar_type & time,
				const scalar_type & dt,
				jacobian_type & J,
				const state_t & yn,
				const state_t & ynm1) const
  {
    J =  this->jacobian(yn, time);
    constexpr auto one = ::pressio::utils::constants::one<scalar_type>();
    J.coeffs() *= -dt;
    J.coeffRef(0,0) += one;
    J.coeffRef(1,1) += one;
    J.coeffRef(2,2) += one;
  }

};



struct Bdf1Solver
{
  using app_t		= MyApp;
  using sc_t		= typename app_t::scalar_type;
  using nstate_t	= typename app_t::state_type;
  using nveloc_t	= typename app_t::velocity_type;
  using njacobian_t	= typename app_t::jacobian_type;

  using state_t		= ::pressio::containers::Vector<nstate_t>;
  using res_t		= ::pressio::containers::Vector<nveloc_t>;
  using jac_t		= ::pressio::containers::Matrix<njacobian_t>;

  using stepper_t = ::pressio::ode::ImplicitStepper<::pressio::ode::ImplicitEnum::Euler,
						    state_t, res_t, jac_t, app_t>;

  using lin_solver_name = ::pressio::solvers::linear::iterative::Bicgstab;
  using lin_solver_t = ::pressio::solvers::iterative::EigenIterative<lin_solver_name, jac_t>;
  using nonlin_solver_t = ::pressio::solvers::NewtonRaphson<stepper_t, lin_solver_t, sc_t>;

  app_t appObj_ = {};
  state_t y_ = {};
  stepper_t stepperObj_;
  const sc_t dt_ = 0.01;

  Bdf1Solver(const state_t & yIn)
    : appObj_{}, y_{yIn}, stepperObj_{y_, appObj_}
  {}

  void run(int steps)
  {
    lin_solver_t linSolverObj;
    nonlin_solver_t solverO(stepperObj_, y_, linSolverObj);
    ::pressio::ode::integrateNSteps(stepperObj_, y_, 0.0, dt_, steps, solverO);
  };
};


struct CustomBdf1Solver
{
  using app_t		= MyApp;
  using sc_t		= typename app_t::scalar_type;
  using nstate_t	= typename app_t::state_type;
  using nresid_t	= typename app_t::residual_type;
  using njacobian_t	= typename app_t::jacobian_type;

  using state_t		= ::pressio::containers::Vector<nstate_t>;
  using res_t		= ::pressio::containers::Vector<nresid_t>;
  using jac_t		= ::pressio::containers::Matrix<njacobian_t>;

  using my_custom_order = ::pressio::ode::types::StepperOrder<1>;
  using my_num_states	= ::pressio::ode::types::StepperTotalNumberOfStates<2>;
  using stepper_t = ::pressio::ode::ImplicitStepper<::pressio::ode::ImplicitEnum::Arbitrary,
						    state_t, res_t, jac_t, app_t,
						    my_custom_order, my_num_states>;

  using traits = ::pressio::ode::details::traits<stepper_t>;
  using rp_t = typename traits::residual_policy_t;
  static_assert( !std::is_void<rp_t>::value, "");

  using std_r1 = ::pressio::ode::policy::ImplicitResidualStandardPolicyForArbitraryStepper<state_t, app_t, res_t>;
  using std_r2 = typename traits::standard_res_policy_t;
  static_assert( std::is_same<rp_t, std_r1>::value, "");

  using lin_solver_name = ::pressio::solvers::linear::iterative::Bicgstab;
  using lin_solver_t = ::pressio::solvers::iterative::EigenIterative<lin_solver_name, jac_t>;
  using nonlin_solver_t = ::pressio::solvers::NewtonRaphson<stepper_t, lin_solver_t, sc_t>;

  app_t appObj_		= {};
  state_t y_		= {};
  stepper_t stepperObj_;
  const sc_t dt_	= 0.01;

  CustomBdf1Solver(const state_t & yIn)
    : appObj_{}, y_{yIn}, stepperObj_{y_, appObj_}
  {}

  void run(int steps)
  {
    lin_solver_t linSolverObj;
    nonlin_solver_t solverO(stepperObj_, y_, linSolverObj);
    ::pressio::ode::integrateNSteps(stepperObj_, y_, 0.0, dt_, steps, solverO);
  };


  void runWithStepSizeManagerLambda(int steps)
  {
    lin_solver_t linSolverObj;
    nonlin_solver_t solverO(stepperObj_, y_, linSolverObj);
    using step_t = ::pressio::ode::types::step_t;
    const auto dtSetterLambda = [=](const step_t & step, const sc_t & time, sc_t & dt){
				  std::cout << " SETTING DT " << std::endl;
				  dt = dt_;
				};
    ::pressio::ode::integrateNSteps(stepperObj_, y_, 0.0, steps, solverO, dtSetterLambda);
  };

  void runWithStepSizeManagerLambdaWrong(int steps)
  {
    lin_solver_t linSolverObj;
    nonlin_solver_t solverO(stepperObj_, y_, linSolverObj);
    using step_t = ::pressio::ode::types::step_t;
    const auto dtSetterLambda = [=](const step_t & step, const sc_t & time, sc_t & dt){
				  std::cout << " SETTING DT " << std::endl;
				  dt = dt_*2.;
				};
    ::pressio::ode::integrateNSteps(stepperObj_, y_, 0.0, steps, solverO, dtSetterLambda);
  };

};



TEST(ode_implicit, arbitraryStepperRunEulerConstDt)
{
  ::pressio::containers::Vector<Eigen::VectorXd> y0(3);
  *y0.data() << 1,2,3;

  for (int N = 1; N < 10; N++){
    CustomBdf1Solver S1(y0);
    S1.run(N);
    std::cout << std::setprecision(14) << *S1.y_.data() << "\n";

    Bdf1Solver S2(y0);
    S2.run(N);
    std::cout << std::setprecision(14) << *S2.y_.data() << "\n";

    EXPECT_DOUBLE_EQ( S1.y_[0], S2.y_[0]);
    EXPECT_DOUBLE_EQ( S1.y_[1], S2.y_[1]);
    EXPECT_DOUBLE_EQ( S1.y_[2], S2.y_[2]);
  }
}


TEST(ode_implicit, arbitraryStepperRunEulerDtSetter)
{
  ::pressio::containers::Vector<Eigen::VectorXd> y0(3);
  *y0.data() << 1,2,3;

  for (int N = 1; N < 10; N++){
    CustomBdf1Solver S1(y0);
    S1.runWithStepSizeManagerLambda(N);
    std::cout << std::setprecision(14) << *S1.y_.data() << "\n";

    Bdf1Solver S2(y0);
    S2.run(N);
    std::cout << std::setprecision(14) << *S2.y_.data() << "\n";

    EXPECT_DOUBLE_EQ( S1.y_[0], S2.y_[0]);
    EXPECT_DOUBLE_EQ( S1.y_[1], S2.y_[1]);
    EXPECT_DOUBLE_EQ( S1.y_[2], S2.y_[2]);
  }
}


TEST(ode_implicit, arbitraryStepperRunEulerDtSetterWithWrongDt)
{
  ::pressio::containers::Vector<Eigen::VectorXd> y0(3);
  *y0.data() << 1,2,3;

  for (int N = 1; N < 10; N++){
    CustomBdf1Solver S1(y0);
    // use here a dt that we know wont work because it is wrong
    S1.runWithStepSizeManagerLambdaWrong(N);
    std::cout << std::setprecision(14) << *S1.y_.data() << "\n";

    Bdf1Solver S2(y0);
    S2.run(N);
    std::cout << std::setprecision(14) << *S2.y_.data() << "\n";

    ASSERT_TRUE( S1.y_[0] != S2.y_[0]);
    ASSERT_TRUE( S1.y_[1] != S2.y_[1]);
    ASSERT_TRUE( S1.y_[2] != S2.y_[2]);
  }
}