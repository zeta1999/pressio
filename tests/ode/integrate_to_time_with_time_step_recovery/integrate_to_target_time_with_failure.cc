
#include "pressio_ode.hpp"

template<typename ode_state_type>
struct MyFakeStepper
{
  template<typename solver_type>
  void doStep(ode_state_type & odeState,
	      const double & t,
	      const double & dt,
	      const pressio::ode::types::step_t & step,
	      solver_type & solver)
  {
    static_assert
      (::pressio::ode::concepts::legitimate_solver_for_implicit_stepper<
      solver_type, decltype(*this), ode_state_type>::value,
      "Invalid solver");

    if (step==3 and dt==0.1)
      throw pressio::eh::time_step_failure();
    if (step==5 and (dt==0.1 or dt==0.05))
      throw pressio::eh::time_step_failure();

    for (int i=0; i<odeState.extent(0); i++) odeState[i] += dt;
  }
};

struct MyFakeSolver
{
  template<typename system_t, typename state_t>
  void solve(const system_t & sys, state_t & state){}
};


int main(int argc, char *argv[])
{
  using scalar_t = double;
  using vec_t = Eigen::VectorXd;
  using ode_state_t = pressio::containers::Vector<vec_t>;

  ode_state_t y(3);
  y[0] = 1.0; y[1] = 2.0; y[2] = 3.0;

  MyFakeStepper<ode_state_t> stepper;
  MyFakeSolver solver;

  std::string checkStr= "PASSED";

  auto dtManager = [](const ::pressio::ode::types::step_t & step,
		      const double & time,
		      double & dt,
		      double & minDt,
		      double & dtRedFactor)
		{
		  dt = 0.1;
		  minDt = 0.01;
		  dtRedFactor=2.;
		};

  auto collector = [&checkStr](const ::pressio::ode::types::step_t & step,
		      const double & time,
		      const ode_state_t & y)
		   {
		     if (step==1){
		       if( std::abs(y[0]-1.1) > 1e-10 or
			   std::abs(y[1]-2.1) > 1e-10 or
			   std::abs(y[2]-3.1) > 1e-10)
			 checkStr = "FAILED";

		       if (std::abs(time-0.1) > 1e-10) checkStr="FAILED";
		     }
		     if (step==3){
		       if( std::abs(y[0]-1.25) > 1e-10 or
			   std::abs(y[1]-2.25) > 1e-10 or
			   std::abs(y[2]-3.25) > 1e-10)
			 checkStr = "FAILED";

		       if (std::abs(time-0.25) > 1e-10) checkStr="FAILED";
		     }
		     if (step==5){
		       if( std::abs(y[0]-1.375) > 1e-10 or
			   std::abs(y[1]-2.375) > 1e-10 or
			   std::abs(y[2]-3.375) > 1e-10)
			 checkStr = "FAILED";

		       if (std::abs(time-0.375) > 1e-10) checkStr="FAILED";
		     }
		     // std::cout << step << " "
		     // 	       << y[0] << " "
		     // 	       << y[1] << std::endl;
		   };

  pressio::ode::advanceToTargetTimeWithTimeStepRecovery
    (stepper, y, 0., 0.5, solver, dtManager, collector);

  if( std::abs(y[0]-1.575) > 1e-10 or
      std::abs(y[1]-2.575) > 1e-10 or
      std::abs(y[2]-3.575) > 1e-10)
    checkStr = "FAILED";

  std::cout << *y.data() << std::endl;
  std::cout << checkStr << std::endl;
  return 0;
}
