#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

// TODO: Set the timestep length and duration
size_t N = 10;
double dt = 0.1;

// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.
const double Lf = 2.67;

// set reference velocity
double ref_v = 30;

// set start position for every element of variables
size_t x_start = 0;
size_t y_start = x_start + N;
size_t psi_start = y_start + N;
size_t v_start = psi_start + N;
size_t cte_start = v_start + N;
size_t epsi_start = cte_start + N;
size_t delta_start = epsi_start + N;
size_t a_start = delta_start + N - 1;

class FG_eval {
 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  FG_eval(Eigen::VectorXd coeffs) { this->coeffs = coeffs; }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {
    // TODO: implement MPC
    // `fg` a vector of the cost constraints, `vars` is a vector of variable values (state & actuators)
    // NOTE: You'll probably go back and forth between this function and
    // the Solver function below.
	
	// note: "vars" gives every thing we need and we should change the values in "fg" ( to call optimizer to minimize these values ) 
	//============================ part1. cost value: fg[0] ==========================//
	  AD<double> loss = 0.0;
	  
	  // 1.1: cte and epsi, for every point there is always cte and epsi, so there are N iterations
	  // in order to minimize cte and epsi
	  for (int i = 0; i < int(N); i++)
	  {
		  loss += 2000 * CppAD::pow(vars[cte_start + i], 2) + 1000 * CppAD::pow(vars[epsi_start + i], 2);
	  }

	  // 1.2: velocity, for every point there is always velocity, so there are N iterations
	  // in order to minimize the difference between car's velocity and reference velocity
	  for (int i = 0; i < int(N); i++)
	  {
		  loss += CppAD::pow((vars[v_start + i] - ref_v), 2);
	  }

	  // 1.3: control input, between two points we need once time control input, so there are N-1 iterations
	  // in order to make the control input always in a samll magnitude
	  for (int i = 0; i < int(N-1); i++)
	  {
		  loss += 5 * CppAD::pow((vars[a_start + i]), 2) + 5 * CppAD::pow((vars[delta_start + i]), 2);
	  }

	  // 1.4: the change of control input, for every twice control input we can calculate differece for one time, so there are N-2 iterations
	  // in order to make the control input change steadily.
	  for (int i = 0; i < int(N - 2); i++)
	  {
		  loss += 200 * CppAD::pow((vars[a_start + i + 1] - vars[a_start + i]), 2) + 50 * CppAD::pow((vars[delta_start + i + 1] - vars[delta_start + i]), 2);
	  }

	  fg[0] = loss;
	  //============================ part1. cost value: fg[0] ==========================//


	  //============================ part2. constrain of motion model: fg[1::] ==========================//
	  fg[1 + x_start] = vars[x_start];
	  fg[1 + y_start] = vars[y_start];
	  fg[1 + psi_start] = vars[psi_start];
	  fg[1 + v_start] = vars[v_start];
	  fg[1 + cte_start] = vars[cte_start];
	  fg[1 + epsi_start] = vars[epsi_start];

	  for (int t = 1; t < int(N); t++) {
		  AD<double> x1 = vars[x_start + t];
		  AD<double> x0 = vars[x_start + t - 1];

		  AD<double> y1 = vars[y_start + t];
		  AD<double> y0 = vars[y_start + t - 1];

		  AD<double> psi1 = vars[psi_start + t];
		  AD<double> psi0 = vars[psi_start + t - 1];

		  AD<double> delta0 = vars[delta_start + t - 1];

		  AD<double> v1 = vars[v_start + t];
		  AD<double> v0 = vars[v_start + t - 1];

		  AD<double> a0 = vars[a_start + t - 1];

		  AD<double> epsi0 = vars[epsi_start + t - 1];
		  AD<double> epsi1 = vars[epsi_start + t];

		  AD<double> cte0 = vars[cte_start + t - 1];
		  AD<double> cte1 = vars[cte_start + t];

		  AD<double> f0 = coeffs[0] + coeffs[1] * x0 + coeffs[2] * x0 * x0 + coeffs[3] * x0 * x0 * x0;
		  AD<double> psides0 = CppAD::atan((3 * coeffs[3] * x0 * x0 + 2 * coeffs[2] * x0 + coeffs[1]));

		  fg[1 + x_start + t] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
		  fg[1 + y_start + t] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
		  fg[1 + psi_start + t] = psi1 - (psi0 - v0 / Lf * delta0 * dt);
		  // Note: In the simulator, a positive value implies a right turn
		  // so we change the equation here
		  fg[1 + v_start + t] = v1 - (v0 + a0 * dt);
		  fg[1 + cte_start + t] =
			  cte1 - ((f0 - y0) + (v0 * CppAD::sin(epsi0) * dt));
		  fg[1 + epsi_start + t] =
			  epsi1 - ((psi0 - psides0) + v0 * delta0 / Lf * dt);
	  }
	  //============================ part2. constrain of motion model: fg[1::] ==========================//
  
  }
};

//
// MPC class definition implementation.
//
MPC::MPC() {}
MPC::~MPC() {}

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
  bool ok = true;
  typedef CPPAD_TESTVECTOR(double) Dvector;

  // TODO: Set the number of model variables (includes both states and inputs).
  // For example: If the state is a 4 element vector, the actuators is a 2
  // element vector and there are 10 timesteps. The number of variables is:
  //
  // 4 * 10 + 2 * 9
  size_t n_vars = N * 6 + (N - 1) * 2; //state is a 6 element vector : x, y, psi, v, cte, epsi; actuators is a 2 element vector : a, delta;
  // TODO: Set the number of constraints
  size_t n_constraints = N * 6; // for every element in state, there is always a constraint.

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);
  for (int i = 0; i < int(n_vars); i++) {
    vars[i] = 0;
  }
  double x = state[0];
  double y = state[1];
  double psi = state[2];
  double v = state[3];
  double cte = state[4];
  double epsi = state[5];

  vars[x_start] = x;
  vars[y_start] = y;
  vars[psi_start] = psi;
  vars[v_start] = v;
  vars[cte_start] = cte;
  vars[epsi_start] = epsi;

  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);
  // TODO: Set lower and upper limits for variables.

  // in order to initailize vars_lowerbound and vars_upperbound, set a very high upper limit and very low lower limit.
  for (int i = 0; i < int(delta_start); i++) {
	  vars_lowerbound[i] = -1.0e19;
	  vars_upperbound[i] = 1.0e19;
  }

  // delta always in [ -25 grad, 25 grad]
  for (int i = delta_start; i < int(a_start); i++)
  {
	  vars_lowerbound[i] = -0.436332 * Lf;
	  vars_upperbound[i] = 0.436332 * Lf;
  }
  // a always in [-1.0,1.0]
  for (int i = a_start; i < int(n_vars); i++)
  {
	  vars_lowerbound[i] = -1.0;
	  vars_upperbound[i] = 1.0;
  }

  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (int i = 0; i < int(n_constraints); i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }
  constraints_lowerbound[x_start] = x;
  constraints_lowerbound[y_start] = y;
  constraints_lowerbound[psi_start] = psi;
  constraints_lowerbound[v_start] = v;
  constraints_lowerbound[cte_start] = cte;
  constraints_lowerbound[epsi_start] = epsi;

  constraints_upperbound[x_start] = x;
  constraints_upperbound[y_start] = y;
  constraints_upperbound[psi_start] = psi;
  constraints_upperbound[v_start] = v;
  constraints_upperbound[cte_start] = cte;
  constraints_upperbound[epsi_start] = epsi;
  // object that computes objective and constraints
  FG_eval fg_eval(coeffs);

  //
  // NOTE: You don't have to worry about these options
  //
  // options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          0.5\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  // Cost
  auto cost = solution.obj_value;
  std::cout << "Cost " << cost << std::endl;

  // TODO: Return the first actuator values. The variables can be accessed with
  // `solution.x[i]`.
  //
  // {...} is shorthand for creating a vector, so auto x1 = {1.0,2.0}
  // creates a 2 element double vector.
  vector<double> result;
  result.push_back(solution.x[delta_start]);
  result.push_back(solution.x[a_start]);
  cout << "delta:" <<solution.x[delta_start] << "  a:" << solution.x[a_start] <<endl ;
  for (int i = 0; i < int(N - 1); i++)
  {
	  result.push_back(solution.x[x_start + i + 1]);
	  result.push_back(solution.x[y_start + i + 1]);
  }

  return result;
}
