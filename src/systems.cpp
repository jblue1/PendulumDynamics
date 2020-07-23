#include "systems.hpp"
#include <stdio.h>

const double g = 9.81;

/**
 * Class constructor
 */
param_forced_pend::param_forced_pend(double *pend_params) {
  A = pend_params[0];
  L = pend_params[1];
  d = pend_params[2];
  omega = pend_params[3];
  b = pend_params[4];
  m = pend_params[5];
  k = pend_params[6];

  x.push_back(0.0);
  x.push_back(0.0);
}

// implement streaming observers

// constructor
param_forced_pend::streaming_observer_txt::streaming_observer_txt(
    std::ofstream &out)
    : write_out(out) {}

/**
 * Operator overload used by odeint to write out data to file.
 */
void param_forced_pend::streaming_observer_txt::operator()(const state_type &x,
                                                           double t) {
  write_out << t;
  for (size_t i = 0; i < x.size(); i++) {
    write_out << " " << x[i];
  }
  write_out << "\n";
}

// constructor
param_forced_pend::streaming_observer_arr::streaming_observer_arr(
    double *data, double dt, double trans_time, double simul_time)
    : count(0), data_buffer(data), trans_time(trans_time),
      num_points((int)simul_time + 1), points_per_sec((int)(1.0 / dt)) {}

/**
 * Operator overload used by odeint to write out data. Becuase HDF5 can only
 * write data stored contigiously in memory, solutions from odeint which
 * would more easily be stored in a N x 3 boost:array (the columns being time,
 * x1 and x2) is instead stored in a 3N x 1 array.
 */
void param_forced_pend::streaming_observer_arr::operator()(const state_type &x,
                                                           double t) {
  if (t < 1e-10) {
    // write initial conditions to array
    data_buffer[0] = t;
    data_buffer[num_points] = x[0];
    data_buffer[2 * (num_points)] = x[1];
  } else if (t >= trans_time) {
    if (count % points_per_sec == 0) {
      data_buffer[count / points_per_sec + 1] = t;
      for (size_t i = 0; i < x.size(); i++) {
        data_buffer[(i + 1) * (num_points) + count / points_per_sec + 1] = x[i];
      }
    }
    count++;
  }
}

/**
 * set state of pendulum to given theta and theta_dot
 */
void param_forced_pend::set_state(double theta, double theta_dot) {
  x[0] = theta;
  x[1] = theta_dot;
}

/**
 * Set pendulum parameters to given values
 */
void param_forced_pend::set_pend_params(double *params) {
  A = params[0];
  L = params[1];
  d = params[2];
  omega = params[3];
  b = params[4];
  m = params[5];
  k = params[6];
}

/**
 * Helper function for the parametrically forced pendulum with external
 * magnetic forcing. Calculates -0.5*r^2, where r is the distance between the
 * pendulum bob and the magnet
 */
double param_forced_pend::f(double t, double theta) {
  assert(A + L <= d); // make sure the pendulum bob doesn't go below the magnet
  double z = cos(omega * t);
  return -0.5 *
         (pow(L, 2) + pow(d, 2) + pow((A * z), 2) + 2 * A * L * z * cos(theta) -
          2 * A * d * z - 2 * d * L * cos(theta));
}

/**
 * Operator overload. Right hand side of equations of motion for the Pendulum
 * (dxdt = f(x, t)).
 */
void param_forced_pend::operator()(const state_type &x, state_type &dxdt,
                                   double t) {
  double z = cos(omega * t);
  dxdt[0] = x[1];
  dxdt[1] =
      -k * x[1] - (sin(x[0]) / L) * (g - A * z * pow(omega, 2) +
                                     (b / m) * (A * z - d) * exp(f(t, x[0])));
}

/**
 * Numerically solve eqns of motion for the pendulum until t_fin
 */
double param_forced_pend::solve(double dt, double abs_err, double rel_err,
                                double t_fin) {
  using namespace boost::numeric::odeint;
  typedef runge_kutta_fehlberg78<state_type> stepper_type;
  int num_steps = (int)(t_fin / dt);
  auto stepper = make_controlled(abs_err, rel_err, stepper_type());
  integrate_n_steps(stepper, boost::ref(*this), x, 0.0, dt, num_steps);
  return x[0];
}

/**
 * Numerically solve eqns of motion for pendulum until t_fin, write solution
 * to txt file.
 */
double param_forced_pend::solve(double dt, double abs_err, double rel_err,
                                double t_fin, std::ofstream &out) {
  using namespace boost::numeric::odeint;
  typedef runge_kutta_fehlberg78<state_type> stepper_type;
  int num_steps = (int)(t_fin / dt);
  auto stepper = make_controlled(abs_err, rel_err, stepper_type());

  integrate_n_steps(stepper, boost::ref(*this), x, 0.0, dt, num_steps,
                    streaming_observer_txt(out));
  return x[0];
}

/**
 * Numerically solve eqns of motion for pendulum until trans_time+simul_time.
 * Write solution at 1 second intervals startign after trans_time to an
 * array.
 */
double param_forced_pend::solve(double dt, double abs_err, double rel_err,
                                double trans_time, double simul_time,
                                double *data) {
  using namespace boost::numeric::odeint;
  typedef runge_kutta_fehlberg78<state_type> stepper_type;
  double t_fin = simul_time + trans_time - 1;
  double num_steps = (int)(t_fin / dt);
  auto stepper = make_controlled(abs_err, rel_err, stepper_type());

  integrate_n_steps(stepper, boost::ref(*this), x, 0.0, dt, num_steps,
                    streaming_observer_arr(data, dt, trans_time, simul_time));
  return x[0];
}
