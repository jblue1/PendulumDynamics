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
  std::fill(x.begin(), x.end(), 0.0);
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
    write_out << " " << std::setprecision(14) << x[i];
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

pendulum_lyap::pendulum_lyap(double pend_params[7])
    : param_forced_pend(&pend_params[7]) {
  A = pend_params[0];
  L = pend_params[1];
  d = pend_params[2];
  omega = pend_params[3];
  b = pend_params[4];
  m = pend_params[5];
  k = pend_params[6];
  std::fill(x.begin(), x.end(), 0.0);
  std::fill(lyap_exps.begin(), lyap_exps.end(), 0.0);
}

/**
 * Helper function for the parametrically forced pendulum with external
 * magnetic forcing. Calculates -0.5*r^2, where r is the distance between the
 * pendulum bob and the magnet
 */
double pendulum_lyap::f(double theta, double phi) {
  double z = cos(phi);
  return -0.5 *
         (pow(L, 2) + pow(d, 2) + pow((A * z), 2) + 2 * A * L * z * cos(theta) -
          2 * A * d * z - 2 * d * L * cos(theta));
}

/**
 * Helper function for calculating the jacobian
 */
double pendulum_lyap::df2dtheta(double theta, double phi) {
  double z = cos(phi);
  return -(cos(theta) / L) * (g - A * pow(omega, 2) * z +
                              (b / m) * (A * z - d) * exp(f(theta, phi))) -
         (b / (2 * m * L)) * pow(sin(theta), 2) * (2 * A * L * z - 2 * d * L) *
             (A * z - d) * exp(f(theta, phi));
}

/**
 * Helper function for calculating the jacobian
 */
double pendulum_lyap::df2dphi(double theta, double phi) {
  return -(sin(theta) / L) *
         (A * pow(omega, 2) * sin(phi) +
          (b / m) * (-A * sin(phi) * exp(f(theta, phi)) +
                     (A * cos(phi) - d) * exp(f(theta, phi)) *
                         (pow(A, 2) * cos(phi) * sin(phi) +
                          A * L * sin(phi) * cos(theta) - A * d * sin(phi))));
}

void pendulum_lyap::rhs(const state_type &x, state_type &dxdt, double t) {
  double z = cos(omega * t);
  dxdt[0] = x[1];
  dxdt[1] =
      -k * x[1] - (sin(x[0]) / L) * (g - A * z * pow(omega, 2) +
                                     (b / m) * (A * z - d) * exp(f(t, x[0])));
  dxdt[2] = omega;
}

/**
 * Operator overload. Right hand side of equations of motion for the Pendulum
 * (dxdt = f(x, t)).
 */
void pendulum_lyap::operator()(const state_type &x, state_type &dxdt,
                               double t) {
  double z = cos(omega * t);
  dxdt[0] = x[1];
  dxdt[1] = -k * x[1] -
            (sin(x[0]) / L) * (g - A * z * pow(omega, 2) +
                               (b / m) * (A * z - d) * exp(f(x[0], t * omega)));
  dxdt[2] = omega;

  for (size_t l = 0; l < num_lyap; ++l) {
    const double *pert = x.begin() + 3 + l * 3;
    double *dpert = dxdt.begin() + 3 + l * 3;

    dpert[0] = pert[1];
    dpert[1] = df2dtheta(x[0], x[2]) * pert[0] - k * pert[1] +
               df2dphi(x[0], x[2]) * pert[2];
    dpert[2] = 0;
  }
}

void pendulum_lyap::set_state(double theta, double theta_dot, double phase) {
  x[0] = theta;
  x[1] = theta_dot;
  x[2] = phase;
}

void pendulum_lyap::set_state(state_type state) { x = state; }

/**
 * Initialize orthonormal peterbation vectors
 */
void pendulum_lyap::init_pert_vecs() {
  for (size_t i = 0; i < num_lyap; ++i) {
    x[n + n * i + i] = 1.0; // initialize orthagonal displacement vectors
  }
}

/**
 * Numerically solve eqns of motion for the pendulum until t_fin
 */
double pendulum_lyap::solve(double dt, double abs_err, double rel_err,
                            double t_fin) {
  using namespace boost::numeric::odeint;
  typedef runge_kutta_fehlberg78<state_type> stepper_type;
  int num_steps = (int)(t_fin / dt);
  auto stepper = make_controlled(abs_err, rel_err, stepper_type());
  integrate_n_steps(stepper, boost::ref(*this), x, 0.0, dt, num_steps);
  return x[0];
}

void pendulum_lyap::gram_schmidt() {
  state_type::iterator start = x.begin() + n;
  state_type::iterator vec1_beg = start, vec1_end = start + n;

  double norm[num_lyap];
  std::fill(norm, norm + num_lyap, 0.0);

  norm[0] = sqrt(std::inner_product(vec1_beg, vec1_end, vec1_beg, 0.0));
  normalize(vec1_beg, vec1_end, norm[0]);

  lyap_exps[0] += log(norm[0]);

  vec1_beg += n;
  vec1_end += n;

  for (int i = 1; i < num_lyap; i++, vec1_beg += n, vec1_end += n) {
    for (int j = 0; j < i; j++) {
      subtr_proj(vec1_beg, start + j * n);
    }
    norm[i] = sqrt(std::inner_product(vec1_beg, vec1_end, vec1_beg, 0.0));
    lyap_exps[i] += log(norm[i]);
    normalize(vec1_beg, vec1_end, norm[i]);
  }
}

boost::array<double, 12> pendulum_lyap::get_state() { return x; }

boost::array<double, 3> pendulum_lyap::get_exps() { return lyap_exps; }

void pendulum_lyap::print_theta() { std::cout << "phase: " << x[2] << "\n"; }
