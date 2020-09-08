#include "systems.hpp"
//#include <boost/numeric/odeint.hpp>
#include <cmath>
#include <iostream>

using namespace std;
using namespace boost::numeric::odeint;

/**
 * Describe parameters taken and print errors if paramters are incorrectly
 * given.
 * @param out output stream for messages
 * @param msg Error message to print
 */
static void usage(std::ostream &out, const char *msg) {
  out << msg << "\n"
      << "\n";
  out << "  Usage:\n";
  out << "        solve_system  A  theta_init  theta_dot_init\n";
  out << "  A              - Driving Amplitude\n";
  out << "  theta_init     - initial angle\n";
  out << "  theta_dot_init - initial angular velocity\n";
  exit(1);
}

int main(int argc, char const *argv[]) {
  if (argc != 4) {
    usage(std::cerr, "Incorrect Number of parameters given.");
  }

  // define parameters for system
  const double g = 9.81;
  const double L = g / pow(M_PI, 2);
  const double d = 4.0;
  const double omega = 2 * M_PI;
  const double b = 50.0;
  const double m = 0.1;
  const double k = 0.2;

  // store command line arguments
  const double A = atof(argv[1]);
  const double theta_init = atof(argv[2]);
  const double theta_dot_init = atof(argv[3]);

  double pend_params[7] = {A, L, d, omega, b, m, k};

  // define parameters for ODE solver
  const double dt = 0.1;
  const double abs_err = 1e-10;
  const double rel_err = 1e-10;
  double phase = 0.0;

  pendulum_lyap pend(pend_params);

  pend.set_state(theta_init, theta_dot_init, phase);

  // solve for 5000 seconds to make sure system is on the attractor
  pend.solve(dt, abs_err, rel_err, 5000.0);

  pend.init_pert_vecs();

  int t = 0.0;
  for (int i = 0; i < 1000; i++) {
    pend.solve(dt, abs_err, rel_err, 10.0);
    t += 10;

    pend.gram_schmidt();

    boost::array<double, 3> exps = pend.get_exps();

    if (t % 1000 == 0) {
      std::cout << "t: " << t << "\n";
      std::cout << "exps: " << exps[0] / t << " " << exps[1] / t << " "
                << exps[2] / t << "\n"
                << "\n";
    }
  }

  return 0;
}
