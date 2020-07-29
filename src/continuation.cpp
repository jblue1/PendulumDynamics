#include "systems.hpp"
#include <cmath>
#include <filesystem>
#include <string>

using namespace boost::numeric::odeint;
namespace fs = std::filesystem;

const double g = 9.81; // gravitational acceleration

/**
 * Describe parameters taken and print errors if paramters are incorrectly
 * given.
 * @param  out   output stream for messages
 * @param  msg   Error message to print
 */
static void usage(std::ostream &out, const char *msg) {
  out << msg << "\n"
      << "\n";
  out << "  Usage:\n";
  out << "        solve_system dir time A theta_init theta_dot_init\n";
  out << "  dir            - path to directory to save file to\n";
  out << "  time           - Time to simulate system\n";
  out << "  A_init         - Initial driving Amplitude\n";
  out << "  A_step         - Amplitude stepsize\n";
  out << "  num_steps      - Number of amplitudes\n";
  out << "  theta_init     - initial angle\n";
  out << "  theta_dot_init - initial angular velocity\n";
  exit(1);
}

/**
 * Write parameters out to a file
 * @param  out      out stream for file to write to
 * @param  params   array with parameter values
 */
static void write_params(std::ostream &out, const double (&params)[15]) {
  out << "L = " << params[0] << "\n";
  out << "d = " << params[1] << "\n";
  out << "omega = " << params[2] << "\n";
  out << "b = " << params[3] << "\n";
  out << "m = " << params[4] << "\n";
  out << "k = " << params[5] << "\n";
  out << "theta_init = " << params[6] << "\n";
  out << "theta_dot_init = " << params[7] << "\n";
  out << "dt = " << params[8] << "\n";
  out << "abs_err = " << params[9] << "\n";
  out << "rel_err = " << params[10] << "\n";
  out << "t_fin = " << params[11] << "\n";
  out << "A_init = " << params[12] << "\n";
  out << "A_step = " << params[13] << "\n";
  out << "A_num_steps = " << params[14] << "\n";
}

int main(int argc, char const *argv[]) {

  if (argc != 8) {
    usage(std::cerr, "Incorrect Number of parameters given.");
  }

  // store command line arguments
  const std::string dir_name(argv[1]);
  const double t_fin = atoi(argv[2]);
  double A = atof(argv[3]);
  const double A_step = atof(argv[4]);
  const int num_steps = atof(argv[5]);
  const double theta_init = atof(argv[6]);
  const double theta_dot_init = atof(argv[7]);

  // define parameters for system
  const double L = g / pow(M_PI, 2);
  const double d = 4.0;
  const double omega = 2 * M_PI;
  const double b = 50.0;
  const double m = 0.1;
  const double k = 0.2;
  double pend_params[7] = {A, L, d, omega, b, m, k};

  // define parameters for ODE solver
  const double dt = 0.1;
  const double abs_err = 1e-10;
  const double rel_err = 1e-10;

  // create directory for saving data
  if (!fs::is_directory(dir_name) || !fs::exists(dir_name)) {
    fs::create_directory(dir_name);
  }

  else {
    usage(std::cerr,
          "Provide a path that does not lead to an existing directory or file");
  }

  const double params[15] = {
      L,  d,       omega,   b,     m, k,      theta_init,       theta_dot_init,
      dt, abs_err, rel_err, t_fin, A, A_step, (double)num_steps};

  // write parameters out to file
  std::ofstream write_out(dir_name + "/params.txt");
  write_params(write_out, params);
  write_out.close();

  // instantiate pendulum
  param_forced_pend pend(pend_params);

  // set initial state
  pend.set_state(theta_init, theta_dot_init);

  pend.solve(dt, abs_err, rel_err, 3000.0); // make sure on the attractor

  for (int i=0; i<num_steps; i++) {
    A += i*A_step;
    pend_params[0] = A;
    pend.set_pend_params(pend_params);

    std::string istr = std::to_string(i);
    std::ofstream write_data(dir_name + "/data" + istr + ".txt");
    pend.solve(dt, abs_err, rel_err, t_fin, write_data);
    write_data.close();

  }

  return 0;
}
