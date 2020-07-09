#include <boost/numeric/odeint.hpp>
#include <string>
#include <fstream>
#include <cmath>
#include <filesystem>
#include "systems.hpp"
#include "streaming_observers.hpp"


using namespace boost::numeric::odeint;
namespace fs = std::filesystem;

// type of container used by ODE solver to store state of the system
typedef std::vector<double> state_type;

const double g = 9.81; //gravitational acceleration


/**
 * Describe parameters taken and print errors if paramters are incorrectly
 * given.
 * @param out output stream for messages
 * @param msg Error message to print
 */
static void usage(std::ostream &out, const char* msg){
    out << msg << "\n" << "\n";
    out << "  Usage:\n";
    out << "        solve_system file A\n";
    out << "  dir            - path to directory to save file to\n";
    out << "  time           - Time to simulate system\n";
    out << "  A              - Driving Amplitude\n";
    out << "  theta_init     - initial angle\n";
    out << "  theta_dot_init - initial angular velocity\n";
    exit(1);
}

static void write_params(std::ostream &out, const double (&params)[6]) {
  out << "L = " << params[0] << "\n";
  out << "d = " << params[1] << "\n";
  out << "omega = " << params[2] << "\n";
  out << "b = " << params[3] << "\n";
  out << "m = " << params[4] << "\n";
  out << "k = " << params[5] << "\n";
}



int main(int argc, char const *argv[]) {

  if (argc != 6) {
    usage(std::cerr, "Incorrect Number of parameters given.");
  }

  // store command line arguments
  const std::string dir_name(argv[1]);
  const int t_fin = atoi(argv[2]);
  const double A = atof(argv[3]);
  const double theta_init = atof(argv[4]);
  const double theta_dot_init = atof(argv[5]);

  // create directory to save data to
  if (!fs::is_directory(dir_name) || !fs::exists(dir_name)) {
    fs::create_directory(dir_name);
  }

  else {
    usage(std::cerr, "Provide a path that does not lead to an existing directory or file");
  }

  // define parameters for system
  const double L = g/pow(M_PI, 2);
  const double d = 4.0;
  const double omega = 2*M_PI;
  const double b = 50.0;
  const double m = 0.1;
  const double k = 0.2;
  const double params[6] = {L, d, omega, b, m, k};

  //write parameters out to file
  std::ofstream write_out(dir_name + "/params.txt");
  write_params(write_out, params);

  std::ofstream write_data(dir_name + "/data.txt");

  // define parameters for ODE solver
  const double abs_err = 1e-10;
  const double rel_err = 1e-10;
  const int points_per_sec = 10;
  const double dt = 1.0/points_per_sec;
  const int num_points = points_per_sec * t_fin + 1;

  std::cout << "num points: " << num_points << "\n";

  // create vector dictating the times at which we want solutions
  std::vector<double> times(num_points);
  times[0] = 0.0;
  for( size_t i=1 ; i<times.size() ; ++i ){
    times[i] = dt*i + 1000.0;
    std::cout << times[i] << "\n";
  }

  typedef runge_kutta_fehlberg78<state_type> error_stepper_type;
  error_stepper_type stepper;

  //instantiate state vector
  state_type x(2);
  x[0] = theta_init;
  x[1] = theta_dot_init;


  integrate_times(make_controlled( abs_err , rel_err , error_stepper_type() ),
                  param_forced_pend(A, L, d, omega, b, m, k),
                  x , times, dt , streaming_observer_csv(write_data));
  write_data.close();


  return 0;
}
