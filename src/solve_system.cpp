#include <boost/numeric/odeint.hpp>
#include <stdio.h>
#include <string>
#include <fstream>
#include <cmath>

using namespace boost::numeric::odeint;

// type of container used by ODE solver to store state of the system
typedef std::vector<double> state_type;

const double g = 9.81; //gravitational acceleration



/**
 * Class containing equations and parameters associated with a parametrically
 * forced pendulum oscillating vertically about the origin with a repulsive
 * magnet located some distance below along the y axis.
 */
class param_forced_pend {
protected:
  double A; // Amplitude of pivot oscillations
  double L; // Length of pendulum
  double d; // Distance from origin to magnet
  double omega; // Angular frequency of oscillations
  double b; // Controls strength of magnet
  double m; // Mass of pendulum
  double k; // Linear friction coefficient


  /**
   * Helper function for the parametrically forced pendulum with external
   * magnetic forcing. Calculates -0.5*r^2, where r is the distance between the
   * pendulum bob and the magnet
   * @param  t     Time
   * @param  theta Angle
   * @param  A     Amplitude of pivot oscillations
   * @param  L     Length of Pendulum
   * @param  d     Distance from origin to magnetic
   * @param  omega Angular freq of pendulum oscillations
   * @return       -0.5r^2
   */
  double f(double t, double theta, double A, double L, double d, double omega){
    double z = cos(omega * t);
    return -0.5*(pow(L, 2) + pow(d, 2) + pow((A*z), 2) + 2*A*L*z*cos(theta) -
            2*A*d*z - 2*d*L*cos(theta));
  }

public:
  /**
   * Class constructor. Params described above.
   */
  param_forced_pend(double A, double L, double d, double omega, double b,
                    double m, double k){
    this->A=A;
    this->L=L;
    this->d=d;
    this->omega=omega;
    this->b=b;
    this->m=m;
    this->k=k;
  }

  /**
   * Operator overload. Right hand side of equations of motion for the Pendulum
   * (dxdt = f(x, t)).
   */
  void operator()(state_type &x, state_type &dxdt, double t){
    double z = cos(omega*t);
    dxdt[0] = x[1];
    dxdt[1] = -k*x[1] - (sin(x[0])/L)*(g - A*z*pow(omega, 2) +
              (b/m)*(A*z-d)*exp(f(t, x[0], A, L, d, omega)));

  }
};

/**
 * Structure used by ODE solver to write out solutions each step.
 */
struct streaming_observer {
  std::ofstream &write_out;


  streaming_observer(std::ofstream &out) : write_out(out) {} // constructor

  /**
   * Operator overload called by odeint to write solutions to file
   */
  void operator()(const state_type &x, double t) const {
    write_out << t;
    for (size_t i=0; i< x.size(); i++) {
      write_out << "," << x[i];
    }
    write_out << "\n";
  }
};
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
    out << "  file  - name of output file\n";
    out << "  A     - amplitude of pivot oscillations\n";
    out << "  time  - Time to simulate system\n";
    exit(1);
}


int main(int argc, char const *argv[]) {

  if (argc != 4) {
    usage(std::cerr, "Incorrect Number of parameters given.");
  }

  // store command line arguments
  const char *file = argv[1];
  const double A = atof(argv[2]);
  const double t_fin = atof(argv[3]);

  // define parameters for system
  const double L = g/pow(M_PI, 2);
  const double d = 4.0;
  const double omega = 2*M_PI;
  const double b = 50.0;
  const double m = 0.1;
  const double k = 0.2;

  // define parameters for ODE solver
  const double abs_err = 1e-10;
  const double rel_err = 1e-10;
  const double points_per_sec = 100.0;
  const double dt = 1/points_per_sec;
  const int num_points = lrint(t_fin * points_per_sec);


  // instantiate streaming_observer
  std::ofstream write_out(file);
  assert(write_out.is_open());

  state_type x(2);
  x[1] = 1.0;
  x[0] = 1.0;

  // create vector dictating the times at which we want solutions
  std::vector<double> times(num_points );
  for( size_t i=0 ; i<times.size() ; ++i )
    times[i] = dt*i;


  typedef runge_kutta_fehlberg78<state_type> error_stepper_type;
  error_stepper_type stepper;

  integrate_times(make_controlled( abs_err , rel_err , error_stepper_type() ),
                  param_forced_pend(A, L, d, omega, b, m, k),
                  x , times , dt , streaming_observer(write_out));

  write_out.close();
  return 0;
}
