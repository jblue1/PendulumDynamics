#include <boost/numeric/odeint.hpp>
#include <stdio.h>
#include <string>
#include <fstream>
#include <cmath>

using namespace boost::numeric::odeint;

typedef std::vector<double> state_type;
const double g = 9.81;


class param_forced_pend {
protected:
  double A;
  double L;
  double d;
  double omega;
  double b;
  double m;
  double k;

  double f(double t, double theta, double A, double L, double d, double omega){
    double z = cos(omega * t);
    return -0.5*(pow(L, 2) + pow(d, 2) + pow((A*z), 2) + 2*A*L*z*cos(theta) -
            2*A*d*z - 2*d*L*cos(theta));
  }

public:
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

  void operator()(state_type &x, state_type &dxdt, double t){
    double z = cos(omega*t);
    dxdt[0] = x[1];
    dxdt[1] = -k*x[1] - (sin(x[0])/L)*(g - A*z*pow(omega, 2) +
              (b/m)*(A*z-d)*exp(f(t, x[0], A, L, d, omega)));

  }
};


struct streaming_observer {
  std::ofstream &write_out;

  streaming_observer(std::ofstream &out) : write_out(out) {}

  void operator()(const state_type &x, double t) const {
    write_out << t;
    for (size_t i=0; i< x.size(); i++) {
      write_out << "," << x[i];
    }
    write_out << "\n";
  }
};

static void usage(std::ostream &out, const char* msg){
    out << msg << "\n" << "\n";
    out << "  Usage:\n";
    out << "        solve_system file A\n";
    out << "  file  - name of output file\n";
    out << "  A     - amplitude of pivot oscillations\n";
    out << "  time  - Time to simulate system\n";
    out << "  1|2   - 1 for 8th order RK and 2 for 5th order\n";
    out << "  err   - value for absolute and relative error tolerance in the ODE solver\n\n";
    exit(1);
}


int main(int argc, char const *argv[]) {
  if (argc != 6) {
    usage(std::cerr, "Incorrect Number of parameters given.");
  }



  const char *file = argv[1];
  const double A = atof(argv[2]);
  const double t_fin = atof(argv[3]);
  const int solver = atoi(argv[4]);
  const double err = atof(argv[5]);
  const double L = g/pow(M_PI, 2);
  const double d = 4.0;
  const double omega = 2*M_PI;
  const double b = 50.0;
  const double m = 0.1;
  const double k = 0.2;
  const double abs_err = err;
  const double rel_err = err;

  const double dt = 0.01;
  std::ofstream write_out(file);
  assert(write_out.is_open());

  state_type x(2);
  x[1] = 1.0;
  x[0] = 1.0;


  if (solver == 1) {
    typedef runge_kutta_fehlberg78<state_type> error_stepper_type;
    error_stepper_type stepper;
    std::cout << "Using Fehlberg78\n";
    integrate_const(make_controlled<error_stepper_type>(abs_err, rel_err),
                    param_forced_pend(A, L, d, omega, b, m, k),x, 0.0, t_fin, dt,
                    streaming_observer(write_out));
    }

  else if (solver == 2){
    typedef runge_kutta_cash_karp54<state_type> error_stepper_type;
    error_stepper_type stepper;
    std::cout << "Using Cash-Karp54\n";
    integrate_const(make_controlled<error_stepper_type>(abs_err, rel_err),
                    param_forced_pend(A, L, d, omega, b, m, k),x, 0.0, t_fin, dt,
                    streaming_observer(write_out));
  }

  else {
    usage(std::cerr, "Please enter 1 or 2 to pick a solver\n");
  }

  write_out.close();
  return 0;
}
