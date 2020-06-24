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
    return -0.5*(pow(L, 2) + pow(d, 2) + pow((A*z), 2) + 2*A*L*z*cos(theta) - 2*A*d*z - 2*d*L*cos(theta));
  }

public:
  param_forced_pend(double A, double L, double d, double omega, double b, double m, double k){
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
    dxdt[1] = -k*x[1] - (sin(x[0])/L)*(g - A*z*pow(omega, 2) + (b/m)*(A*z-d)*exp(f(t, x[0], A, L, d, omega)));

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


int main(int argc, char const *argv[]) {

  typedef runge_kutta_cash_karp54< state_type >  error_stepper_type;
  error_stepper_type stepper;

  const double A = 0.02;
  const double L = g/pow(M_PI, 2);
  const double d = 4.0;
  const double omega = 2*M_PI;
  const double b = 50.0;
  const double m = 0.1;
  const double k = 0.2;
  const double abs_err = 1e-10;
  const double rel_err = 1e-10;


  const double dt = 0.01;

  std::ofstream write_out("test.csv");
  assert(write_out.is_open());

  double x0_step = 2*M_PI / 100;
  double x1_step = 6 / 5;
  for (size_t i = 0; i < 5; i++) {
    for (size_t j=0; j < 5; j++) {
      state_type x(2);
      x[0] = -M_PI + i*x0_step;
      x[1] = -3 + j*x1_step;
      integrate_const(make_controlled<error_stepper_type>(abs_err, rel_err), param_forced_pend(A, L, d, omega, b, m, k),
      x, 0.0, 500.0, dt);
    }
  }


  write_out.close();
  return 0;
}
