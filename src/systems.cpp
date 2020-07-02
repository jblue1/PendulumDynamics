#include "systems.hpp"
#include <fstream>
#include <stdio.h>
#include <boost/numeric/odeint.hpp>

const double g = 9.81;

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
double param_forced_pend::f(double t, double theta, double A, double L,
                            double d, double omega)
  {

  double z = cos(omega * t);
  return -0.5*(pow(L, 2) + pow(d, 2) + pow((A*z), 2) + 2*A*L*z*cos(theta) -
                                      2*A*d*z - 2*d*L*cos(theta));
}
/**
 * Class constructor
 */
param_forced_pend::param_forced_pend(double A, double L, double d, double omega,
                                     double b, double m, double k)
  {
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
  void param_forced_pend::operator()(state_type &x, state_type &dxdt, double t)
  {
    double z = cos(omega*t);
    dxdt[0] = x[1];
    dxdt[1] = -k*x[1] - (sin(x[0])/L)*(g - A*z*pow(omega, 2) +
              (b/m)*(A*z-d)*exp(f(t, x[0], A, L, d, omega)));
  }
