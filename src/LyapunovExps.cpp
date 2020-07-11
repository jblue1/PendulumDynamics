#include <iostream>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include "boost/numeric/odeint/gram_schmidt.hpp"
#include <cmath>

using namespace std;
using namespace boost::numeric::odeint;


// define parameters for system
const double g = 9.81;
const double L = g/pow(M_PI, 2);
const double d = 4.0;
const double omega = 2*M_PI;
const double b = 50.0;
const double m = 0.1;
const double k = 0.2;
const double A = 0.0;


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
double f(double theta, double phi)
  {
  double z = cos(phi);
  return -0.5*(pow(L, 2) + pow(d, 2) + pow((A*z), 2) + 2*A*L*z*cos(theta) -
                                      2*A*d*z - 2*d*L*cos(theta));
  }

/**
 * For system of ODEs
 *
 * dtheta/dt = f1(theta, theta_dot, phi)
 * dtheta_dot/dt = f2(theta, theta_dot, phi)
 * dphi/dt = f3(theta, theta_dot, phi)
 *
 * calculates df2/dtheta for use in calculating the peterbations
 */
double df2dtheta(double theta, double phi)
{
  double z = cos(phi);
  return -(cos(theta)/L)*(g-A*pow(omega, 2)*z + (b/m)*(A*z-d)*exp(f(theta, phi)))
  - (b/(2*m*L))*pow(sin(theta), 2)*(2*A*L*z - 2*d*L)*(A*z - d)*exp(f(theta, phi));
}


/**
 * See df2dtheta for explanation. Calculates df2dphi.
 */
double df2dphi(double theta, double phi)
{
  return
  -(sin(theta)/L) *
  (
    A*pow(omega, 2) * sin(phi) +
    (b/m) *
    (
      - A*sin(phi)*exp(f(theta, phi))
      + (A*cos(phi)-d) * exp(f(theta, phi)) *
      (
        pow(A, 2)*cos(phi)*sin(phi) + A*L*sin(phi)*cos(theta) - A*d*sin(phi)
      )
    )
  );
}

/**
 * Implementation of parametrically driven pendulum with external forcing.
 */
struct pendulum
{
  /**
   * Calulates RHS of dxdt for equations of motion of the pendulum. Template deals
   * with peturbed system using different state types.
   * @param  x      State vector of the system
   * @param  dxdt   Vector of derivatives
   * @param  t      Time
   */
  template<class State, class Deriv>
  void operator()(const State &x_, Deriv &dxdt_, double t) const
  {
    typename boost::range_iterator< const State >::type x = boost::begin( x_ );
    typename boost::range_iterator< Deriv >::type dxdt = boost::begin( dxdt_ );

    double z = cos(x[2]);
    dxdt[0] = x[1];
    dxdt[1] = -k*x[1] - (sin(x[0])/L)*(g - A*z*pow(omega, 2) +
              (b/m)*(A*z-d)*exp(f(x[0], x[2])));
    dxdt[2] = omega;
  }
};

//system with peterbations
const size_t n=3; // dimension of phase space
const size_t num_of_lyap = 3; // number of displacement vectors

// size of state vector, now including displacement vectors
const size_t N = n + n*num_of_lyap;

typedef boost::array< double , N > state_type;
typedef boost::array< double , num_of_lyap > lyap_type;

/**
 * Update whole system including peturbations.
 * @param  x    State vector
 * @param  dxdt Derivatives vector
 * @param  t    Time
 */
void pendulum_with_lyap( const state_type &x , state_type &dxdt , double t )
{
    pendulum()( x , dxdt , t );

    for( size_t l=0 ; l<num_of_lyap ; ++l )
    {
        const double *pert = x.begin() + 3 + l * 3;
        double *dpert = dxdt.begin() + 3 + l * 3;

        dpert[0] = pert[1];
        dpert[1] = df2dtheta(x[0], x[2])*pert[0] - k*pert[1] +
                   df2dphi(x[0], x[2])*pert[2];
        dpert[2] = 0;
    }
}

int main(int argc, char const *argv[]) {
  state_type x; // initialize state vector
  lyap_type lyap; // initialize exponents

  fill( x.begin() , x.end() , 0.0 );
  x[0] = -.5;
  x[1] = -1.5;
  x[2] = 0.0;

  const double dt = 1e-3;

  typedef runge_kutta_fehlberg78<state_type,
                                 double,
                                 state_type,
                                 double,
                                 range_algebra> error_stepper_type;
  error_stepper_type error_stepper;
  cout << "Taking transient steps \n";

  // transient steps to make sure system is on the attractor
  integrate_n_steps( error_stepper, pendulum() , std::make_pair( x.begin() , x.begin() + n ) , 0.0 , dt , 100 );
  cout << "Finished transient steps \n";

  fill( x.begin()+n , x.end() , 0.0 );

  for( size_t i=0 ; i<num_of_lyap ; ++i )
  {
    x[n+n*i+i] = 1.0; // initialize orthagonal displacement vectors
  }

  fill( lyap.begin() , lyap.end() , 0.0 );

  double t = 0.0;
  size_t count = 0;

  cout << "Begining exponent calculation \n";
  while( true )
  {

      t = integrate_n_steps(make_controlled<error_stepper_type>(1e-12, 1e-12),
      pendulum_with_lyap , x , t , dt , 100 ); // take 100 integration steps

      gram_schmidt< num_of_lyap >( x , lyap , n );
	  if (count % 10000 == 0)
      {
		double sum = 0;
		cout << "time: " << t;
		for( size_t i=0 ; i<num_of_lyap ; ++i )
	    {		  
			cout << "\t" << lyap[i] / t ;
			sum += lyap[i];
        }
		cout << "\t" << "sum:  " << sum;
		cout << endl;
      }

      ++count;
   }
  return 0;
}
