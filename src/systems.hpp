#ifndef SYSTEMSREF
#define SYSTEMSREF
#include <boost/numeric/odeint.hpp>
#include <fstream>

typedef std::vector<double> state_type;

class param_forced_pend
{
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
  double f(double t, double theta, double A, double L, double d, double omega);

public:
  /**
   * Class constructor. Params described above.
   */
  param_forced_pend(double A, double L, double d, double omega, double b,
                    double m, double k);


  /**
  * Operator overload. Right hand side of equations of motion for the pendulum
  * (dxdt = f(x, t)).
  */
  void operator()(state_type &x, state_type &dxdt, double t);
};



#endif
