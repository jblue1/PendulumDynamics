#ifndef SYSTEMSREF
#define SYSTEMSREF
#include <boost/numeric/odeint.hpp>
#include <fstream>

typedef std::vector<double> state_type; // used for integration with odeint

/**
 * Class containing data and methods for numerically integrating a
 * parametrically driven pendulum with external forcing
 */
class param_forced_pend {
private:
  // paramaters for pendulum
  double A;     // Amplitude of pivot oscillations
  double L;     // Length of pendulum
  double d;     // Distance from origin to magnet
  double omega; // Angular frequency of oscillations
  double b;     // Controls strength of magnet
  double m;     // Mass of pendulum
  double k;     // Linear friction coefficient

  state_type x; // state of system x[0] = theta, x[1] = theta_dot

  /**
   * Struct used to write data from odeint to a csv file
   */
  struct streaming_observer_txt {
    std::ofstream &write_out;

    streaming_observer_txt(std::ofstream &out); // constructor

    /**
     * Operator overload to write out data
     */
    void operator()(const state_type &x, double t);
  };

  /**
   * Struct used to write data from odeint to an array
   */
  struct streaming_observer_arr {
    int count;           // times observer has been called
    double *data_buffer; // array to store data in
    double trans_time;   // array of params about the solver
    int num_points;      // number of data points to save
    int points_per_sec;  // points per second given by the solver

    streaming_observer_arr(double *data, double dt, double trans_time,
                           double simul_time); // constructor

    /**
     * Operator overload called by odeint to write solutions to array
     */
    void operator()(const state_type &x, double t);
  };

public:
  /**
   * Class constructor.
   * pend_params[0] - A
   * pend_params[1] - L
   * pend_params[2] - d
   * pend_params[3] - omega
   * pend_params[4] - b
   * pend_params[5] - m
   * pend_params[6] - k
   */
  param_forced_pend(double pend_params[7]);

  /**
   * Set the state of the pendulum
   * @param  theta      Angle (rad)
   * @param  theta_dot  Angular velocity (rad/s)
   */
  void set_state(double theta, double theta_dot);

  /**
   * Set pendulum parameters. Arguments same as in the class constructor.
   */
  void set_pend_params(double params[7]);

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
  double f(double t, double theta);

  /**
   * Operator overload. Right hand side of equations of motion for the pendulum
   * (dxdt = f(x, t)).
   */
  void operator()(const state_type &x, state_type &dxdt, double t);

  /**
   * Numerically solves equations of motion, storing the final state in x. Does
   * not call an observer.
   * @param   dt       Initial step size
   * @param   abs_err  See below
   * @param   rel_err  See below
   * @param   t_fin    Final time for integration
   * @return  Angle of the final state
   * NOTE: abs_err and rel_err are used in this equation
   *     abs_err + rel_err * ( |x| + * dt * |dxdt| )
   * If the error calculated by the solver for a given step is less than
   * the above expression, the step is accepted, otherwise a smaller step size
   * is chosen.
   */
  double solve(double dt, double abs_err, double rel_err, double t_fin);

  /**
   * Function overload. Solves equations of motion and writes data out to a txt
   * file.
   * @param   dt       Initial step size
   * @param   abs_err  See note above
   * @param   rel_err  See note above
   * @param   t_fin    Final time for integration
   * @param   out      ofstream for text file to write out to.
   * @return  Angle of the final state
   */
  double solve(double dt, double abs_err, double rel_err, double t_fin,
               std::ofstream &out);

  /**
   * Function overload. Solves equations of motion saves it in an array to later
   * be written to an h5 file.
   * @param   dt          Initial step size
   * @param   abs_err     See note above
   * @param   rel_err     See note above
   * @param   trans_time  Time to simulate before taking data to ensure system
   * is on the attractor.
   * @param   simul_time  Time after trans_time to simulate system
   * @param   data        Array to store data in
   * @return  Angle of the final state
   */
  double solve(double dt, double abs_err, double rel_err, double trans_time,
               double simul_time, double *data);
};

#endif
