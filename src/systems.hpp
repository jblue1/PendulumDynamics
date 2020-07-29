#ifndef SYSTEMSREF
#define SYSTEMSREF
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include <fstream>

/**
 * Class containing data and methods for numerically integrating a
 * parametrically driven pendulum with external forcing
 */
class param_forced_pend {
  typedef boost::array<double, 2>
      state_type; // used for integration with odeint

protected:
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

class pendulum_lyap : public param_forced_pend {
  static const size_t n = 3;        // dimension of phase space
  static const size_t num_lyap = 3; // number of displacement vectors

  // size of state vector, now including displacement vectors
  static const size_t N = n + n * num_lyap;

  typedef boost::array<double, N> state_type;
  typedef boost::array<double, num_lyap> lyap_type;

  // containers for lyapunov exponents and state of the system
  lyap_type lyap_exps;
  state_type x;


public:
  pendulum_lyap(double pend_params[7]); // constructor

  /**
   * Helper function for the parametrically forced pendulum with external
   * magnetic forcing. Calculates -0.5*r^2, where r is the distance between the
   * pendulum bob and the magnet
   * @param  A     Amplitude of pivot oscillations
   * @param  theta Angle
   * @param  phi   Phase of forcing (omega*t)
   * @return       -0.5r^2
   */
  double f(double theta, double phi);

  /**
   * Helper function for  calculating the jacobian
   * @param  A     Amplitude of pivot oscillations
   * @param  theta Angle
   * @param  phi   Phase of forcing (omega*t)
   * @return df2dtheta
   *
   * For the system of ODEs
   * dtheta/dt = f1(theta, theta_dot, phi)
   * dtheta_dot/dt = f2(theta, theta_dot, phi)
   * dphi/dt = f3(theta, theta_dot, phi)
   *
   * calculates df2/dtheta for use in calculating the peterbations
   */
  double df2dtheta(double theta, double phi);

  /**
   * See df2dtheta for explanation. Calculates df2dphi.
   * @param  A     Amplitude of pivot oscillations
   * @param  theta Angle
   * @param  phi   Phase of forcing (omega*t)
   * @rerurn df2dphi
   */
  double df2dphi(double theta, double phi);

  /**
   * Right hand side of equations of motion for the pendulum
   * (dxdt = f(x, t)).
   */
  void rhs(const state_type &x, state_type &dxdt, double t);

  /**
   * Contains equations of motion for pendulum and "tangent space" equations
   * for the peterbations.
   */
  void operator()(const state_type &x, state_type &dxdt, double t);

  /**
   * Set state of pendulum
   * @param  theta      Angle
   * @param  theta_dot  Angular velocity
   * @param  phase      Phase of forcing (omega*t)
   */
  void set_state(double theta, double theta_dot, double phase);

  void set_state(state_type state);

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
   * Initialize orthonormal peterbation vectors
   */
  void init_pert_vecs();

  state_type get_state();

  lyap_type get_exps();

  void print_theta();

  /**
   * Normalize vector to unit length
   * @param  first  Iterator pointing to start of vector
   * @param  last   Iterator pointing to end of vector
   */
  template <class Iterator, class T>
  void normalize(Iterator first, Iterator last, T norm) {
    for (Iterator ptr = first; ptr < last; ptr++) {
      *ptr /= norm;
    }
  }

  /**
   * For vectors u and v, subtracts the projection of v onto u from v and
   * stores it in place, ie
   * v = v - (<u, v>/<u, u>)*u
   * @param  first1  Iterator for start of v
   * @param  first2  Iterator for start of u
   */
  template <class Iterator> void subtr_proj(Iterator first1, Iterator first2) {
    double numer = std::inner_product(first1, first1 + n, first2, 0.0);
    double denom = std::inner_product(first2, first2 + n, first2, 0.0);
    for (Iterator ptr = first1; ptr < first1 + n; ptr++, first2++) {
      *ptr -= (numer / denom) * *first2;
    }
  }

  /**
   * Perform a modified gram schmidt on the peterbation vectors, store the sizes
   * of the vectors prior to normalizing.
   */
  void gram_schmidt();
};

#endif
