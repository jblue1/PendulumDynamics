#ifndef SYSTEMSREF
#define SYSTEMSREF
#include <boost/numeric/odeint.hpp>
#include <fstream>

typedef std::vector<double> state_type;

class param_forced_pend
{
private:
  // paramaters for pendulum
  double A; // Amplitude of pivot oscillations
  double L; // Length of pendulum
  double d; // Distance from origin to magnet
  double omega; // Angular frequency of oscillations
  double b; // Controls strength of magnet
  double m; // Mass of pendulum
  double k; // Linear friction coefficient

  state_type x; // state of system x[0] = theta, x[1] = theta_dot


  /**
  * Struct used to write data from odeint to a csv file
  */
  struct streaming_observer_txt
  {
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
   struct streaming_observer_arr
   {
     int count;
     double *data_buffer;
     double *solve_params;
     int num_points;
     int points_per_sec;
     int points_saved;



     streaming_observer_arr(double *data, double *solve_params); // constructor

     /**
      * Operator overload called by odeint to write solutions to array
      */
     void operator()(const state_type &x, double t);
    };



public:


  /**
  * Class constructor. Params described above.
  */
  param_forced_pend(double pend_params[7]);



  void set_state(double theta, double theta_dot);

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

  double solve(double solve_params[5]);

  double solve(double solve_params[5], std::ofstream &out);

  double solve(double solve_params[5], double *data);

};



#endif
