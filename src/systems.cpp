#include "systems.hpp"
#include <stdio.h>

const double g = 9.81;

/**
 * Class constructor
 */
 param_forced_pend::param_forced_pend(double *pend_params)
  {
    A = pend_params[0];
    L = pend_params[1];
    d = pend_params[2];
    omega = pend_params[3];
    b = pend_params[4];
    m = pend_params[5];
    k = pend_params[6];

    x.push_back(0.0);
    x.push_back(0.0);
  }

  //implement streaming observers

  // constructor
  param_forced_pend::streaming_observer_txt::streaming_observer_txt(std::ofstream &out):
    write_out(out){}

  /**
  * Operator overload used by odeint to write out data to file.
  */
  void param_forced_pend::streaming_observer_txt::operator()(const state_type &x, double t)
  {
    write_out << t;
    for (size_t i=0; i < x.size(); i++) {
      write_out << " " << x[i];
    }
    write_out << "\n";
  }

  // constructor
  param_forced_pend::streaming_observer_arr::streaming_observer_arr(double *data, double *params) :
    count(0),
    data_buffer(data),
    solve_params(params),
    num_points((int) solve_params[4] + 1),
    points_per_sec((int) (1.0/params[0])),
    points_saved(0)
    {}

  /**
  * Operator overload used by odeint to write out data. Becuase HDF5 can only
  * write data stored contigiously in memory, solutions from odeint which
  * would more easily be stored in a N x 3 boost:array (the columns being time,
  * x1 and x2) is instead stored in a 3N x 1 array.
  */
  void param_forced_pend::streaming_observer_arr::operator()(const state_type &x, double t)
  {
    if (t < 1e-10)
    {
      data_buffer[0] = t;
      data_buffer[num_points] = x[0];
      data_buffer[2*(num_points)] = x[1];
      points_saved++;
    }
    else if (t >= solve_params[3])
    {
      //std::cout << "time: " << t << "\n";
      if (count % points_per_sec == 0)
      {
        data_buffer[count/points_per_sec + 1] = t;
        for (size_t i=0; i< x.size(); i++)
        {
          data_buffer[(i+1)*(num_points) + count/points_per_sec + 1] = x[i];
        }
        points_saved++;
      }
      count++;
    }

  }


  void param_forced_pend::set_state(double theta, double theta_dot)
  {
    x[0] = theta;
    x[1] = theta_dot;
  }

  void param_forced_pend::set_pend_params(double *params)
  {
    A = params[0];
    L = params[1];
    d = params[2];
    omega = params[3];
    b = params[4];
    m = params[5];
    k = params[6];
  }

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
  double param_forced_pend::f(double t, double theta)
    {
    assert(A + L <= d); // make sure the pendulum bob doesn't go below the magnet
    double z = cos(omega * t);
    return -0.5*(pow(L, 2) + pow(d, 2) + pow((A*z), 2) + 2*A*L*z*cos(theta) -
                                        2*A*d*z - 2*d*L*cos(theta));
    }

  /**
  * Operator overload. Right hand side of equations of motion for the Pendulum
  * (dxdt = f(x, t)).
  */
  // template<class State, class Deriv>
  void param_forced_pend::operator()(const state_type &x, state_type &dxdt, double t)
  {
    double z = cos(omega*t);

    // typename boost::range_iterator< const State >::type x = boost::begin( x_ );
    // typename boost::range_iterator< Deriv >::type dxdt = boost::begin( dxdt_ );

    dxdt[0] = x[1];
    dxdt[1] = -k*x[1] - (sin(x[0])/L)*(g - A*z*pow(omega, 2) +
              (b/m)*(A*z-d)*exp(f(t, x[0])));
  }

  /**
  *
  */
  double param_forced_pend::solve(double solve_params[5])
  {
    using namespace boost::numeric::odeint;
    typedef runge_kutta_fehlberg78<state_type> stepper_type;
    double dt = solve_params[0];
    double abs_err = solve_params[1];
    double rel_err = solve_params[2];
    double t_fin = solve_params[4];
    auto stepper = make_controlled(abs_err, rel_err, stepper_type());
    integrate_n_steps(stepper, boost::ref( *this ), x, 0.0, dt, (int) (t_fin/dt));
    return x[0];
  }


  /**
  *
  */
  double param_forced_pend::solve(double solve_params[5], std::ofstream &out)
  {
    using namespace boost::numeric::odeint;
    typedef runge_kutta_fehlberg78<state_type> stepper_type;
    double dt = solve_params[0];
    double abs_err = solve_params[1];
    double rel_err = solve_params[2];
    double t_fin = solve_params[4];
    auto stepper = make_controlled(abs_err, rel_err, stepper_type());

    integrate_n_steps(stepper, boost::ref( *this ), x, 0.0, dt, (int) (t_fin/dt), streaming_observer_txt(out));
    return x[0];
  }


  /**
  *
  */
  double param_forced_pend::solve(double solve_params[5], double *data)
  {
    using namespace boost::numeric::odeint;
    typedef runge_kutta_fehlberg78<state_type> stepper_type;
    double dt = solve_params[0];
    double abs_err = solve_params[1];
    double rel_err = solve_params[2];
    double t_fin = solve_params[3] + solve_params[4] - 1;
    auto stepper = make_controlled(abs_err, rel_err, stepper_type());

    integrate_n_steps(stepper, boost::ref( *this ), x, 0.0, dt, (int) (t_fin/dt), streaming_observer_arr(data, solve_params));
    return x[0];
  }
