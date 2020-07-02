#ifndef OBSERVERSREF
#define OBSERVERSREF

#include <boost/numeric/odeint.hpp>
#include <fstream>

typedef std::vector<double> state_type;


struct streaming_observer_csv
{
  std::ofstream &write_out;
  streaming_observer_csv(std::ofstream &out);
  void operator()(const state_type &x, double t);

};

struct streaming_observer_h5
{
  int count;
  double *data_buffer;
  const int DY;

  streaming_observer_h5(double *data, const int y);

  /**
   * Operator overload called by odeint to write solutions to file
   */
  void operator()(const state_type &x, double t);
};

#endif
