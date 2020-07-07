#include "streaming_observers.hpp"
#include <fstream>

typedef std::vector<double> state_type;

streaming_observer_csv::streaming_observer_csv(std::ofstream &out):
  write_out(out){}

void streaming_observer_csv::operator()(const state_type &x, double t)
{
  write_out << t;
  for (size_t i=0; i < x.size(); i++) {
    write_out << "," << x[i];
  }
  write_out << "\n";
}

// constructor
streaming_observer_h5::streaming_observer_h5(double *data, const int y) :
  count(0),
  data_buffer(data),
  DY(y) // y dimension of data array
  {}

/**
* Operator overload used by odeint to write out data. Becuase HDF5 can only
* write data stored contigiously in memory, solutions from odeint which
* would more easily be stored in a N x 3 boost:array (the columns being time,
* x1 and x2) is instead stored in a 3N x 1 array.
*/
void streaming_observer_h5::operator()(const state_type &x, double t)
{
  data_buffer[count] = t;
  for (size_t i=0; i< x.size(); i++) {
    data_buffer[(i+1)*DY + count] = x[i];
  }
  count++;
}
