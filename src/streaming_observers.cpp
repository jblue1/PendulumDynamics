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

streaming_observer_h5::streaming_observer_h5(double *data, const int y) :
  count(0),
  data_buffer(data),
  DY(y)
  {}

void streaming_observer_h5::operator()(const state_type &x, double t)
{
  data_buffer[count] = t;
  for (size_t i=0; i< x.size(); i++) {
    data_buffer[(i+1)*DY + count] = x[i];
  }
  count++;
}
