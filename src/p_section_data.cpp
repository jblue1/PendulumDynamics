#include <boost/numeric/odeint.hpp>
#include <stdio.h>
#include <string>
#include <fstream>
#include <cmath>
#include <filesystem>
#include "H5Cpp.h"
#include "systems.hpp"
#include "streaming_observers.hpp"


using namespace boost::numeric::odeint;
using namespace H5;
namespace fs = std::filesystem;

// type of container used by ODE solver to store state of the system
typedef std::vector<double> state_type;

const double g = 9.81; //gravitational acceleration

/**
 * Describe parameters taken and print errors if paramters are incorrectly
 * given.
 * @param out output stream for messages
 * @param msg Error message to print
 */
static void usage(std::ostream &out, const char* msg){
    out << msg << "\n" << "\n";
    out << "  Usage:\n";
    out << "        solve_system file A\n";
    out << "  dir        - path to directory to save file to\n";
    out << "  time       - Time to simulate system\n";
    out << "  A_start    - starting driving amplitude\n";
    out << "  A_step     - driving amplitude stepsize\n";
    out << "  num_steps  - number of A values\n";
    exit(1);
}

static void write_params(std::ostream &out, const double (&params)[6]) {
  out << "L = " << params[0] << "\n";
  out << "d = " << params[1] << "\n";
  out << "omega = " << params[2] << "\n";
  out << "b = " << params[3] << "\n";
  out << "m = " << params[4] << "\n";
  out << "k = " << params[5] << "\n";
}



int main(int argc, char const *argv[]) {


  if (argc != 6) {
    usage(std::cerr, "Incorrect Number of parameters given.");
  }

  // store command line arguments
  const std::string dir_name(argv[1]);
  const double t_fin = atof(argv[2]);
  const double A_start = atof(argv[3]);
  const double A_step = atof(argv[4]);
  const int num_steps = atoi(argv[5]);

  // create directory to save data to
  if (!fs::is_directory(dir_name) || !fs::exists(dir_name)) {
    fs::create_directory(dir_name);
  }

  else {
    usage(std::cerr, "Provide a path that does not lead to an existing directory or file");
  }

  // define parameters for system
  const double L = g/pow(M_PI, 2);
  const double d = 4.0;
  const double omega = 2*M_PI;
  const double b = 50.0;
  const double m = 0.1;
  const double k = 0.2;
  const double params[6] = {L, d, omega, b, m, k};

  //write parameters out to file
  std::ofstream write_out(dir_name + "/params.txt");
  write_params(write_out, params);

  // define parameters for ODE solver
  const double abs_err = 1e-10;
  const double rel_err = 1e-10;
  const double points_per_sec = 100.0;
  const double dt = 1/points_per_sec;
  const int num_points = lrint(t_fin + 1);



  // Create HDF5 file
  const H5std_string FILE_NAME(dir_name + "/data.h5");
  H5File file(FILE_NAME, H5F_ACC_EXCL); // open will fail if file already exists


  // Constants needed to create dataspaces
  const hsize_t DX = 3*num_points;
  const int RANK = 1;
  hsize_t dataspace_dims[1] = {DX};


  // create vector dictating the times at which we want solutions
  std::vector<double> times(num_points);
  times[0] = 0.0;
  for( size_t i=1 ; i<times.size() ; ++i ){
    times[i] = i - 0.5; //storing data at every half second
  }

  typedef runge_kutta_fehlberg78<state_type> error_stepper_type;
  error_stepper_type stepper;

  double A;

  const double step_theta = 2*M_PI/199;
  const double step_theta_dot = 6.0/4.0;

  for (size_t i=0; i < num_steps; i++){

    A = (i)*A_step + A_start;
    std::cout << "A: " << A << "\n";
    std::string Astr = std::to_string(A);
    Astr.pop_back();
    Astr.pop_back();

    // Create group inside file
    Group group(file.createGroup("/group" + Astr));
    int count = 0;
    for (size_t j=0; j<100; j++) {
      for (size_t p=0; p<5; p++) {

        // create DataSpace
        DataSpace dataspace(RANK, dataspace_dims);

        //create dataset
        std::string countstr = std::to_string(count);
        std::string dset = "dset" + countstr;
        const H5std_string DATASET_NAME(dset);
        DataSet dataset = group.createDataSet(DATASET_NAME, PredType::NATIVE_DOUBLE, dataspace);

        //instantiate state vector
        state_type x(2);
        x[0] = -M_PI + step_theta*j;
        x[1] = -3 + step_theta_dot*p;
        //instantiate data array
        double *data = new double [3 * (num_points + 1)];

        integrate_times(make_controlled( abs_err , rel_err , error_stepper_type() ),
                        param_forced_pend(A, L, d, omega, b, m, k),
                        x , times, dt , streaming_observer_h5(data, num_points));
        dataset.write(data, PredType::NATIVE_DOUBLE);
        count++;
        //free memory of data array
        delete [] data;
      }
    }
  }

  return 0;
}
