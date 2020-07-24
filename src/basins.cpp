#include "H5Cpp.h"
#include "systems.hpp"
#include <assert.h>
#include <cmath>
#include <filesystem>
#include <mpi.h>
#include <string>

using namespace H5;
namespace fs = std::filesystem;

const double g = 9.81; // gravitational acceleration

/**
 * Describe parameters taken and print errors if parameters are incorrectly
 * given.
 * @param out output stream for messages
 * @param msg Error message to print
 */
static void usage(std::ostream &out, const char *msg) {
  out << msg << "\n"
      << "\n";
  out << "  Usage:\n";
  out << "        solve_system dir time A_start A_step num_steps\n";
  out << "  dir     - path to directory to save file to\n";
  out << "  A       - Driving amplitude\n";
  out << "  Time    - Time to simulate system\n";
  exit(1);
}

/**
 * Write parameters out to a file
 * @param  out      out stream for file to write to
 * @param  params   array with parameter values
 */
static void write_params(std::ostream &out, const double (&params)[12]) {
  out << "L = " << params[0] << "\n";
  out << "d = " << params[1] << "\n";
  out << "omega = " << params[2] << "\n";
  out << "b = " << params[3] << "\n";
  out << "m = " << params[4] << "\n";
  out << "k = " << params[5] << "\n";
  out << "A = " << params[6] << "\n";
  out << "dt = " << params[7] << "\n";
  out << "abs_err = " << params[8] << "\n";
  out << "rel_err = " << params[9] << "\n";
  out << "trans_time = " << params[10] << "\n";
  out << "simul_time = " << params[11] << "\n";
}

int main(int argc, char *argv[]) {

  if (argc != 4) {
    usage(std::cerr, "Incorrect Number of parameters given.");
  }

  // store command line arguments
  const std::string dir_name(argv[1]);
  const double A = atof(argv[2]);
  const double simul_time = atof(argv[3]);

  // define parameters for system
  const double L = g / pow(M_PI, 2);
  const double d = 4.0;
  const double omega = 2 * M_PI;
  const double b = 50.0;
  const double m = 0.1;
  const double k = 0.2;
  double pend_params[7] = {A, L, d, omega, b, m, k};

  // define parameters for ODE solver
  const double dt = 0.1;
  const double abs_err = 1e-10;
  const double rel_err = 1e-10;
  const double trans_time = 5000.0;

  // create directory for saving data
  if (!fs::is_directory(dir_name) || !fs::exists(dir_name)) {
    fs::create_directory(dir_name);
  }

  else {
    usage(std::cerr,
          "Provide a path that does not lead to an existing directory or file");
  }
  const double params[12] = {L, d,  omega,   b,       m,          k,
                             A, dt, abs_err, rel_err, trans_time, simul_time};

  // write parameters out to file
  std::ofstream write_out(dir_name + "/params.txt");
  write_params(write_out, params);
  write_out.close();

  // Constants needed to create dataspaces
  int num_points = (int)(simul_time) + 1;
  const hsize_t DX = 3 * num_points;
  const int RANK = 1;
  hsize_t dataspace_dims[1] = {DX};

  const double step_theta = 2 * M_PI / 103.0;
  const double step_theta_dot = 6.0 / 99.0;

  // Create HDF5 file
  const H5std_string FILE_NAME(dir_name + "/data.h5");
  H5File file(FILE_NAME, H5F_ACC_EXCL); // open will fail if file already exists
  int count = 0;
  for (size_t j = 0; j < 104; j++) {
    std::cout << "\n";
    std::cout << "\n";
    std::cout << "j: " << j << "\n";
    std::cout << "\n";
    std::cout << "\n";
    for (size_t p = 0; p < 100; p++) {
      std::cout << "p: " << p << "\n";

      // create DataSpace
      DataSpace dataspace(RANK, dataspace_dims);

      // create dataset
      std::string countstr = std::to_string(count);
      std::string dset = "dset" + countstr;
      const H5std_string DATASET_NAME(dset);
      DataSet dataset = file.createDataSet(
          DATASET_NAME, PredType::NATIVE_DOUBLE, dataspace);

      // instantiate pendulum object
      param_forced_pend pend(pend_params);

      // set initial state
      pend.set_state(-M_PI + step_theta * j, -3 + step_theta_dot * p);

      // instantiate data array
      double *data = new double[3 * (num_points)];

      // solve system
      pend.solve(dt, abs_err, rel_err, trans_time, simul_time, data);

      dataset.write(data, PredType::NATIVE_DOUBLE);
      dataset.close();

      // free memory of data array
      delete[] data;
      count++;
    }
  }



  return 0;
}
