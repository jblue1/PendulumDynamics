#include "H5Cpp.h"
#include "systems.hpp"
#include <cmath>
#include <filesystem>
#include <assert.h>
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
  out << "  dir        - path to directory to save file to\n";
  out << "  time       - Time to simulate system\n";
  out << "  A_start    - starting driving amplitude\n";
  out << "  A_step     - driving amplitude stepsize\n";
  out << "  num_steps  - number of A values (if you are using parallel "
         "processing\n";
  out << "               this needs to be a multiple of the number of "
         "procceses)\n";
  exit(1);
}

/**
 * Write parameters out to a file
 * @param  out      out stream for file to write to
 * @param  params   array with parameter values
 */
static void write_params(std::ostream &out, const double (&params)[14]) {
  out << "L = " << params[0] << "\n";
  out << "d = " << params[1] << "\n";
  out << "omega = " << params[2] << "\n";
  out << "b = " << params[3] << "\n";
  out << "m = " << params[4] << "\n";
  out << "k = " << params[5] << "\n";
  out << "A_start = " << params[6] << "\n";
  out << "A_step = " << params[7] << "\n";
  out << "num_steps = " << params[8] << "\n";
  out << "dt = " << params[9] << "\n";
  out << "abs_err = " << params[10] << "\n";
  out << "rel_err = " << params[11] << "\n";
  out << "trans_time = " << params[12] << "\n";
  out << "simul_time = " << params[13] << "\n";
}

int main(int argc, char *argv[]) {

  MPI_Init(&argc, &argv);
  int num_procs;
  int rank;
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::cout << "Hello world from process " << rank + 1 << " of " << num_procs
            << "\n";

  if (rank == 0) {
    if (argc != 6) {
      usage(std::cerr, "Incorrect Number of parameters given.");
    }
  }

  // store command line arguments
  const std::string dir_name(argv[1]);
  const double simul_time = atof(argv[2]);
  const double A_start = atof(argv[3]);
  const double A_step = atof(argv[4]);
  const int num_steps = atoi(argv[5]);

  if (rank == 0) {
    if (num_steps % num_procs != 0) {
      usage(std::cerr,
            "Please choose a num_steps divisible by the number of processes.");
    }
  }

  // define parameters for system
  const double L = g / pow(M_PI, 2);
  const double d = 4.0;
  const double omega = 2 * M_PI;
  const double b = 50.0;
  const double m = 0.1;
  const double k = 0.2;
  double A = 0.0;
  double pend_params[7] = {A, L, d, omega, b, m, k};

  // define parameters for ODE solver
  const double dt = 0.1;
  const double abs_err = 1e-10;
  const double rel_err = 1e-10;
  const double trans_time = 3000.0;

  if (rank == 0) {
    // create directory to save data to
    if (!fs::is_directory(dir_name) || !fs::exists(dir_name)) {
      fs::create_directory(dir_name);
    }

    else {
      usage(
          std::cerr,
          "Provide a path that does not lead to an existing directory or file");
    }
    const double params[14] = {L,          d,         omega,
                               b,          m,         k,
                               A_start,    A_step,    double(num_steps),
                               dt,         abs_err,   rel_err,
                               trans_time, simul_time};

    // write parameters out to file
    std::ofstream write_out(dir_name + "/params.txt");
    write_params(write_out, params);
    write_out.close();
  }

  // Constants needed to create dataspaces
  int num_points = (int)(simul_time) + 1;
  const hsize_t DX = 3 * num_points;
  const int RANK = 1;
  hsize_t dataspace_dims[1] = {DX};

  const double step_theta = 2 * M_PI / 99.0;
  const double step_theta_dot = 6.0 / 4.0;

  MPI_Barrier(MPI_COMM_WORLD);

  // Create HDF5 file
  std::string rankstr = std::to_string(rank);
  const H5std_string FILE_NAME(dir_name + "/proc" + rankstr + "data.h5");
  H5File file(FILE_NAME, H5F_ACC_EXCL); // open will fail if file already exists

  int interval = num_steps / num_procs;
  for (int i = rank * interval; i < (rank + 1) * interval; i++) {
    A = (i)*A_step + A_start;
    pend_params[0] = A;
    std::cout << "A: " << A << "\n";
    std::string istr = std::to_string(i);

    // Create group inside file
    Group group(file.createGroup("/group" + istr));
    int count = 0;
    for (size_t j = 0; j < 100; j++) {
      for (size_t p = 0; p < 5; p++) {

        // create DataSpace
        DataSpace dataspace(RANK, dataspace_dims);

        // create dataset
        std::string countstr = std::to_string(count);
        std::string dset = "dset" + countstr;
        const H5std_string DATASET_NAME(dset);
        DataSet dataset = group.createDataSet(
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
        count++;

        // free memory of data array
        delete[] data;
      }
    }
  }
  MPI_Finalize();
  return 0;
}
