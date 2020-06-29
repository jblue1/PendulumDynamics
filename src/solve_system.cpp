#include <boost/numeric/odeint.hpp>
#include <stdio.h>
#include <string>
#include <fstream>
#include <cmath>
#include <filesystem>
#include "H5Cpp.h"


using namespace boost::numeric::odeint;
using namespace H5;
namespace fs = std::filesystem;

// type of container used by ODE solver to store state of the system
typedef std::vector<double> state_type;

const double g = 9.81; //gravitational acceleration



/**
 * Class containing equations and parameters associated with a parametrically
 * forced pendulum oscillating vertically about the origin with a repulsive
 * magnet located some distance below along the y axis.
 */
class param_forced_pend {
protected:
  double A; // Amplitude of pivot oscillations
  double L; // Length of pendulum
  double d; // Distance from origin to magnet
  double omega; // Angular frequency of oscillations
  double b; // Controls strength of magnet
  double m; // Mass of pendulum
  double k; // Linear friction coefficient


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
  double f(double t, double theta, double A, double L, double d, double omega){
    double z = cos(omega * t);
    return -0.5*(pow(L, 2) + pow(d, 2) + pow((A*z), 2) + 2*A*L*z*cos(theta) -
            2*A*d*z - 2*d*L*cos(theta));
  }

public:
  /**
   * Class constructor. Params described above.
   */
  param_forced_pend(double A, double L, double d, double omega, double b,
                    double m, double k){
    this->A=A;
    this->L=L;
    this->d=d;
    this->omega=omega;
    this->b=b;
    this->m=m;
    this->k=k;
  }

  /**
   * Operator overload. Right hand side of equations of motion for the Pendulum
   * (dxdt = f(x, t)).
   */
  void operator()(state_type &x, state_type &dxdt, double t){
    double z = cos(omega*t);
    dxdt[0] = x[1];
    dxdt[1] = -k*x[1] - (sin(x[0])/L)*(g - A*z*pow(omega, 2) +
              (b/m)*(A*z-d)*exp(f(t, x[0], A, L, d, omega)));

  }
};

/**
 * Structure used by ODE solver to write out solutions each step.
 */
struct streaming_observer {
  int count;
  double *data_buffer;
  const int DY;

  streaming_observer(double *data, const int y) :
                                                 data_buffer(data),
                                                 count(0),
                                                 DY(y)
                                                 {} // constructor

  /**
   * Operator overload called by odeint to write solutions to file
   */
  void operator()(const state_type &x, double t)  {
    data_buffer[count] = t;
    for (size_t i=0; i< x.size(); i++) {
      data_buffer[(i+1)*DY + count] = x[i];
    }
    count++;
  }
};

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
    out << "  dir   - path to directory to save file to\n";
    out << "  time  - Time to simulate system\n";
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


  if (argc != 3) {
    usage(std::cerr, "Incorrect Number of parameters given.");
  }

  // store command line arguments
  const std::string dir_name(argv[1]);
  const double t_fin = atof(argv[2]);

  std::cout << "Here" << "\n";

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
  const int num_points = lrint(t_fin * points_per_sec);



  // Create HDF5 file
  const H5std_string FILE_NAME(dir_name + "/data.h5");
  H5File file(FILE_NAME, H5F_ACC_EXCL); // open will fail if file already exists


  // Constants needed to create dataspaces
  const hsize_t DX = 3*num_points;
  const int RANK = 1;
  hsize_t dataspace_dims[1] = {DX};


  // Constants needed to create attributes
  const H5std_string ATTR_NAME("Initial Conditions");
  hsize_t attr_dims[1] = {2};


  // create vector dictating the times at which we want solutions
  std::vector<double> times(num_points );
  for( size_t i=0 ; i<times.size() ; ++i ){
    times[i] = dt*i;
  }

  typedef runge_kutta_fehlberg78<state_type> error_stepper_type;
  error_stepper_type stepper;


  const double A_step = 0.0001;
  double A;

  const double step_theta = 2*M_PI/199;
  const double step_theta_dot = 6.0/4.0;

  for (size_t i=0; i < 100; i++){

    A = i*A_step + 0.01;
    std::cout << "A: " << A << "\n";
    std::string Astr = std::to_string(A);
    Astr.pop_back();
    Astr.pop_back();

    // Create group inside file
    Group group(file.createGroup("/group" + Astr));
    int count = 0;
    for (size_t j=0; j<100; j++) {
      for (size_t k=0; k<5; k++) {

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
        x[1] = -3 + step_theta_dot*k;

        // create attribute to store initial Conditions
        // stored as array with [theta_init, theta_dot_init]
        double attr_data[2] = {x[0], x[1]};
        DataSpace attr_dataspace = DataSpace(1, attr_dims);
        Attribute attribute = dataset.createAttribute(ATTR_NAME, PredType::NATIVE_DOUBLE, attr_dataspace);
        attribute.write(PredType::NATIVE_DOUBLE, attr_data);

        // instantiate data array
        double *data = new double [3 * num_points];

        integrate_times(make_controlled( abs_err , rel_err , error_stepper_type() ),
                        param_forced_pend(A, L, d, omega, b, m, k),
                        x , times , dt , streaming_observer(data, num_points));
        dataset.write(data, PredType::NATIVE_DOUBLE);
        count++;
        //free memory of data array
        delete [] data;
      }
    }
  }

  return 0;
}
