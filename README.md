# PendulumDynamics
This is a library used to simulate the motion of a parametrically driven pendulum with external forcing. See TREND_Poster.pdf for a summary of the work.
### Tested on:
 - MacOS 10.15.4
 - Ubuntu 20.04 LTS
 
 ### Requirements
  - [BOOST](https://www.boost.org/)
  - [HDF5](https://www.hdfgroup.org/)
  - [OpenMPI](https://www.open-mpi.org/)
  
In order for the Makefile to work, you need to set the following environmental variables (assuming mpic++ is your mpi wrapper)

```Makefile
export BOOST_ROOT=path/to/your/boost/installation
export HDF5_CXX=mpic++
export HDF5_CLINKER=mpic++
```
