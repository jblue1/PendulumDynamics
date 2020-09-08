#include "../src/systems.hpp"
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE systems
#include <boost/array.hpp>
#include <boost/test/unit_test.hpp>
#include <cmath>

// test helper function f
BOOST_AUTO_TEST_CASE(f_test) {
  // true values calculated with mathematica
  double pend_params[7] = {0.5, 1.0, 2.0, 2.0, 0.0, 0.0, 0.0};
  param_forced_pend pendulum(pend_params);
  double b = pendulum.f(0.0, M_PI);
  BOOST_CHECK_CLOSE(b, -25.0 / 8.0, 1e-12f);
  b = pendulum.f(0.0, -M_PI);
  BOOST_CHECK_CLOSE(b, -25.0 / 8.0, 1e-12f);
  b = pendulum.f(1.0, -M_PI);
  BOOST_CHECK_CLOSE(b, -5.145867528516738, 1e-12f);
}

BOOST_AUTO_TEST_CASE(df2dtheta) {
  double pend_params[7] = {0.5, 1.0, 2.0, 2.0, 50.0, 0.1, 0.2};
  double t = 0.0;
  pendulum_lyap pendulum(pend_params);
  double b = pendulum.df2dtheta(0, t * 2.0);
  BOOST_CHECK_CLOSE(b, 654.0626769384467, 1e-12);
  pendulum.set_state(M_PI, 0, 0);
  t = 1.0;
  b = pendulum.df2dtheta(M_PI, t * 2.0);
  BOOST_CHECK_CLOSE(b, 4.213024011792606, 1e-12);
}

BOOST_AUTO_TEST_CASE(df2dphi) {
  double pend_params[7] = {0.5, 1.0, 2.0, 2.0, 50.0, 0.1, 0.2};
  double t = 0.0;
  pendulum_lyap pendulum(pend_params);
  double b = pendulum.df2dphi(0, t * 2.0);
  BOOST_CHECK_CLOSE(b, 0, 1e-12);
  pendulum.set_state(M_PI, 0, 0);
  t = 1.0;
  b = pendulum.df2dphi(M_PI / 2, t * 2.0);
  BOOST_CHECK_CLOSE(b, -48.49698041762917, 1e-12);
}

BOOST_AUTO_TEST_CASE(set_state) {
  double pend_params[7] = {0.5, 1.0, 2.0, 2.0, 50.0, 0.1, 0.2};
  pendulum_lyap pendulum(pend_params);
  double x0 = -1;
  double x1 = 3;
  double x2 = 0;
  pendulum.set_state(x0, x1, x2);
  boost::array<double, 12> x = pendulum.get_state();
  BOOST_TEST(x[0] == x0);
  BOOST_TEST(x[1] == x1);
  BOOST_TEST(x[2] == x2);
}

BOOST_AUTO_TEST_CASE(init_pert_vec) {
  double pend_params[7] = {0.5, 1.0, 2.0, 2.0, 50.0, 0.1, 0.2};
  pendulum_lyap pendulum(pend_params);
  pendulum.set_state(1.0, 1.0, 5.0);
  pendulum.init_pert_vecs();
  boost::array<double, 12> x = pendulum.get_state();
  BOOST_TEST(x[0] == 1.0);
  BOOST_TEST(x[1] == 1.0);
  BOOST_TEST(x[2] == 5.0);
  BOOST_TEST(x[3] == 1);
  BOOST_TEST(x[7] == 1);
  BOOST_TEST(x[11] == 1);
  BOOST_TEST(x[10] == 0);
}

BOOST_AUTO_TEST_CASE(normalize) {
  double pend_params[7] = {0.5, 1.0, 2.0, 2.0, 50.0, 0.1, 0.2};
  pendulum_lyap pendulum(pend_params);
  boost::array<double, 3> vec = {1, 2, 2};
  boost::array<double, 3>::iterator first = vec.begin();
  boost::array<double, 3>::iterator last = vec.end();
  double norm = sqrt(std::inner_product(first, last, first, 0.0));
  pendulum.normalize(first, last, norm);
  BOOST_CHECK_CLOSE(vec[0], 1.0 / 3.0, 1e-12);
  BOOST_CHECK_CLOSE(vec[1], 2.0 / 3.0, 1e-12);
  BOOST_CHECK_CLOSE(vec[2], 2.0 / 3.0, 1e-12);
}

BOOST_AUTO_TEST_CASE(gram_schmidt) {
  double pend_params[7] = {0.5, 1.0, 2.0, 2.0, 50.0, 0.1, 0.2};
  pendulum_lyap pendulum(pend_params);
  boost::array<double, 12> state1 = {0.0,  0.0, 0.0, 1.0, 2.0, 2.0,
                                     -1.0, 0.0, 2.0, 0.0, 0.0, 1.0};
  pendulum.set_state(state1);
  pendulum.gram_schmidt();
  state1 = pendulum.get_state();
  boost::array<double, 12> state2 = {
      0.0,        0.0,        0.0,       1.0 / 3.0, 2.0 / 3.0,  2.0 / 3.0,
      -2.0 / 3.0, -1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0, -2.0 / 3.0, 1.0 / 3.0};
  for (int i = 3; i < 12; i++) {
    BOOST_TEST_CONTEXT("index " << i) {
      BOOST_CHECK_CLOSE(state1[i], state2[i], 1e-12);
    }
  }

  state1 = {0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 2.0, 0.0, 2.0, 0.0, 1.0};
  state2 = {0.0,
            0.0,
            0.0,
            1.0 / sqrt(3.0),
            1.0 / sqrt(3.0),
            1.0 / sqrt(3.0),
            0.0,
            1.0 / sqrt(2.0),
            -1.0 / sqrt(2.0),
            2.0 / sqrt(6.0),
            -1.0 / sqrt(6.0),
            -1.0 / sqrt(6.0)};

  pendulum.set_state(state1);
  pendulum.gram_schmidt();
  state1 = pendulum.get_state();

  for (int i = 3; i < 12; i++) {
    BOOST_TEST_CONTEXT("index " << i) {
      BOOST_CHECK_CLOSE(state1[i], state2[i], 1e-12);
    }
  }


}

// test operator overload
BOOST_AUTO_TEST_CASE(derivatives) {
  boost::array<double, 2> x = {0, 0};
  boost::array<double, 2> dxdt = {5.0, 5.0};
  double t = 0.0;
  double pend_params[7] = {0.5, 1.0, 2.0, 2.0, 50.0, 0.1, 0.2};
  param_forced_pend pendulum(pend_params);
  pendulum(x, dxdt, t);
  BOOST_CHECK_SMALL(dxdt[0], 1e-10);
  BOOST_CHECK_SMALL(dxdt[1], 1e-10);

  x[0] = 1.0;
  x[1] = 10.0;
  pendulum(x, dxdt, t);
  BOOST_CHECK_CLOSE(dxdt[0], x[1], 1e-12f);
  BOOST_CHECK_CLOSE(dxdt[1], 270.9059275869695, 1e-12f);

  t = M_PI;
  pendulum(x, dxdt, t);
  BOOST_CHECK_CLOSE(dxdt[0], x[1], 1e-12f);
  BOOST_CHECK_CLOSE(dxdt[1], 270.9059275869695, 1e-12f);

  x[0] = -1.0 - 2 * M_PI;
  x[1] = -10.0;
  pendulum(x, dxdt, t);
  BOOST_CHECK_CLOSE(dxdt[0], x[1], 1e-12f);
  BOOST_CHECK_CLOSE(dxdt[1], -270.9059275869695, 1e-12f);
}

// test Integration
BOOST_AUTO_TEST_CASE(solve) {
  // check undamped case
  double A = 0;
  double L = 1.0;
  double d = 2.0;
  double omega = 3.0;
  double b = 0.0;
  double m = 0.1;
  double k = 0;
  double pend_params[7] = {A, L, d, omega, b, m, k};

  double theta_init = 0.001;
  double theta_dot_init = 0.001;
  double dt = 0.1;
  double abs_err = 1e-10;
  double rel_err = 1e-10;
  double t_fin = 10.0;
  double omega_0 = sqrt(9.81 / L);
  double theta_fin = theta_init * cos(omega_0 * t_fin) +
                     (theta_dot_init / omega_0) * sin(omega_0 * t_fin);

  param_forced_pend pend(pend_params);
  pend.set_state(theta_init, theta_dot_init);

  double x = pend.solve(dt, abs_err, rel_err, t_fin);
  BOOST_CHECK_CLOSE(theta_fin, x, 0.01f);

  // check damped case
  k = 0.2;
  pend_params[6] = k;
  t_fin = 5.0;
  pend.set_pend_params(pend_params);
  pend.set_state(theta_init, theta_dot_init);
  double beta = k / 2.0;
  double omega_1 = sqrt(pow(omega_0, 2) - pow(beta, 2));

  double theta_fin2 =
      exp(-beta * t_fin) * (theta_init * cos(omega_1 * t_fin) +
                            (theta_dot_init / omega_1) * sin(omega_1 * t_fin));
  x = pend.solve(dt, abs_err, rel_err, t_fin);

  BOOST_CHECK_CLOSE(theta_fin2, x, 0.5f);
}

BOOST_AUTO_TEST_CASE(solve_lyap) {
  double A = 0.01;
  double L = 1.0;
  double d = 2.0;
  double omega = 3.0;
  double b = 10.0;
  double m = 0.1;
  double k = 0.2;
  double pend_params[7] = {A, L, d, omega, b, m, k};
  double theta_init = 1.0;
  double theta_dot_init = 2.0;
  double dt = 0.1;
  double abs_err = 1e-10;
  double rel_err = 1e-10;
  double t_fin = 100.0;

  param_forced_pend pf_pend(pend_params);
  pendulum_lyap lyap_pend(pend_params);

  pf_pend.set_state(theta_init, theta_dot_init);
  lyap_pend.set_state(theta_init, theta_dot_init, 0);

  double theta_fin1 = pf_pend.solve(dt, abs_err, rel_err, t_fin);
  lyap_pend.set_state(theta_init, theta_dot_init, 1);
  double theta_fin2 = lyap_pend.solve(dt, abs_err, rel_err, t_fin);

  BOOST_CHECK_CLOSE(theta_fin1, theta_fin2, 1e-12);
}

// test that the different streaming_observer_txt works
BOOST_AUTO_TEST_CASE(streaming_observers_txt) {
  double A = 0.1;
  double L = 1.0;
  double d = 2.0;
  double omega = 3.0;
  double b = 10.0;
  double m = 0.1;
  double k = 0.2;
  double pend_params[7] = {A, L, d, omega, b, m, k};
  double theta_init = 1.0;
  double theta_dot_init = 2.0;
  double dt = 0.1;
  double abs_err = 1e-10;
  double rel_err = 1e-10;
  double t_fin = 3500.0;

  int num_data_points = (int)t_fin / dt;

  param_forced_pend pend(pend_params);
  pend.set_state(theta_init, theta_dot_init);
  double x1 = pend.solve(dt, abs_err, rel_err, t_fin);

  pend.set_state(theta_init, theta_dot_init);
  std::ofstream write_out("test.txt");
  double x2 = pend.solve(dt, abs_err, rel_err, t_fin, write_out);
  write_out.close();

  int count = 0;
  double x, y, z;
  double t;
  std::ifstream read_txt("test.txt");
  while (read_txt >> x >> y >> z) {
    if (count == 35000) {
      t = x;
    }
    count++;
  }
  BOOST_CHECK_CLOSE(t, t_fin, 1e-10); // make sure final time is same as t_fin
  BOOST_TEST(num_data_points ==
             count - 1); // make sure correct number of rows written
  read_txt.close();

  BOOST_CHECK_CLOSE(x1, x2, 1e-12); // make sure solutions are the same
  remove("test.txt");
}

// make sure streaming observer_txt works
BOOST_AUTO_TEST_CASE(streaming_observer_arr) {
  double A = 0.1;
  double L = 1.0;
  double d = 2.0;
  double omega = 3.0;
  double b = 10.0;
  double m = 0.1;
  double k = 0.2;
  double pend_params[7] = {A, L, d, omega, b, m, k};
  double theta_init = 1.0;
  double theta_dot_init = 2.0;
  double dt = 0.1;
  double abs_err = 1e-10;
  double rel_err = 1e-10;
  double trans_time = 3000.0;
  double simul_time = 500.0;
  int num_data_points = (int)(simul_time) + 1;

  param_forced_pend pend(pend_params);
  pend.set_state(theta_init, theta_dot_init);
  double x1 = pend.solve(dt, abs_err, rel_err, trans_time + simul_time - 1);

  double *data1 = new double[3 * (num_data_points)];
  pend.set_state(theta_init, theta_dot_init);
  double x2 = pend.solve(dt, abs_err, rel_err, trans_time, simul_time, data1);

  // make sure initial conditions are being recorded
  BOOST_TEST(data1[0] == 0.0);
  BOOST_TEST(data1[num_data_points] == theta_init);
  BOOST_TEST(data1[2 * num_data_points] == theta_dot_init);

  // make sure streaming observer is recording the correct solution
  BOOST_CHECK_CLOSE(data1[11], 3010.0, 1e-12);
  BOOST_CHECK_CLOSE(x1, data1[2 * num_data_points - 1], 1e-12);
  delete[] data1;

  // try with different dt, trans_time and simul_time
  dt = 0.1;
  trans_time = 3000.0;
  simul_time = 4000.0;
  num_data_points = (int)(simul_time) + 1;

  pend.set_state(theta_init, theta_dot_init);
  x1 = pend.solve(dt, abs_err, rel_err, trans_time + simul_time - 1);

  double *data2 = new double[3 * num_data_points];

  pend.set_state(theta_init, theta_dot_init);
  x2 = pend.solve(dt, abs_err, rel_err, trans_time, simul_time, data2);

  // make sure initial conditions are being recorded
  BOOST_TEST(data2[0] == 0.0);
  BOOST_TEST(data2[num_data_points] == theta_init);
  BOOST_TEST(data2[2 * num_data_points] == theta_dot_init);

  // make sure streaming observer is recording the correct solution
  BOOST_CHECK_CLOSE(data2[11], 3010, 1e-12);
  BOOST_CHECK_CLOSE(x1, data2[2 * num_data_points - 1], 1e-12);

  delete[] data2;
}
