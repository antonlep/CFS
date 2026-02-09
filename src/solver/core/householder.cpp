#include "householder.h"
#include <Eigen/Dense>
#include <iostream>
#include <vector>

int sign(double x) {
  if (x < 0)
    return -1;
  return 1;
};

Eigen::VectorXd householder(Eigen::MatrixXd K_init, Eigen::VectorXd F) {
  std::vector<Eigen::MatrixXd> H_list;
  const int n = K_init.rows();
  int r = K_init.rows();
  int c = K_init.cols();
  Eigen::MatrixXd K = K_init;
  for (int i = c; i > 0; i--) {
    Eigen::VectorXd k = K.col(0);
    Eigen::VectorXd v =
        k + sign(k(0)) * k.norm() * Eigen::MatrixXd::Identity(r, 1);
    Eigen::MatrixXd H = Eigen::MatrixXd::Identity(r, r) -
                        2 * v * v.transpose() / (v.transpose() * v);
    H_list.push_back(H);
    Eigen::MatrixXd HK = H * K;
    K = HK.block(1, 1, r - 1, c - 1);
    r = K.rows();
    c = K.cols();
  }

  Eigen::MatrixXd R = Eigen::MatrixXd::Identity(n, n);
  for (auto it = H_list.rbegin(); it != H_list.rend(); ++it) {
    const Eigen::MatrixXd &h_temp = *it;
    int offset = n - h_temp.rows();
    Eigen::MatrixXd h = Eigen::MatrixXd::Identity(n, n);
    h.block(offset, offset, h_temp.rows(), h_temp.cols()) = h_temp;
    R = R * h;
  }
  R *= K_init;

  Eigen::MatrixXd Q = Eigen::MatrixXd::Identity(n, n);
  for (Eigen::MatrixXd h_temp : H_list) {
    int offset = n - h_temp.rows();
    Eigen::MatrixXd h = Eigen::MatrixXd::Identity(n, n);
    h.block(offset, offset, h_temp.rows(), h_temp.cols()) = h_temp;
    Q = Q * h;
  }
  Eigen::MatrixXd y = Q.transpose() * F;

  Eigen::VectorXd u(n);
  for (int i = n - 1; i >= 0; i--) {
    double sum = y(i);
    for (int j = i + 1; j < n; j++) {
      sum -= R(i, j) * u(j);
    }
    u(i) = sum / R(i, i);
  }

  return u;
}