#include "Eigen/Core"
#include <Eigen/Dense>
#include <iostream>
using Eigen::MatrixXd;
using Eigen::VectorXd;

struct Node {
  float x;
  float y;
};

float area(Node i, Node j, Node k) {
  return 0.5 * (i.x * (j.y - k.y) + j.x * (k.y - i.y) + k.x * (i.y - j.y));
}

MatrixXd stiffness_matrix(float E, float v, float t, Node nodes[3]) {
  Node i = nodes[0];
  Node j = nodes[1];
  Node k = nodes[2];
  float A = area(i, j, k);

  MatrixXd B(3, 6);
  MatrixXd D(3, 3);
  MatrixXd K(6, 6);

  B << j.y - k.y, 0, k.y - i.y, 0, i.y - j.y, 0, 0, k.x - j.x, 0, i.x - k.x, 0,
      j.x - i.x, k.x - j.x, j.y - k.y, i.x - k.x, k.y - i.y, j.x - i.x,
      i.y - j.y;
  B *= 1 / (2 * A);

  D << 1, v, 0, v, 1, 0, 0, 0, (1 - v) / 2;
  D *= E / (1 - v * v);

  K = t * A * B.transpose() * D * B;
  return K;
}

MatrixXd update_global_stiffness(MatrixXd K, MatrixXd k, int order[3]) {
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      K(order[i] * 2, order[j] * 2) += k(i * 2, j * 2);
      K(order[i] * 2, order[j] * 2 + 1) += k(i * 2, j * 2 + 1);
      K(order[i] * 2 + 1, order[j] * 2) += k(i * 2 + 1, j * 2);
      K(order[i] * 2 + 1, order[j] * 2 + 1) += k(i * 2 + 1, j * 2 + 1);
    }
  }
  return K;
}

VectorXd element_stress(float E, float v, float t, Node nodes[3], int order[3],
                        VectorXd u_all) {
  std::vector<int> ind_order;
  for (int i = 0; i < 3; i++) {
    ind_order.push_back(order[i] * 2);
    ind_order.push_back(order[i] * 2 + 1);
  };
  VectorXd u = u_all(ind_order);

  Node i = nodes[0];
  Node j = nodes[1];
  Node k = nodes[2];
  float A = area(i, j, k);

  MatrixXd B(3, 6);
  MatrixXd D(3, 3);
  MatrixXd K(6, 6);

  B << j.y - k.y, 0, k.y - i.y, 0, i.y - j.y, 0, 0, k.x - j.x, 0, i.x - k.x, 0,
      j.x - i.x, k.x - j.x, j.y - k.y, i.x - k.x, k.y - i.y, j.x - i.x,
      i.y - j.y;
  B *= 1 / (2 * A);

  D << 1, v, 0, v, 1, 0, 0, 0, (1 - v) / 2;
  D *= E / (1 - v * v);

  VectorXd s = D * B * u;
  return s;
}

int add(int a, int b) { return a + b; }

float plate(float h, float w, float t, float f) {
  float E = 210E6;
  float v = 0.3;
  f = 0.5 * f;
  Node p1 = {0, 0};
  Node p2 = {0, 10};
  Node p3 = {20, 10};
  Node p4 = {20, 0};
  Node nodes1[3] = {p1, p3, p2};
  Node nodes2[3] = {p1, p4, p3};
  int order1[3] = {0, 2, 1};
  int order2[3] = {0, 3, 2};
  int size = 4;
  MatrixXd K = MatrixXd::Zero(size * 2, size * 2);

  MatrixXd K1 = stiffness_matrix(E, v, t, nodes1);
  K = update_global_stiffness(K, K1, order1);

  MatrixXd K2 = stiffness_matrix(E, v, t, nodes2);
  K = update_global_stiffness(K, K2, order2);

  std::vector<int> ind_nonzero{4, 5, 6, 7};
  MatrixXd K_nonzero;
  K_nonzero = K(ind_nonzero, ind_nonzero);

  VectorXd F(4);
  F << f, 0, f, 0;
  VectorXd u_nonzero(size);
  u_nonzero = K_nonzero.colPivHouseholderQr().solve(F);
  VectorXd u_zero(size);
  u_zero = VectorXd::Zero(size);
  VectorXd u(size * 2);
  u << u_zero, u_nonzero;
  std::cout << "u" << std::endl;
  std::cout << u << std::endl;
  std::cout << std::endl;

  VectorXd S1 = element_stress(E, v, t, nodes1, order1, u);
  std::cout << "stress element 1" << std::endl;
  std::cout << S1 << std::endl;
  std::cout << std::endl;

  VectorXd S2 = element_stress(E, v, t, nodes2, order2, u);
  std::cout << "stress element 2" << std::endl;
  std::cout << S2 << std::endl;
  std::cout << std::endl;

  return u[0];
}
