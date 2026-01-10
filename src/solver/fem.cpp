#include <Eigen/Dense>

struct Node {
  float x;
  float y;
};

float area(Node i, Node j, Node k) {
  return 0.5 * (i.x * (j.y - k.y) + j.x * (k.y - i.y) + k.x * (i.y - j.y));
}

auto stiffness_matrix(float E, float v, Node i, Node j, Node k) {
  float A = area(i, j, k);
  Eigen::MatrixXd B(3, 6);
  B << j.y - k.y, 0, k.y - i.y, 0, i.y - j.y, 0, 0, k.x - j.x, 0, i.x - k.x, 0,
      j.x - i.x, k.x - j.x, j.y - k.y, i.x - k.x, k.y - i.y, j.x - i.x,
      i.y - j.y;
  B = 1 / (2 * A) * B;
  // float B[3][6] = 1/(2*A)*{
  //     {j.y - k.y, 0, k.y - i.y, 0, i.y - j.y, 0},
  //     {0, k.x - j.x, 0, i.x - k.x, 0, j.x - i.x},
  //     {k.x - j.x, j.y - k.y, i.x - k.x, k.y - i.y, j.x - i.x, i.y - j.y}
  // };
  // float D[3][3] = E/(1-v**2)*{{1, v, 0}, {v, 1, 0}, {0, 0, (1-v)/2}};
  return i.x;
}

int add(int a, int b) { return a + b; }

float plate(float h, float w, float t, float f) {
  float E = 210E9;
  float v = 0.3;
  // float F = 0.5 * t * h;
  Node p1 = {0, 0};
  Node p2 = {0, 10};
  Node p3 = {20, 10};
  Node p4 = {0, 20};

  return stiffness_matrix(E, v, p1, p2, p3);
}

// Eigen::VectorXd solve_poisson(const Mesh &mesh) {
//   // assemble matrix A and vector b
//   // solve Ax = b
//   return x;
// }
