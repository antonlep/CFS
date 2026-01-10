// #include <Eigen/Dense>

struct Point {
  float x;
  float y;
};

auto stiffness_matrix(Point i, Point j, Point k) { return i.x; }

int add(int a, int b) { return a + b; }

float plate(float h, float w, float t, float f) {
  float e = 210E9;
  float v = 0.3;
  float F = 0.5 * t * h;
  Point p = {h, w};

  return stiffness_matrix(p, p, p);
}

// Eigen::VectorXd solve_poisson(const Mesh &mesh) {
//   // assemble matrix A and vector b
//   // solve Ax = b
//   return x;
// }
