#include "fem_element.h"
#include <Eigen/Dense>

static double area(const Node &a, const Node &b, const Node &c) {
  return 0.5 * (a.x * (b.y - c.y) + b.x * (c.y - a.y) + c.x * (a.y - b.y));
}

Eigen::Matrix<double, 6, 6> element_stiffness(double E, double nu, double t,
                                              const Nodes &nodes,
                                              const Element &e) {
  const Node &i = nodes[e[0]];
  const Node &j = nodes[e[1]];
  const Node &k = nodes[e[2]];

  double A = area(i, j, k);

  Eigen::Matrix<double, 3, 6> B;
  B << j.y - k.y, 0, k.y - i.y, 0, i.y - j.y, 0, 0, k.x - j.x, 0, i.x - k.x, 0,
      j.x - i.x, k.x - j.x, j.y - k.y, i.x - k.x, k.y - i.y, j.x - i.x,
      i.y - j.y;

  B /= (2 * A);

  Eigen::Matrix3d D;
  D << 1, nu, 0, nu, 1, 0, 0, 0, (1 - nu) / 2;

  D *= E / (1 - nu * nu);

  return t * A * B.transpose() * D * B;
}

Eigen::Vector3d element_stress(double E, double nu, const Nodes &nodes,
                               const Element &e, const Eigen::VectorXd &u) {
  Eigen::Matrix<double, 6, 1> ue;
  for (size_t i = 0; i < 3; i++) {
    ue(2 * i) = u(2 * e[i]);
    ue(2 * i + 1) = u(2 * e[i] + 1);
  }

  double A = area(nodes[e[0]], nodes[e[1]], nodes[e[2]]);

  Eigen::Matrix<double, 3, 6> B;
  B << nodes[e[1]].y - nodes[e[2]].y, 0, nodes[e[2]].y - nodes[e[0]].y, 0,
      nodes[e[0]].y - nodes[e[1]].y, 0, 0, nodes[e[2]].x - nodes[e[1]].x, 0,
      nodes[e[0]].x - nodes[e[2]].x, 0, nodes[e[1]].x - nodes[e[0]].x,
      nodes[e[2]].x - nodes[e[1]].x, nodes[e[1]].y - nodes[e[2]].y,
      nodes[e[0]].x - nodes[e[2]].x, nodes[e[2]].y - nodes[e[0]].y,
      nodes[e[1]].x - nodes[e[0]].x, nodes[e[0]].y - nodes[e[1]].y;

  B /= (2 * A);

  Eigen::Matrix3d D;
  D << 1, nu, 0, nu, 1, 0, 0, 0, (1 - nu) / 2;

  D *= E / (1 - nu * nu);

  return D * B * ue;
}

static void compute_B_LST(const Nodes &nodes, const Element &e, double xi,
                          double eta, Eigen::Matrix<double, 3, 12> &B) {
  // Shape function derivatives in natural coordinates
  double L1 = 1.0 - xi - eta;
  double L2 = xi;
  double L3 = eta;

  double dN_dxi[6] = {-3 + 4 * L1,   4 * L2 - 1, 0,
                      4 * (L1 - L2), 4 * L3,     -4 * L3};

  double dN_deta[6] = {-3 + 4 * L1, 0,      4 * L3 - 1,
                       -4 * L2,     4 * L2, 4 * (L1 - L3)};

  // Jacobian
  Eigen::Matrix2d J = Eigen::Matrix2d::Zero();

  for (int i = 0; i < 6; ++i) {
    const Node &n = nodes[e[i]];
    J(0, 0) += dN_dxi[i] * n.x;
    J(0, 1) += dN_dxi[i] * n.y;
    J(1, 0) += dN_deta[i] * n.x;
    J(1, 1) += dN_deta[i] * n.y;
  }

  Eigen::Matrix2d invJ = J.inverse();

  B.setZero();

  for (int i = 0; i < 6; ++i) {

    Eigen::Vector2d dN_nat(dN_dxi[i], dN_deta[i]);
    Eigen::Vector2d dN = invJ * dN_nat;

    double dNx = dN(0);
    double dNy = dN(1);

    B(0, 2 * i) = dNx;
    B(1, 2 * i + 1) = dNy;
    B(2, 2 * i) = dNy;
    B(2, 2 * i + 1) = dNx;
  }
}

std::array<Eigen::Vector3d, 3>
compute_gauss_stress_lst(double E, double nu, const Nodes &nodes,
                         const Element &e, const Eigen::VectorXd &u) {
  Eigen::Matrix<double, 12, 1> ue;

  for (int i = 0; i < 6; ++i) {
    ue(2 * i) = u(2 * e[i]);
    ue(2 * i + 1) = u(2 * e[i] + 1);
  }

  Eigen::Matrix3d D;
  D << 1, nu, 0, nu, 1, 0, 0, 0, (1 - nu) / 2;

  D *= E / (1 - nu * nu);

  const double xi[3] = {1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0};
  const double eta[3] = {1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0};

  std::array<Eigen::Vector3d, 3> stresses;

  for (int gp = 0; gp < 3; ++gp) {

    Eigen::Matrix<double, 3, 12> B;
    compute_B_LST(nodes, e, xi[gp], eta[gp], B);

    stresses[gp] = D * B * ue;
  }

  return stresses;
}

std::array<Eigen::Vector3d, 6>
extrapolate_to_nodes(const std::array<Eigen::Vector3d, 3> &sg) {
  std::array<Eigen::Vector3d, 6> sn;

  Eigen::Matrix3d A;
  A << 1, 1.0 / 6.0, 1.0 / 6.0, 1, 2.0 / 3.0, 1.0 / 6.0, 1, 1.0 / 6.0,
      2.0 / 3.0;

  Eigen::Matrix3d Ainv = A.inverse();

  const double xi[6] = {0, 1, 0, 0.5, 0.5, 0};
  const double eta[6] = {0, 0, 1, 0, 0.5, 0.5};

  for (int comp = 0; comp < 3; ++comp) {

    Eigen::Vector3d s;
    for (int i = 0; i < 3; ++i)
      s(i) = sg[i](comp);

    Eigen::Vector3d coeff = Ainv * s;

    for (int i = 0; i < 6; ++i)
      sn[i](comp) = coeff(0) + coeff(1) * xi[i] + coeff(2) * eta[i];
  }

  return sn;
}

Eigen::Matrix<double, 12, 12> element_stiffness_lst(double E, double nu,
                                                    double t,
                                                    const Nodes &nodes,
                                                    const Element &e) {
  Eigen::Matrix<double, 12, 12> Ke = Eigen::Matrix<double, 12, 12>::Zero();

  // Elasticity matrix (plane stress)
  Eigen::Matrix3d D;
  D << 1, nu, 0, nu, 1, 0, 0, 0, (1 - nu) / 2;

  D *= E / (1 - nu * nu);

  // 3-point Gauss integration
  const double xi[3] = {1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0};
  const double eta[3] = {1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0};
  const double w = 1.0 / 6.0;

  for (int gp = 0; gp < 3; ++gp) {
    double L1 = 1.0 - xi[gp] - eta[gp];
    double L2 = xi[gp];
    double L3 = eta[gp];

    // Shape function derivatives in natural coordinates
    double dN_dxi[6] = {-(4 * L1 - 1), 4 * L2 - 1, 0.0,
                        4 * (L1 - L2), 4 * L3,     -4 * L3};

    double dN_deta[6] = {-(4 * L1 - 1), 0.0,    4 * L3 - 1,
                         -4 * L2,       4 * L2, 4 * (L1 - L3)};

    // Build Jacobian
    Eigen::Matrix2d J = Eigen::Matrix2d::Zero();

    for (int i = 0; i < 6; ++i) {
      const Node &n = nodes[e[i]];
      J(0, 0) += dN_dxi[i] * n.x;
      J(0, 1) += dN_dxi[i] * n.y;
      J(1, 0) += dN_deta[i] * n.x;
      J(1, 1) += dN_deta[i] * n.y;
    }

    double detJ = J.determinant();
    Eigen::Matrix2d invJ = J.inverse();

    // Build B matrix (3x12)
    Eigen::Matrix<double, 3, 12> B = Eigen::Matrix<double, 3, 12>::Zero();

    for (int i = 0; i < 6; ++i) {
      Eigen::Vector2d dN_nat(dN_dxi[i], dN_deta[i]);
      Eigen::Vector2d dN = invJ * dN_nat;

      double dNx = dN(0);
      double dNy = dN(1);

      B(0, 2 * i) = dNx;
      B(1, 2 * i + 1) = dNy;
      B(2, 2 * i) = dNy;
      B(2, 2 * i + 1) = dNx;
    }

    Ke += t * (B.transpose() * D * B) * detJ * w;
  }

  return Ke;
}