#include "fem_solve.h"
#include "fem_element.h"
#include "householder.h"
#include <Eigen/Dense>
#include <cassert>
#include <chrono>
#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <unordered_set>

auto mdof = [](size_t node, size_t comp) -> int {
  return static_cast<int>(node * 2 + comp);
};

// ───────────────────────────────────────────────────────────────
//  Helper: validate that a matrix/vector has no NaN or Inf
// ───────────────────────────────────────────────────────────────
template <typename Derived>
static void check_finite(const Eigen::MatrixBase<Derived> &M,
                         const std::string &label) {
  if (!M.allFinite()) {
    std::ostringstream oss;
    oss << "ERROR: " << label << " contains NaN or Inf!\n" << M << "\n";
    throw std::runtime_error(oss.str());
  }
}

// ───────────────────────────────────────────────────────────────
//  Helper: validate SolverInput before any solve
// ───────────────────────────────────────────────────────────────
static void validate_input(const SolverInput &input) {
  if (input.nodes.empty())
    throw std::runtime_error("ERROR: No nodes in input.");

  if (input.elements.empty())
    throw std::runtime_error("ERROR: No elements in input.");

  const size_t nnodes = input.nodes.size();

  // Check element connectivity bounds
  for (size_t ei = 0; ei < input.elements.size(); ++ei) {
    const auto &e = input.elements[ei];
    for (size_t i = 0; i < 6; ++i) {
      if (e[i] >= nnodes) {
        std::ostringstream oss;
        oss << "ERROR: Element " << ei << " local node " << i
            << " has global index " << e[i] << " which is out of range [0, "
            << nnodes - 1 << "].";
        throw std::runtime_error(oss.str());
      }
    }
  }

  // Check boundary edge node bounds
  for (size_t bi = 0; bi < input.boundary_edges.size(); ++bi) {
    const auto &edge = input.boundary_edges[bi];
    for (size_t nid : {edge.n1, edge.n2, edge.n3}) {
      if (nid >= nnodes) {
        std::ostringstream oss;
        oss << "ERROR: Boundary edge " << bi << " references node " << nid
            << " which is out of range [0, " << nnodes - 1 << "].";
        throw std::runtime_error(oss.str());
      }
    }
    if (edge.h < 0.0) {
      std::ostringstream oss;
      oss << "ERROR: Boundary edge " << bi
          << " has negative convection coefficient h = " << edge.h;
      throw std::runtime_error(oss.str());
    }
  }

  // Check boundary edge node bounds
  for (size_t bi = 0; bi < input.traction_edges.size(); ++bi) {
    const auto &edge = input.traction_edges[bi];
    for (size_t nid : {edge.n1, edge.n2, edge.n3}) {
      if (nid >= nnodes) {
        std::ostringstream oss;
        oss << "ERROR: Boundary edge " << bi << " references node " << nid
            << " which is out of range [0, " << nnodes - 1 << "].";
        throw std::runtime_error(oss.str());
      }
    }
  }

  // Check fixed DOF bounds
  const size_t mech_ndof = input.nodes.size() * 2;
  if (input.fixed_dofs.size() != input.fixed_values.size()) {
    throw std::runtime_error(
        "ERROR: fixed_dofs and fixed_values have different sizes.");
  }
  for (size_t i = 0; i < input.fixed_dofs.size(); ++i) {
    if (input.fixed_dofs[i] >= mech_ndof) {
      std::ostringstream oss;
      oss << "ERROR: fixed_dofs[" << i << "] = " << input.fixed_dofs[i]
          << " is out of range [0, " << mech_ndof - 1 << "].";
      throw std::runtime_error(oss.str());
    }
  }

  // Check force vector size
  if (input.forces.size() != mech_ndof) {
    std::ostringstream oss;
    oss << "ERROR: Force vector size (" << input.forces.size()
        << ") != 2 * number of nodes (" << mech_ndof << ").";
    throw std::runtime_error(oss.str());
  }

  // Check node coordinates are finite
  for (size_t i = 0; i < input.nodes.size(); ++i) {
    if (!std::isfinite(input.nodes[i].x) || !std::isfinite(input.nodes[i].y)) {
      std::ostringstream oss;
      oss << "ERROR: Node " << i << " has non-finite coordinates ("
          << input.nodes[i].x << ", " << input.nodes[i].y << ").";
      throw std::runtime_error(oss.str());
    }
  }
}

// ═══════════════════════════════════════════════════════════════
//  THERMAL SOLVER
// ═══════════════════════════════════════════════════════════════
Eigen::VectorXd solve_temperature(const SolverInput &input, double k) {

  if (k <= 0.0)
    throw std::runtime_error("ERROR: Thermal conductivity k must be > 0.");

  if (input.boundary_edges.empty())
    std::cerr << "WARNING: No boundary edges — thermal system may be "
                 "singular (no BCs).\n";

  size_t ndof = input.nodes.size();
  Eigen::SparseMatrix<double> Ktt(eidx(ndof), eidx(ndof));
  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(input.elements.size() * 36 +
                   input.boundary_edges.size() * 9);
  Eigen::VectorXd Qt = Eigen::VectorXd::Zero(eidx(ndof));

  const double s_gp[2] = {-1.0 / std::sqrt(3.0), 1.0 / std::sqrt(3.0)};
  const double w_gp[2] = {1.0, 1.0};

  // ── Convection boundary edges ──
  for (size_t bi = 0; bi < input.boundary_edges.size(); ++bi) {
    const auto &edge = input.boundary_edges[bi];
    size_t nodes[3] = {edge.n1, edge.n2, edge.n3};

    Eigen::Matrix3d Ke = Eigen::Matrix3d::Zero();
    Eigen::Vector3d Fe = Eigen::Vector3d::Zero();

    const Node &x1 = input.nodes[nodes[0]];
    const Node &x2 = input.nodes[nodes[2]];

    double dx = x2.x - x1.x;
    double dy = x2.y - x1.y;
    double L = std::sqrt(dx * dx + dy * dy);

    if (L < 1e-15) {
      std::ostringstream oss;
      oss << "ERROR: Boundary edge " << bi << " has zero length (nodes "
          << nodes[0] << " and " << nodes[2] << " coincide).";
      throw std::runtime_error(oss.str());
    }

    double detJ = L / 2.0;

    for (int gp = 0; gp < 2; ++gp) {
      double s = s_gp[gp];
      double w = w_gp[gp];

      double N[3];
      N[0] = 0.5 * s * (s - 1.0);
      N[1] = 1.0 - s * s;
      N[2] = 0.5 * s * (s + 1.0);

      for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
          Ke(eidx(i), eidx(j)) += edge.h * N[i] * N[j] * detJ * w;
        }
        Fe(eidx(i)) += edge.h * edge.Tinf * N[i] * detJ * w;
      }
    }

    check_finite(Ke, "Boundary edge " + std::to_string(bi) + " Ke");
    check_finite(Fe, "Boundary edge " + std::to_string(bi) + " Fe");

    for (size_t i = 0; i < 3; ++i) {
      size_t gi = nodes[i];
      for (size_t j = 0; j < 3; ++j) {
        size_t gj = nodes[j];
        triplets.emplace_back(gi, gj, Ke(eidx(i), eidx(j)));
      }
      Qt(eidx(gi)) += Fe(eidx(i));
    }
  }

  // ── Conduction element stiffness ──
  for (size_t ei = 0; ei < input.elements.size(); ++ei) {
    const auto &e = input.elements[ei];
    auto Ke = element_stiffness_heat(k, 1.0, input.nodes, e);

    if (!Ke.allFinite()) {
      std::ostringstream oss;
      oss << "ERROR: Thermal element " << ei << " stiffness contains NaN/Inf.\n"
          << "  Nodes: ";
      for (size_t i = 0; i < 6; ++i)
        oss << e[i] << "(" << input.nodes[e[i]].x << "," << input.nodes[e[i]].y
            << ") ";
      oss << "\n  Check for degenerate/inverted element geometry.";
      throw std::runtime_error(oss.str());
    }

    for (size_t i = 0; i < 6; ++i)
      for (size_t j = 0; j < 6; ++j)
        triplets.emplace_back(e[i], e[j], Ke(eidx(i), eidx(j)));
  }

  // ── Assemble and solve ──
  Ktt.setFromTriplets(triplets.begin(), triplets.end());

  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
  solver.compute(Ktt);

  if (solver.info() != Eigen::Success) {
    throw std::runtime_error(
        "ERROR: LDLT factorization failed for thermal system. "
        "Matrix may be singular or not positive definite.");
  }

  // Check pivots
  Eigen::VectorXd D = solver.vectorD();
  for (Eigen::Index i = 0; i < D.size(); ++i) {
    if (!std::isfinite(D(i))) {
      std::ostringstream oss;
      oss << "ERROR: Thermal LDLT pivot D(" << i << ") = " << D(i)
          << " is not finite.";
      throw std::runtime_error(oss.str());
    }
    if (D(i) <= 0.0) {
      std::ostringstream oss;
      oss << "ERROR: Thermal LDLT pivot D(" << i << ") = " << D(i)
          << " is non-positive. Matrix is not SPD.\n"
          << "  Likely cause: inverted element or missing boundary conditions.";
      throw std::runtime_error(oss.str());
    }
  }

  Eigen::VectorXd T = solver.solve(Qt);

  if (solver.info() != Eigen::Success)
    throw std::runtime_error("ERROR: Thermal back-substitution failed.");

  if (!T.allFinite())
    throw std::runtime_error(
        "ERROR: Temperature solution contains NaN or Inf.");

  return T;
}

// ═══════════════════════════════════════════════════════════════
//  DISPLACEMENT SOLVER
// ═══════════════════════════════════════════════════════════════
Eigen::VectorXd solve_displacement(const SolverInput &input,
                                   const Eigen::VectorXd &T, double E,
                                   double nu, double t, double alpha,
                                   double T0) {

  if (E <= 0.0)
    throw std::runtime_error("ERROR: Young's modulus E must be > 0.");
  if (nu <= -1.0 || nu >= 0.5)
    throw std::runtime_error("ERROR: Poisson's ratio nu must be in (-1, 0.5).");
  if (t <= 0.0)
    throw std::runtime_error("ERROR: Thickness t must be > 0.");

  if (T.size() != static_cast<int>(input.nodes.size())) {
    std::ostringstream oss;
    oss << "ERROR: Temperature vector size (" << T.size()
        << ") != number of nodes (" << input.nodes.size() << ").";
    throw std::runtime_error(oss.str());
  }

  size_t ndof = input.nodes.size() * 2;
  Eigen::SparseMatrix<double> K(eidx(ndof), eidx(ndof));
  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(input.elements.size() * 144);

  Eigen::VectorXd F =
      Eigen::Map<const Eigen::VectorXd>(input.forces.data(), eidx(ndof));

  // ── Mechanical stiffness ──
  for (size_t ei = 0; ei < input.elements.size(); ++ei) {
    const auto &e = input.elements[ei];
    auto Ke = element_stiffness_lst(E, nu, t, input.nodes, e);

    if (!Ke.allFinite()) {
      std::ostringstream oss;
      oss << "ERROR: Mechanical element " << ei
          << " stiffness contains NaN/Inf.\n"
          << "  Nodes: ";
      for (size_t i = 0; i < 6; ++i)
        oss << e[i] << "(" << input.nodes[e[i]].x << "," << input.nodes[e[i]].y
            << ") ";
      throw std::runtime_error(oss.str());
    }

    for (size_t i = 0; i < 6; i++) {
      for (size_t j = 0; j < 6; j++) {
        for (size_t a = 0; a < 2; a++) {
          for (size_t b = 0; b < 2; b++) {
            int row = mdof(e[i], a);
            int col = mdof(e[j], b);
            triplets.emplace_back(row, col,
                                  Ke(eidx(2 * i + a), eidx(2 * j + b)));
          }
        }
      }
    }
  }

  K.setFromTriplets(triplets.begin(), triplets.end());

  assert(K.rows() == K.cols());
  assert(static_cast<size_t>(K.rows()) == ndof);

  // ── Thermal force ──
  for (size_t ei = 0; ei < input.elements.size(); ++ei) {
    const auto &e = input.elements[ei];
    Eigen::Matrix<double, 6, 1> Te;
    for (size_t i = 0; i < 6; ++i)
      Te(eidx(i)) = T(eidx(e[i]));

    auto Fe_th = element_thermal_force(E, nu, alpha, T0, input.nodes, e, Te);

    if (!Fe_th.allFinite()) {
      std::ostringstream oss;
      oss << "ERROR: Element " << ei << " thermal force contains NaN/Inf.";
      throw std::runtime_error(oss.str());
    }

    for (size_t i = 0; i < 6; ++i) {
      for (size_t a = 0; a < 2; ++a) {
        int global_row = mdof(e[i], a);
        Eigen::Index local_index = eidx(2 * i + a);
        F(global_row) += Fe_th(local_index);
      }
    }
  }

  const double s_gp[2] = {-1.0 / std::sqrt(3.0), 1.0 / std::sqrt(3.0)};
  const double w_gp[2] = {1.0, 1.0};

  // ── Convection boundary edges ──
  for (size_t bi = 0; bi < input.traction_edges.size(); ++bi) {
    const auto &edge = input.traction_edges[bi];
    size_t nodes[3] = {edge.n1, edge.n2, edge.n3};

    // Edge length (straight edge: n1 → n3)
    const Node &x1 = input.nodes[nodes[0]];
    const Node &x2 = input.nodes[nodes[2]];

    double dx = x2.x - x1.x;
    double dy = x2.y - x1.y;
    double L = std::sqrt(dx * dx + dy * dy);

    if (L < 1e-15) {
      std::ostringstream oss;
      oss << "ERROR: Traction edge " << bi << " has zero length (nodes "
          << nodes[0] << " and " << nodes[2] << " coincide).";
      throw std::runtime_error(oss.str());
    }

    double detJ = L / 2.0;

    for (int gp = 0; gp < 2; ++gp) {
      double s = s_gp[gp];
      double w = w_gp[gp];

      // 1D quadratic shape functions on [-1, +1]
      // N1 at s=-1 (node n1), N3 at s=+1 (node n3), N2 at s=0 (mid-node n2)

      double N[3];
      N[0] = 0.5 * s * (s - 1.0);
      N[1] = 1.0 - s * s;
      N[2] = 0.5 * s * (s + 1.0);

      for (size_t i = 0; i < 3; ++i) {
        Eigen::Index dof_x = eidx(nodes[i] * 2);
        Eigen::Index dof_y = eidx(nodes[i] * 2 + 1);
        F(dof_x) += N[i] * edge.tx * detJ * w;
        F(dof_y) += N[i] * edge.ty * detJ * w;
      }
    }
  }

  // ── Dirichlet BCs ──
  if (input.fixed_dofs.empty())
    std::cerr << "WARNING: No Dirichlet BCs — mechanical system will be "
                 "singular (rigid body modes).\n";

  for (size_t i = 0; i < input.fixed_dofs.size(); i++) {
    size_t d = input.fixed_dofs[i];
    F -= K.col(eidx(d)) * input.fixed_values[i];
  }

  std::unordered_set<size_t> fixed(input.fixed_dofs.begin(),
                                   input.fixed_dofs.end());

  std::vector<Eigen::Index> map(ndof, -1);
  int counter_free = 0;

  for (size_t i = 0; i < ndof; i++) {
    if (!fixed.count(i))
      map[i] = counter_free++;
  }

  if (counter_free == 0)
    throw std::runtime_error("ERROR: All DOFs are fixed — nothing to solve.");

  // ── Extract free-free system ──
  std::vector<Eigen::Triplet<double>> triplets_ff;

  for (int k = 0; k < K.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(K, k); it; ++it) {
      auto i = static_cast<size_t>(it.row());
      auto j = static_cast<size_t>(it.col());
      if (map[i] != -1 && map[j] != -1)
        triplets_ff.emplace_back(map[i], map[j], it.value());
    }
  }

  Eigen::VectorXd F_ff(counter_free);
  for (size_t i = 0; i < ndof; i++) {
    if (map[i] != -1)
      F_ff(map[i]) = F(eidx(i));
  }

  check_finite(F_ff, "Reduced force vector F_ff");

  Eigen::SparseMatrix<double> K_ff(counter_free, counter_free);
  K_ff.setFromTriplets(triplets_ff.begin(), triplets_ff.end());

  // ── Solve ──
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
  solver.compute(K_ff);

  if (solver.info() != Eigen::Success) {
    throw std::runtime_error(
        "ERROR: LDLT factorization failed for displacement system.");
  }

  // Check pivots
  Eigen::VectorXd D = solver.vectorD();
  for (int i = 0; i < D.size(); ++i) {
    if (!std::isfinite(D(i))) {
      std::ostringstream oss;
      oss << "ERROR: Displacement LDLT pivot D(" << i << ") = " << D(i)
          << " is not finite.";
      throw std::runtime_error(oss.str());
    }
    if (D(i) <= 0.0) {
      std::ostringstream oss;
      oss << "ERROR: Displacement LDLT pivot D(" << i << ") = " << D(i)
          << " is non-positive.\n"
          << "  Likely cause: insufficient boundary conditions "
             "(rigid body mode) or inverted element.";
      throw std::runtime_error(oss.str());
    }
  }

  Eigen::VectorXd u_free = solver.solve(F_ff);

  if (solver.info() != Eigen::Success)
    throw std::runtime_error("ERROR: Displacement back-substitution failed.");

  if (!u_free.allFinite())
    throw std::runtime_error(
        "ERROR: Displacement solution contains NaN or Inf.");

  // ── Scatter back to full vector ──
  Eigen::VectorXd u = Eigen::VectorXd::Zero(eidx(ndof));
  for (size_t i = 0; i < ndof; i++) {
    if (map[i] != -1)
      u(eidx(i)) = u_free(map[i]);
  }
  for (size_t i = 0; i < input.fixed_dofs.size(); i++) {
    u(eidx(input.fixed_dofs[i])) = input.fixed_values[i];
  }

  return u;
}

// ═══════════════════════════════════════════════════════════════
//  STRESS RECOVERY
// ═══════════════════════════════════════════════════════════════
std::vector<Eigen::Vector3d> solve_nodal_stress(const SolverInput &input,
                                                const Eigen::VectorXd &u,
                                                const Eigen::VectorXd &T,
                                                double E, double nu,
                                                double alpha, double T0) {

  if (u.size() != static_cast<int>(input.nodes.size() * 2)) {
    std::ostringstream oss;
    oss << "ERROR: Displacement vector size (" << u.size()
        << ") != 2 * number of nodes (" << input.nodes.size() * 2 << ").";
    throw std::runtime_error(oss.str());
  }

  std::vector<Eigen::Vector3d> nodal_stress(input.nodes.size(),
                                            Eigen::Vector3d::Zero());
  std::vector<int> counter(input.nodes.size(), 0);

  for (size_t ei = 0; ei < input.elements.size(); ++ei) {
    const auto &e = input.elements[ei];

    auto sg = compute_gauss_stress_lst(E, nu, input.nodes, e, u, alpha, T0, &T);

    auto sn = extrapolate_to_nodes(sg);

    for (size_t i = 0; i < 6; ++i) {
      size_t gid = e[i];

      if (!std::isfinite(sn[i](0)) || !std::isfinite(sn[i](1)) ||
          !std::isfinite(sn[i](2))) {
        std::ostringstream oss;
        oss << "ERROR: Element " << ei << " node " << i << " (global " << gid
            << ") has non-finite stress: " << sn[i].transpose();
        throw std::runtime_error(oss.str());
      }

      nodal_stress[gid] += sn[i];
      counter[gid]++;
    }
  }

  for (size_t i = 0; i < nodal_stress.size(); ++i) {
    if (counter[i] > 0)
      nodal_stress[i] /= counter[i];
    else
      std::cerr << "WARNING: Node " << i
                << " is not referenced by any element (orphan node).\n";
  }

  return nodal_stress;
}

// ═══════════════════════════════════════════════════════════════
//  TOP-LEVEL DRIVER
// ═══════════════════════════════════════════════════════════════
SolverOutput solve(const SolverInput &input) {

  // ── Validate input first ──
  validate_input(input);

  const auto &mat = input.material;

  Eigen::VectorXd T;
  if ((input.boundary_edges.empty()) || (input.boundary_edges.size() == 0)) {
    T = Eigen::VectorXd::Constant(eidx(input.nodes.size()), mat.T0);
  } else {
    T = solve_temperature(input, mat.k);
  }
  Eigen::VectorXd u =
      solve_displacement(input, T, mat.E, mat.nu, mat.t, mat.alpha, mat.T0);
  std::vector<Eigen::Vector3d> nodal_stress =
      solve_nodal_stress(input, u, T, mat.E, mat.nu, mat.alpha, mat.T0);

  SolverOutput out;
  out.stress.reserve(nodal_stress.size());
  for (const auto &s : nodal_stress)
    out.stress.push_back({s(0), s(1), s(2)});

  for (size_t i = 0; i < input.nodes.size(); i++) {
    out.disp.push_back({u(eidx(2 * i)), u(eidx(2 * i + 1))});
    out.temperature.push_back(T(eidx(i)));
  }

  return out;
}