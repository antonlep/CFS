#include "fem_solve.h"
#include "fem_element.h"
#include "householder.h"
#include <Eigen/Dense>
#include <cassert>
#include <chrono>
#include <iostream>
#include <unordered_set>

SolverOutput fem_solve(const SolverInput &input) {

  constexpr double E = 210e6;
  constexpr double nu = 0.3;
  constexpr double t = 1.0;

  size_t ndof = input.nodes.size() * 2;

  // Eigen::MatrixXd K = Eigen::MatrixXd::Zero(ndof, ndof);
  Eigen::SparseMatrix<double> K(ndof, ndof);
  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(input.elements.size() * 144);

  for (const auto &e : input.elements) {
    auto ke = element_stiffness_lst(E, nu, t, input.nodes, e);
    for (size_t i = 0; i < 6; i++) {
      for (size_t j = 0; j < 6; j++) {
        int gi = e[i];
        int gj = e[j];
        int row = 2 * gi;
        int col = 2 * gj;
        triplets.emplace_back(row, col, ke(2 * i, 2 * j));
        triplets.emplace_back(row, col + 1, ke(2 * i, 2 * j + 1));
        triplets.emplace_back(row + 1, col, ke(2 * i + 1, 2 * j));
        triplets.emplace_back(row + 1, col + 1, ke(2 * i + 1, 2 * j + 1));
      }
    }
  }
  K.setFromTriplets(triplets.begin(), triplets.end());

  Eigen::VectorXd F =
      Eigen::Map<const Eigen::VectorXd>(input.forces.data(), ndof);

  assert(K.rows() == K.cols());
  assert(K.rows() == ndof);

  // if (!K.allFinite())
  //   throw std::runtime_error("Stiffness matrix contains NaN/Inf");

  // if (!F.allFinite())
  //   throw std::runtime_error("Load vector contains NaN/Inf");

  // Apply Dirichlet BC
  for (size_t i = 0; i < input.fixed_dofs.size(); i++) {
    size_t d = input.fixed_dofs[i];
    F -= K.col(d) * input.fixed_values[i];
  }

  std::unordered_set<size_t> fixed(input.fixed_dofs.begin(),
                                   input.fixed_dofs.end());

  std::vector<int> map(ndof, -1);
  int counter_free = 0;

  for (size_t i = 0; i < ndof; i++) {
    if (!fixed.count(i)) {
      map[i] = counter_free++;
    }
  }

  std::vector<Eigen::Triplet<double>> triplets_ff;

  for (int k = 0; k < K.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(K, k); it; ++it) {
      int i = it.row();
      int j = it.col();

      if (map[i] != -1 && map[j] != -1) {
        triplets_ff.emplace_back(map[i], map[j], it.value());
      }
    }
  }

  Eigen::VectorXd F_ff(counter_free);

  for (size_t i = 0; i < ndof; i++) {
    if (map[i] != -1) {
      F_ff(map[i]) = F(i);
    }
  }

  Eigen::SparseMatrix<double> K_ff(counter_free, counter_free);
  K_ff.setFromTriplets(triplets_ff.begin(), triplets_ff.end());

  std::vector<size_t> free;
  for (size_t i = 0; i < ndof; i++)
    if (!fixed.count(i))
      free.push_back(i);

  Eigen::VectorXi free_idx(free.size());
  for (size_t i = 0; i < free.size(); i++)
    free_idx(i) = free[i];

  // Eigen::VectorXd u_free = householder(K(free_idx, free_idx), F(free_idx));
  // Eigen::VectorXd u_free = K(free_idx, free_idx).llt().solve(F(free_idx));
  Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;
  solver.compute(K_ff);
  Eigen::VectorXd u_free = solver.solve(F_ff);

  Eigen::VectorXd u = Eigen::VectorXd::Zero(ndof);

  for (size_t i = 0; i < ndof; i++) {
    if (map[i] != -1) {
      u(i) = u_free(map[i]);
    }
  }

  for (size_t i = 0; i < input.fixed_dofs.size(); i++) {
    u(input.fixed_dofs[i]) = input.fixed_values[i];
  }

  SolverOutput out;

  std::vector<Eigen::Vector3d> nodal_stress(input.nodes.size(),
                                            Eigen::Vector3d::Zero());
  std::vector<int> counter(input.nodes.size(), 0);

  for (const auto &e : input.elements) {

    auto sg = compute_gauss_stress_lst(E, nu, input.nodes, e, u);
    auto sn = extrapolate_to_nodes(sg);

    for (int i = 0; i < 6; ++i) {
      int gid = e[i];
      nodal_stress[gid] += sn[i];
      counter[gid]++;
    }
  }

  for (size_t i = 0; i < nodal_stress.size(); ++i)
    if (counter[i] > 0)
      nodal_stress[i] /= counter[i];

  out.stress.reserve(nodal_stress.size());

  for (const auto &s : nodal_stress) {
    out.stress.push_back({s(0), s(1), s(2)});
  }

  for (size_t i = 0; i < input.nodes.size(); i++)
    out.disp.push_back({u(2 * i), u(2 * i + 1)});

  return out;
}
