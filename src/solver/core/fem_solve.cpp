#include "fem_solve.h"
#include "fem_element.h"
#include "householder.h"
#include <Eigen/Dense>
#include <cassert>
#include <iostream>
#include <unordered_set>

SolverOutput fem_solve(const SolverInput &input) {
  constexpr double E = 210e6;
  constexpr double nu = 0.3;
  constexpr double t = 1.0;

  size_t ndof = input.nodes.size() * 2;

  Eigen::MatrixXd K = Eigen::MatrixXd::Zero(ndof, ndof);

  for (const auto &e : input.elements) {
    auto ke = element_stiffness_lst(E, nu, t, input.nodes, e);
    for (size_t i = 0; i < 6; i++)
      for (size_t j = 0; j < 6; j++)
        K.block<2, 2>(2 * e[i], 2 * e[j]) += ke.block<2, 2>(2 * i, 2 * j);
  }

  Eigen::VectorXd F =
      Eigen::Map<const Eigen::VectorXd>(input.forces.data(), ndof);

  assert(K.rows() == K.cols());
  assert(K.rows() == ndof);

  if (!K.allFinite())
    throw std::runtime_error("Stiffness matrix contains NaN/Inf");

  if (!F.allFinite())
    throw std::runtime_error("Load vector contains NaN/Inf");

  // Apply Dirichlet BC
  for (size_t i = 0; i < input.fixed_dofs.size(); i++) {
    size_t d = input.fixed_dofs[i];
    F -= K.col(d) * input.fixed_values[i];
  }

  std::unordered_set<size_t> fixed(input.fixed_dofs.begin(),
                                   input.fixed_dofs.end());

  std::vector<size_t> free;
  for (size_t i = 0; i < ndof; i++)
    if (!fixed.count(i))
      free.push_back(i);

  Eigen::VectorXi free_idx(free.size());
  for (size_t i = 0; i < free.size(); i++)
    free_idx(i) = free[i];

  Eigen::VectorXd u_free = householder(K(free_idx, free_idx), F(free_idx));

  Eigen::VectorXd u = Eigen::VectorXd::Zero(ndof);
  for (size_t i = 0; i < free.size(); i++)
    u(free[i]) = u_free(i);
  for (size_t i = 0; i < input.fixed_dofs.size(); i++)
    u(input.fixed_dofs[i]) = input.fixed_values[i];

  if (!u.allFinite()) {
    throw std::runtime_error("Solution vector contains NaN/Inf");
  }

  SolverOutput out;

  // for (const auto &e : input.elements) {
  //   auto stresses = element_stress_lst(E, nu, input.nodes, e, u);
  //   for (const auto &s : stresses) {
  //     out.stress.push_back(s);
  //   }
  // }

  std::vector<Eigen::Vector3d> nodal_stress(input.nodes.size(),
                                            Eigen::Vector3d::Zero());
  std::vector<int> counter(input.nodes.size(), 0);

  for (const auto &e : input.elements) {

    auto sg = compute_gauss_stress_lst(E, nu, input.nodes, e, u);
    auto sn = extrapolate_to_nodes(sg);
    for (int i = 0; i < 3; ++i)
      std::cout << sg[i].transpose() << std::endl;

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
