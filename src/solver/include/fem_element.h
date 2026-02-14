#pragma once
#include "types.h"
#include <Eigen/Dense>

Eigen::Matrix<double, 6, 6> element_stiffness_cst(double E, double nu, double t,
                                                  const Nodes &nodes,
                                                  const Element &e);

Eigen::Vector3d element_stress_cst(double E, double nu, const Nodes &nodes,
                                   const Element &e, const Eigen::VectorXd &u);

Eigen::Matrix<double, 12, 12> element_stiffness_lst(double E, double nu,
                                                    double t,
                                                    const Nodes &nodes,
                                                    const Element &e);

std::array<Eigen::Vector3d, 3>
compute_gauss_stress_lst(double E, double nu, const Nodes &nodes,
                         const Element &e, const Eigen::VectorXd &u);

std::array<Eigen::Vector3d, 6>
extrapolate_to_nodes(const std::array<Eigen::Vector3d, 3> &sg);