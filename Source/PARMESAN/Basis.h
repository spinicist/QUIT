#pragma once

#include <Eigen/Dense>
#include <string>

auto ReadBasis(std::string const &path) -> Eigen::MatrixXd;
