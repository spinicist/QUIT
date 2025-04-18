#pragma once
#include <Eigen/Dense>

auto ParameterGrid(int const nPar,
    Eigen::ArrayXd const &lo,
    Eigen::ArrayXd const &hi,
    Eigen::ArrayXi const &N) -> Eigen::ArrayXXd;
