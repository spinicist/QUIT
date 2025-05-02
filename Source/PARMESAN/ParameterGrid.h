#pragma once
#include <Eigen/Dense>

namespace QI {

auto RegularPars(Eigen::ArrayXd const &lo, Eigen::ArrayXd const &hi, Eigen::ArrayXi const &N)
    -> Eigen::ArrayXXd;

auto RandomPars(Eigen::ArrayXd const &lo, Eigen::ArrayXd const &hi, Eigen::Index const N)
    -> Eigen::ArrayXXd;

} // namespace QI
