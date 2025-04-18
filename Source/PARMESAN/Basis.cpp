#include "Basis.h"

#include "JSON.h"

auto ReadBasis(std::string const &path) -> Eigen::MatrixXd {
    if (path.size()) {
        auto bj = QI::ReadJSON(path);
        return QI::MatrixFromJSON(bj, "basis", 1.0, -1, -1);
    } else {
        return Eigen::MatrixXd();
    }
}
