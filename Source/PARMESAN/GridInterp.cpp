#include "GridInterp.h"
#include "Log.h"

namespace QI {

InterpGrid::InterpGrid(Eigen::ArrayXd const  &x,
                         Eigen::ArrayXd const  &y,
                         Eigen::ArrayXXd const &z) :
    m_x(x),
    m_y(y), m_z(z), nX(m_x.size()), nY(m_y.size()) {
    if (nX != m_z.cols()) {
        Fail("Grid had {} columns expected {}", m_z.cols(), nX);
    }
    if (nY != m_z.rows()) {
        Fail("Grid had {} rows expected {}", m_z.rows(), nY);
    }
}

auto InterpGrid::FindIndices(Eigen::ArrayXd const &a, double const x) -> InterpGrid::IndexPair {
    int i0 = 0, i1 = 0;
    for (; i1 < a.size(); i1++) {
        if (x < a[i1])
            break;
    }
    if (i1 == 0) {
        i1 = 1;
    }
    if (i1 == a.size()) {
        i1 = a.size() - 1;
        i0 = i1 - 1;
    }
    return {i0, i1};
}

double InterpGrid::operator()(const double &x, const double &y) const {
    auto const ix = FindIndices(m_x, x);
    auto const iy = FindIndices(m_y, y);

    double const dx = m_x[ix.hi] - m_x[ix.lo];
    double const lx = (x - m_x[ix.lo]) / dx;
    double const dy = m_y[iy.hi] - m_y[iy.lo];
    double const ly = (y - m_y[iy.lo]) / dy;

    Eigen::Matrix2d Q;
    Q << m_z(ix.lo, iy.lo), m_z(ix.lo, iy.hi), //
    m_z(ix.hi, iy.lo), m_z(ix.hi, iy.hi);

    Eigen::RowVector2d vx;
    vx << (1 - lx), lx;
    Eigen::Vector2d vy;
    vy << (1- ly), ly;
    
    return vx * Q * vy;
}

} // namespace QI