#pragma once

#include "Macro.h"
#include <Eigen/Dense>
#include <functional>
#include <unsupported/Eigen/MatrixFunctions>

template <typename AugmentedMatrix>
auto SolveSteadyState(AugmentedMatrix const &X)
    -> Eigen::Vector<typename AugmentedMatrix::Scalar, AugmentedMatrix::RowsAtCompileTime> {
    // Solve for steady-state and re-augment
    using T               = typename AugmentedMatrix::Scalar;
    const long N          = AugmentedMatrix::RowsAtCompileTime;
    using ReducedMatrix   = Eigen::Matrix<T, N - 1, N - 1>;
    using ReducedVector   = Eigen::Vector<T, N - 1>;
    using AugmentedVector = Eigen::Vector<T, N>;

    ReducedMatrix   Xr   = (X - AugmentedMatrix::Identity()).template topLeftCorner<N - 1, N - 1>();
    ReducedVector   b    = -X.template topRightCorner<N - 1, 1>();
    ReducedVector   m_ss = Xr.partialPivLu().solve(b);
    AugmentedVector m_aug;
    m_aug << m_ss, 1.;
    return m_aug;
}

// Geometric Sum formula to get average
template <typename AugmentedMatrix, typename AugmentedVector>
auto GeometricAvg(AugmentedMatrix const &X,
                  AugmentedMatrix const &Xn,
                  AugmentedVector const &a,
                  int const &            n)
    -> Eigen::Vector<typename AugmentedMatrix::Scalar, AugmentedMatrix::RowsAtCompileTime - 1> {
    using T      = typename AugmentedMatrix::Scalar;
    const long N = AugmentedMatrix::RowsAtCompileTime;
    // using ReducedMatrix = Eigen::Matrix<T, N - 1, N - 1>;
    using ReducedVector = Eigen::Vector<T, N - 1>;

    AugmentedMatrix const LHS = (AugmentedMatrix::Identity() - X);
    ReducedVector const   RHS = ((AugmentedMatrix::Identity() - Xn) * a).template head<N - 1>() -
                              (n * LHS.template topRightCorner<N - 1, 1>());
    ReducedVector const m_gm =
        LHS.template topLeftCorner<N - 1, N - 1>().partialPivLu().solve(RHS) / n;
    return m_gm;
}
