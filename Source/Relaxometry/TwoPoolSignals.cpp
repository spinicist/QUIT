#include "TwoPoolSignals.h"

// #define QI_DEBUG_BUILD
#include "Macro.h"

#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

#include "Helpers.h"

namespace QI {

Eigen::ArrayXcd SPGR2(double const            PD,
                      double const            T1_a,
                      double const            T2_a,
                      double const            T1_b,
                      double const            T2_b,
                      double const            tau_a,
                      double const            f_a,
                      double const            f0,
                      double const            B1,
                      SPGREchoSequence const &spgr) {
    Eigen::Matrix2d  A, eATR;
    Eigen::Vector2d  M0, Mz;
    Eigen::Vector2cd Mxy;
    double           k_ab, k_ba, f_b;
    CalcExchange(tau_a, f_a, f_b, k_ab, k_ba);
    M0 << f_a, f_b;
    A << -(1. / T1_a) - k_ab, k_ba, //
        k_ab, -(1. / T1_b) - k_ba;
    eATR                  = (spgr.TR * A).exp();
    Eigen::Matrix2cd echo = Eigen::Matrix2cd::Zero();
    // T2' absorbed into PD as it effects both components equally
    echo(0, 0) = std::polar(exp(-spgr.TE / T2_a), 2. * M_PI * f0 * spgr.TE);
    echo(1, 1) = std::polar(exp(-spgr.TE / T2_b), 2. * M_PI * f0 * spgr.TE);

    Eigen::Vector2d const RHS = (Eigen::Matrix2d::Identity() - eATR) * M0;
    Eigen::VectorXcd      signal(spgr.size());
    for (int i = 0; i < spgr.size(); i++) {
        const double a = spgr.FA[i] * B1;
        Mz.noalias()   = (Eigen::Matrix2d::Identity() - eATR * cos(a)).partialPivLu().solve(RHS);
        Mxy.noalias()  = echo * PD * Mz * sin(a);
        signal(i)      = Mxy(0) + Mxy(1);
    }
    QI_DBMSG("SPGR2\n");
    QI_DB(PD);
    QI_DB(T1_a);
    QI_DB(T2_a);
    QI_DB(T1_b);
    QI_DB(T2_b);
    QI_DB(tau_a);
    QI_DB(f_a);
    QI_DB(f0);
    QI_DB(B1);
    QI_DB(spgr.TR);
    QI_DB(spgr.TE);
    QI_DBVEC(spgr.FA);
    QI_DBVEC(signal);
    return signal;
}

Eigen::ArrayXcd SSFP2(double const            PD,
                      double const            T1_a,
                      double const            T2_a,
                      double const            T1_b,
                      double const            T2_b,
                      double const            tau_a,
                      double const            f_a,
                      double const            f0,
                      double const            B1,
                      QI::SSFPSequence const &s) {
    const double &TR   = s.TR;
    const double  E1_a = exp(-TR / T1_a);
    const double  E1_b = exp(-TR / T1_b);
    const double  E2_a = exp(-TR / T2_a);
    const double  E2_b = exp(-TR / T2_b);
    double        f_b, k_ab, k_ba;
    CalcExchange(tau_a, f_a, f_b, k_ab, k_ba);
    const double   E_ab  = exp(-TR * k_ab / f_b);
    const double   K1    = E_ab * f_b + f_a;
    const double   K2    = E_ab * f_a + f_b;
    const double   K3    = f_a * (1 - E_ab);
    const double   K4    = f_b * (1 - E_ab);
    Eigen::ArrayXd alpha = B1 * s.FA;
    Eigen::ArrayXd theta = s.PhaseInc + 2. * M_PI * f0 * TR;

    // TR Evolution
    Eigen::MatrixXd M(4, s.size());
    Eigen::Matrix6d LHS;
    Eigen::Vector6d RHS;
    RHS << 0, 0, 0, 0, -E1_b * K3 * f_b + f_a * (-E1_a * K1 + 1),
        -E1_a * K4 * f_a + f_b * (-E1_b * K2 + 1);

    // TE Evolution
    const double sE2_a    = exp(-TR / (2. * T2_a));
    const double sE2_b    = exp(-TR / (2. * T2_b));
    const double sqrtE_ab = exp(-TR * k_ab / (2 * f_b));
    const double K1e      = sqrtE_ab * f_b + f_a;
    const double K2e      = sqrtE_ab * f_a + f_b;
    const double K3e      = f_a * (1 - sqrtE_ab);
    const double K4e      = f_b * (1 - sqrtE_ab);
    const double cte      = cos(M_PI * f0 * TR);
    const double ste      = sin(M_PI * f0 * TR);

    Eigen::Matrix4d echo;
    echo << sE2_a * K1e * cte, sE2_b * K3e * cte, -sE2_a * K1e * ste, -sE2_b * K3e * ste, //
        sE2_a * K4e * cte, sE2_b * K2e * cte, -sE2_a * K4e * ste, -sE2_b * K2e * ste,     //
        sE2_a * K1e * ste, sE2_b * K3e * ste, sE2_a * K1e * cte, sE2_b * K3e * cte,       //
        sE2_a * K4e * ste, sE2_b * K2e * ste, sE2_a * K4e * cte, sE2_b * K2e * cte;

    for (int i = 0; i < s.size(); i++) {
        const double ca  = cos(alpha[i]);
        const double sa  = sin(alpha[i]);
        const double ctr = cos(theta[i]);
        const double str = sin(theta[i]);

        LHS << -E2_a * K1 * ctr + ca, -E2_b * K3 * ctr, E2_a * K1 * str, E2_b * K3 * str, sa, 0,
            -E2_a * K4 * ctr, -E2_b * K2 * ctr + ca, E2_a * K4 * str, E2_b * K2 * str, 0, sa,
            -E2_a * K1 * str, -E2_b * K3 * str, -E2_a * K1 * ctr + 1, -E2_b * K3 * ctr, 0, 0,
            -E2_a * K4 * str, -E2_b * K2 * str, -E2_a * K4 * ctr, -E2_b * K2 * ctr + 1, 0, 0, //
            -sa, 0, 0, 0, -E1_a * K1 + ca, -E1_b * K3,                                        //
            0, -sa, 0, 0, -E1_a * K4, -E1_b * K2 + ca;
        M.col(i) = echo * (LHS.partialPivLu().solve(RHS)).head(4);
    }

    Eigen::ArrayXcd mce(s.size());
    mce.real() = PD * (M.row(0) + M.row(1));
    mce.imag() = PD * (M.row(2) + M.row(3));

    QI_DBMSG("SSFP2\n");
    QI_DB(PD);
    QI_DB(T1_a);
    QI_DB(T2_a);
    QI_DB(T1_b);
    QI_DB(T2_b);
    QI_DB(tau_a);
    QI_DB(f_a);
    QI_DB(f0);
    QI_DB(B1);
    QI_DB(s.TR);
    QI_DBVEC(s.PhaseInc);
    QI_DBVEC(s.FA);
    QI_DBVEC(mce);
    return mce;
}

} // namespace QI
