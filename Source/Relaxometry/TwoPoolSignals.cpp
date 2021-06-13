#include "TwoPoolSignals.h"

#include "Helpers.h"
#include <unsupported/Eigen/MatrixFunctions>

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
    A << -(1. / T1_a) - k_ab, k_ba, k_ab, -(1. / T1_b) - k_ba;
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
        Mxy.noalias()  = echo * Mz * sin(a);
        signal(i)      = PD * (Mxy(0) + Mxy(1));
    }
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
                      QI::SSFPSequence const &ssfp) {
    const double &TR   = ssfp.TR;
    const double  E1_a = exp(-TR / T1_a);
    const double  E1_b = exp(-TR / T1_b);
    const double  E2_a = exp(-TR / T2_a);
    const double  E2_b = exp(-TR / T2_b);
    double        f_b, k_ab, k_ba;
    QI::CalcExchange(tau_a, f_a, f_b, k_ab, k_ba);
    QI_DB(PD);
    QI_DB(f_a);
    QI_DB(f_b);
    QI_DB(tau_a);
    QI_DB(k_ab);
    QI_DB(k_ba);

    // Evolution over TR
    const double         E_ab   = exp(-TR * k_ab / f_b);
    const double         K1     = E_ab * f_b + f_a;
    const double         K2     = E_ab * f_a + f_b;
    const double         K3     = f_a * (1 - E_ab);
    const double         K4     = f_b * (1 - E_ab);
    const Eigen::ArrayXd alpha  = B1 * ssfp.FA;
    const double         theta0 = 2. * M_PI * f0 * TR;
    const Eigen::ArrayXd theta  = theta0 + ssfp.PhaseInc;

    Eigen::Matrix6d LHS;
    Eigen::Vector6d RHS;
    RHS << 0, 0, 0, 0, -E1_b * K3 * f_b + f_a * (-E1_a * K1 + 1),
        -E1_a * K4 * f_a + f_b * (-E1_b * K2 + 1);

    // Evolution to TE
    const double sE2_a = sqrt(E2_a);
    const double sE2_b = sqrt(E2_b);
    const double sE_ab = exp(-TR * k_ab / (2 * f_b));
    const double eK1   = sE_ab * f_b + f_a;
    const double eK2   = sE_ab * f_a + f_b;
    const double eK3   = f_a * (1 - sE_ab);
    const double eK4   = f_b * (1 - sE_ab);
    const double cp    = cos(theta0 / 2.);
    const double sp    = sin(theta0 / 2.);

    Eigen::Matrix4d echo;
    echo << sE2_a * eK1 * cp, sE2_b * eK3 * cp, -sE2_a * eK1 * sp, -sE2_b * eK3 * sp, //
        sE2_a * eK4 * cp, sE2_b * eK2 * cp, -sE2_a * eK4 * sp, -sE2_b * eK2 * sp,     //
        sE2_a * eK1 * sp, sE2_b * eK3 * sp, sE2_a * eK1 * cp, sE2_b * eK3 * cp,       //
        sE2_a * eK4 * sp, sE2_b * eK2 * sp, sE2_a * eK4 * cp, sE2_b * eK2 * cp;

    // Now loop over flip-angles / phase-incs
    Eigen::ArrayXcd signal(ssfp.size());
    for (int i = 0; i < ssfp.size(); i++) {
        const double ca  = cos(alpha[i]);
        const double sa  = sin(alpha[i]);
        const double cta = cos(theta[i]);
        const double ctb = cos(theta[i]);
        const double sta = sin(theta[i]);
        const double stb = sin(theta[i]);

        LHS << -E2_a * K1 * cta + ca, -E2_b * K3 * cta, E2_a * K1 * sta, E2_b * K3 * sta, sa, 0, //
            -E2_a * K4 * ctb, -E2_b * K2 * ctb + ca, E2_a * K4 * stb, E2_b * K2 * stb, 0, sa,    //
            -E2_a * K1 * sta, -E2_b * K3 * sta, -E2_a * K1 * cta + 1, -E2_b * K3 * cta, 0, 0,    //
            -E2_a * K4 * stb, -E2_b * K2 * stb, -E2_a * K4 * ctb, -E2_b * K2 * ctb + 1, 0, 0,    //
            -sa, 0, 0, 0, -E1_a * K1 + ca, -E1_b * K3,                                           //
            0, -sa, 0, 0, -E1_a * K4, -E1_b * K2 + ca;
        Eigen::Vector4d const M = PD * echo * (LHS.partialPivLu().solve(RHS)).head(4);

        signal(i) = std::complex(M[0] + M[1], M[2] + M[3]);
    }
    return signal;
}

} // namespace QI
