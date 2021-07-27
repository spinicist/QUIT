/* This code is adapted from Boost's bessel_j0.hpp
 * The original code is (c) 2006 Xioagang Zhang and subject to BoostLicense
 */

#pragma once

#include <array>

template <class T, class U, size_t count>
U evaluate_rational(const std::array<T, count> &num,
                    const std::array<T, count> &denom,
                    const U &                   z_) {
    U z(z_);
    U s1, s2;
    if (z <= 1.) {
        s1 = U(num[count - 1]);
        s2 = U(denom[count - 1]);
        for (int i = (int)count - 2; i >= 0; --i) {
            s1 *= z;
            s2 *= z;
            s1 += num[i];
            s2 += denom[i];
        }
    } else {
        z  = 1. / z;
        s1 = U(num[0]);
        s2 = U(denom[0]);
        for (unsigned i = 1; i < count; ++i) {
            s1 *= z;
            s2 *= z;
            s1 += num[i];
            s2 += denom[i];
        }
    }
    return s1 / s2;
}

template <typename T> T bessel_j0(T x) {
    static const std::array<double, 7> P1{-4.1298668500990866786e+11,
                                          2.7282507878605942706e+10,
                                          -6.2140700423540120665e+08,
                                          6.6302997904833794242e+06,
                                          -3.6629814655107086448e+04,
                                          1.0344222815443188943e+02,
                                          -1.2117036164593528341e-01};
    static const std::array<double, 7> Q1{2.3883787996332290397e+12,
                                          2.6328198300859648632e+10,
                                          1.3985097372263433271e+08,
                                          4.5612696224219938200e+05,
                                          9.3614022392337710626e+02,
                                          1.0,
                                          0.0};
    static const std::array<double, 8> P2{-1.8319397969392084011e+03,
                                          -1.2254078161378989535e+04,
                                          -7.2879702464464618998e+03,
                                          1.0341910641583726701e+04,
                                          1.1725046279757103576e+04,
                                          4.4176707025325087628e+03,
                                          7.4321196680624245801e+02,
                                          4.8591703355916499363e+01};
    static const std::array<double, 8> Q2{-3.5783478026152301072e+05,
                                          2.4599102262586308984e+05,
                                          -8.4055062591169562211e+04,
                                          1.8680990008359188352e+04,
                                          -2.9458766545509337327e+03,
                                          3.3307310774649071172e+02,
                                          -2.5258076240801555057e+01,
                                          1.0};
    static const std::array<double, 6> PC{2.2779090197304684302e+04,
                                          4.1345386639580765797e+04,
                                          2.1170523380864944322e+04,
                                          3.4806486443249270347e+03,
                                          1.5376201909008354296e+02,
                                          8.8961548424210455236e-01};
    static const std::array<double, 6> QC{2.2779090197304684318e+04,
                                          4.1370412495510416640e+04,
                                          2.1215350561880115730e+04,
                                          3.5028735138235608207e+03,
                                          1.5711159858080893649e+02,
                                          1.0};
    static const std::array<double, 6> PS{-8.9226600200800094098e+01,
                                          -1.8591953644342993800e+02,
                                          -1.1183429920482737611e+02,
                                          -2.2300261666214198472e+01,
                                          -1.2441026745835638459e+00,
                                          -8.8033303048680751817e-03};
    static const std::array<double, 6> QS{5.7105024128512061905e+03,
                                          1.1951131543434613647e+04,
                                          7.2642780169211018836e+03,
                                          1.4887231232283756582e+03,
                                          9.0593769594993125859e+01,
                                          1.0};
    static const double x1 = 2.4048255576957727686e+00, x2 = 5.5200781102863106496e+00,
                        x11 = 6.160e+02, x12 = -1.42444230422723137837e-03, x21 = 1.4130e+03,
                        x22 = 5.46860286310649596604e-04;

    T value, factor, r, rc, rs;

    if (x < 0.) {
        x = -x; // even function
    }
    if (x == 0.) {
        return static_cast<T>(1);
    }
    if (x <= 4.) // x in (0, 4]
    {
        T y    = x * x;
        r      = evaluate_rational(P1, Q1, y);
        factor = (x + x1) * ((x - x11 / 256.) - x12);
        value  = factor * r;
    } else if (x <= 8.0) // x in (4, 8]
    {
        T y    = 1. - (x * x) / 64.;
        r      = evaluate_rational(P2, Q2, y);
        factor = (x + x2) * ((x - x21 / 256.) - x22);
        value  = factor * r;
    } else // x in (8, \infty)
    {
        T y    = 8. / x;
        T y2   = y * y;
        rc     = evaluate_rational(PC, QC, y2);
        rs     = evaluate_rational(PS, QS, y2);
        factor = (1. / sqrt(M_PI)) / sqrt(x);
        //
        // What follows is really just:
        //
        // T z = x - pi/4;
        // value = factor * (rc * cos(z) - y * rs * sin(z));
        //
        // But using the addition formulae for sin and cos, plus
        // the special values for sin/cos of pi/4.
        //
        T sx  = sin(x);
        T cx  = cos(x);
        value = factor * (rc * (cx + sx) - y * rs * (sx - cx));
    }

    return value;
}