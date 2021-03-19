#pragma once

#define EIGEN_USE_THREADS

#include <unsupported/Eigen/CXX11/Tensor>
#include <unsupported/Eigen/CXX11/ThreadPool>

#include "Log.h"

template <typename T> typename T::Scalar Dot(T const &a, T const &b) {
    Eigen::TensorFixedSize<typename T::Scalar, Eigen::Sizes<>> d = (a.conjugate() * b).sum();
    return d();
}

template <typename T> float Norm(T const &a) {
    return sqrt(std::real(Dot(a, a)));
}

template <typename T>
auto ForwardDiff(T const &a, Eigen::Index const d, Eigen::array<long, 3> const dims) {
    using Dims = Eigen::array<long, 3>;
    Dims const sz{dims[0] - 2, dims[1] - 2, dims[2] - 2};
    Dims const st1{1, 1, 1};
    Dims       fwd{1, 1, 1};
    fwd[d] = 2;
    fmt::print(FMT_STRING("dims {} sz {} st1 {} fwd {}\n"),
               fmt::join(dims, ","),
               fmt::join(sz, ","),
               fmt::join(st1, ","),
               fmt::join(fwd, ","));
    return (a.slice(fwd, sz) - a.slice(st1, sz));
}

template <typename T>
auto BackwardDiff(T const &a, Eigen::Index const d, Eigen::array<long, 3> const dims) {
    using Dims = Eigen::array<long, 3>;
    Dims const sz{dims[0] - 2, dims[1] - 2, dims[2] - 2};
    Dims const st1{1, 1, 1};
    Dims       bck{1, 1, 1};
    bck[d] = 0;

    return (a.slice(st1, sz) - a.slice(bck, sz));
}

template <typename T>
auto CentralDiff(T const &a, Eigen::Index const d, Eigen::array<long, 3> const dims) {
    using Dims = Eigen::array<long, 3>;
    Dims const sz{dims[0] - 2, dims[1] - 2, dims[2] - 2};
    Dims const st1{1, 1, 1};
    Dims       fwd{1, 1, 1};
    Dims       bck{1, 1, 1};
    fwd[d] = 2;
    bck[d] = 0;

    return (a.slice(fwd, sz) - a.slice(bck, sz)) / a.slice(st1, sz).constant(2.f);
}

template <typename T>
void Grad(Eigen::Tensor<T, 3> const &a, Eigen::Tensor<T, 4> &g, Eigen::ThreadPoolDevice &dev) {
    using Dims = typename Eigen::Tensor<T, 3>::Dimensions;
    Dims const sz{a.dimension(0) - 2, a.dimension(1) - 2, a.dimension(2) - 2};
    Dims const st1{1, 1, 1};
    fmt::print(FMT_STRING("GRAD size {} st1 {}\n"), fmt::join(sz, ","), fmt::join(st1, ","));
    g.template chip<3>(0).slice(st1, sz).device(dev) = ForwardDiff(a, 0, a.dimensions());
    g.template chip<3>(1).slice(st1, sz).device(dev) = ForwardDiff(a, 1, a.dimensions());
    g.template chip<3>(2).slice(st1, sz).device(dev) = ForwardDiff(a, 2, a.dimensions());
}

template <typename T>
void Grad(Eigen::Tensor<T, 4> const &x, Eigen::Tensor<T, 4> &gx, Eigen::ThreadPoolDevice &dev) {
    using Dims = typename Eigen::Tensor<T, 3>::Dimensions;
    Dims const sz1{x.dimension(0), x.dimension(1), x.dimension(2)};
    Dims const sz2{x.dimension(0) - 2, x.dimension(1) - 2, x.dimension(2) - 2};
    Dims const st1{1, 1, 1};

    gx.template chip<3>(0).slice(st1, sz2).device(dev) =
        BackwardDiff(x.template chip<3>(0), 0, sz1);
    gx.template chip<3>(1).slice(st1, sz2).device(dev) =
        BackwardDiff(x.template chip<3>(1), 1, sz1);
    gx.template chip<3>(2).slice(st1, sz2).device(dev) =
        BackwardDiff(x.template chip<3>(2), 2, sz1);

    gx.template chip<3>(3).slice(st1, sz2).device(dev) =
        (BackwardDiff(x.template chip<3>(0), 1, sz1) +
         BackwardDiff(x.template chip<3>(1), 0, sz1)) /
        gx.template chip<3>(3).slice(st1, sz2).constant(2.f);

    gx.template chip<3>(4).slice(st1, sz2).device(dev) =
        (BackwardDiff(x.template chip<3>(0), 2, sz1) +
         BackwardDiff(x.template chip<3>(2), 0, sz1)) /
        gx.template chip<3>(4).slice(st1, sz2).constant(2.f);

    gx.template chip<3>(5).slice(st1, sz2).device(dev) =
        (BackwardDiff(x.template chip<3>(1), 2, sz1) +
         BackwardDiff(x.template chip<3>(2), 1, sz1)) /
        gx.template chip<3>(5).slice(st1, sz2).constant(2.f);
}

template <typename T>
inline void
Div(Eigen::Tensor<T, 4> const &x, Eigen::Tensor<T, 3> &div, Eigen::ThreadPoolDevice &dev) {
    using Dims = typename Eigen::Tensor<T, 3>::Dimensions;
    Dims const sz1{x.dimension(0), x.dimension(1), x.dimension(2)};
    Dims const sz2{x.dimension(0) - 2, x.dimension(1) - 2, x.dimension(2) - 2};
    Dims const st1{1, 1, 1};
    div.slice(st1, sz2).device(dev) = BackwardDiff(x.template chip<3>(0), 0, sz1) +
                                      BackwardDiff(x.template chip<3>(1), 1, sz1) +
                                      BackwardDiff(x.template chip<3>(2), 2, sz1);
}

template <typename T>
inline void
Div(Eigen::Tensor<T, 4> const &x, Eigen::Tensor<T, 4> &div, Eigen::ThreadPoolDevice &dev) {
    using Dims = typename Eigen::Tensor<T, 3>::Dimensions;
    Dims const sz1{x.dimension(0), x.dimension(1), x.dimension(2)};
    Dims const sz2{x.dimension(0) - 2, x.dimension(1) - 2, x.dimension(2) - 2};
    Dims const st1{1, 1, 1};
    div.template chip<3>(0).slice(st1, sz2).device(dev) =
        ForwardDiff(x.template chip<3>(0), 0, sz1) + ForwardDiff(x.template chip<3>(3), 1, sz1) +
        ForwardDiff(x.template chip<3>(4), 2, sz1);
    div.template chip<3>(1).slice(st1, sz2).device(dev) =
        ForwardDiff(x.template chip<3>(3), 0, sz1) + ForwardDiff(x.template chip<3>(1), 1, sz1) +
        ForwardDiff(x.template chip<3>(5), 2, sz1);
    div.template chip<3>(2).slice(st1, sz2).device(dev) =
        ForwardDiff(x.template chip<3>(4), 0, sz1) + ForwardDiff(x.template chip<3>(5), 1, sz1) +
        ForwardDiff(x.template chip<3>(2), 2, sz1);
}

template <typename T> auto ProjectP(Eigen::Tensor<T, 4> const &p, float const a) {
    Eigen::IndexList<int, int, int, Eigen::type2index<1>> res;
    res.set(0, p.dimension(0));
    res.set(1, p.dimension(1));
    res.set(2, p.dimension(2));
    Eigen::IndexList<Eigen::type2index<1>,
                     Eigen::type2index<1>,
                     Eigen::type2index<1>,
                     Eigen::type2index<3>>
        brd;

    Eigen::Tensor<float, 3> const normp =
        (p * p.conjugate()).sum(Eigen::array<int, 1>{3}).real().sqrt() / a;
    fmt::print(FMT_STRING("normp {}\n"), fmt::join(normp.dimensions(), ","));
    Eigen::Tensor<float, 3> const clamped = (normp > 1.f).select(normp, normp.constant(1.f));
    fmt::print(FMT_STRING("clamped {}\n"), fmt::join(clamped.dimensions(), ","));
    Eigen::Tensor<T, 4> result = p / clamped.reshape(res).broadcast(brd).template cast<T>();
    fmt::print(FMT_STRING("result {}\n"), fmt::join(result.dimensions(), ","));
    return result;
}

template <typename T> auto ProjectQ(Eigen::Tensor<T, 4> const &q, float const a) {
    Eigen::IndexList<int, int, int, Eigen::type2index<1>> res;
    res.set(0, q.dimension(0));
    res.set(1, q.dimension(1));
    res.set(2, q.dimension(2));
    Eigen::IndexList<Eigen::type2index<1>,
                     Eigen::type2index<1>,
                     Eigen::type2index<1>,
                     Eigen::type2index<6>>
        brd;
    using Dims4 = typename Eigen::Tensor<T, 4>::Dimensions;
    auto const q1 =
        q.slice(Dims4{0, 0, 0, 0}, Dims4{q.dimension(0), q.dimension(1), q.dimension(2), 3});
    auto const q2 =
        q.slice(Dims4{0, 0, 0, 3}, Dims4{q.dimension(0), q.dimension(1), q.dimension(2), 3});
    Eigen::Tensor<float, 3> const normq =
        ((q1 * q1.conjugate()).sum(Eigen::array<int, 1>{3}).real() +
         (q2 * q2.conjugate()).sum(Eigen::array<int, 1>{3}).real() * 2.f)
            .sqrt() /
        a;
    fmt::print(FMT_STRING("normq {}\n"), fmt::join(normq.dimensions(), ","));
    Eigen::Tensor<float, 3> const clamped = (normq > 1.f).select(normq, normq.constant(1.f));
    fmt::print(FMT_STRING("clamped {}\n"), fmt::join(clamped.dimensions(), ","));
    Eigen::Tensor<T, 4> result = q / clamped.reshape(res).broadcast(brd).template cast<T>();
    fmt::print(FMT_STRING("result {}\n"), fmt::join(result.dimensions(), ","));
    return result;
}

template <typename T>
Eigen::Tensor<T, 3> tgvdenoise(Eigen::Tensor<T, 3> const &image,
                               long const                 max_its,
                               float const                thresh,
                               float const                alpha,
                               float const                reduction,
                               float const                step_size,
                               bool                       vb,
                               Eigen::ThreadPoolDevice &  dev) {
    using T3                     = Eigen::Tensor<T, 3>;
    using T4                     = Eigen::Tensor<T, 4>;
    typename T3::Dimensions dims = image.dimensions();
    typename T4::Dimensions dims3{dims[0], dims[1], dims[2], 3};
    typename T4::Dimensions dims6{dims[0], dims[1], dims[2], 6};

    float const scale = Norm(image);
    // Primal Variables
    T3 u     = image / image.constant(scale);
    T3 u_    = u;
    T3 u_old = u;
    T4 grad_u(dims3);
    T4 v(dims3);
    T4 v_(dims3);
    T4 v_old(dims3);
    T4 grad_v(dims6); // Symmetric rank-2 tensor stored as xx yy zz (xy + yx) (xz + zx) (yz + zy)
    grad_u.setZero();
    v_.setZero();
    v_old.setZero();
    grad_v.setZero();

    // Dual Variables
    T4 p(dims3);
    T4 q(dims6);
    T4 divq(dims3);
    T3 divp(dims);
    p.setZero();
    divp.setZero();
    q.setZero();
    divq.setZero();

    float const alpha00 = alpha;
    float const alpha10 = alpha / 2.f;
    float const alpha01 = alpha00 * reduction;
    float const alpha11 = alpha10 * reduction;

    // Step lengths
    float const tau_p = 1.f / step_size;
    float const tau_d = 1.f / (step_size / 2.f);

    QI::Info(vb, "TGV Scale {}", scale);

    for (auto ii = 0.f; ii < max_its; ii++) {
        // log.image(u, fmt::format(FMT_STRING("tgv-u-{:02}.nii"), ii));
        // Regularisation factors
        float const prog   = static_cast<float>(ii) / ((max_its == 1) ? 1. : (max_its - 1.f));
        float const alpha0 = std::exp(std::log(alpha01) * prog + std::log(alpha00) * (1.f - prog));
        float const alpha1 = std::exp(std::log(alpha11) * prog + std::log(alpha10) * (1.f - prog));

        // Update p
        fmt::print("UPDATE P\n");
        fmt::print(FMT_STRING("u_ {} grad_u {}\n"),
                   fmt::join(u_.dimensions(), ","),
                   fmt::join(grad_u.dimensions(), ","));
        Grad(u_, grad_u, dev);
        fmt::print("Finished grad\n");
        p.device(dev) =
            p - p.constant(tau_d) * (grad_u + v_); // Paper says +tau, but code says -tau
        p.device(dev) = ProjectP(p, alpha1);

        // Update q
        fmt::print("UPDATE Q\n");
        Grad(v_, grad_v, dev);
        q.device(dev) = q - q.constant(tau_d) * grad_v;
        q.device(dev) = ProjectQ(q, alpha0);

        // Update u
        fmt::print("UPDATE U\n");
        u_old.device(dev) = u;
        Div(p, divp, dev);
        u.device(dev) = ((u - u.constant(tau_p) * divp) + (image * u.constant(tau_p / scale))) /
                        u.constant(1.f + tau_p); // Prox op
        u_.device(dev) = 2.0 * u - u_old;

        // Update v
        fmt::print("UPDATE V\n");
        v_old.device(dev) = v;
        Div(q, divq, dev);
        v.device(dev)  = v - v.constant(tau_p) * (divq - p);
        v_.device(dev) = 2.0 * v - v_old;

        // Check for convergence
        fmt::print("CONVERGENCE\n");
        float const delta = Norm(u - u_old);
        QI::Info(vb, FMT_STRING("TGV {}: ɑ0 {:.2g} ɑ1 {:.2g} δ {}"), ii + 1, alpha0, alpha1, delta);
        if (delta < thresh) {
            QI::Info(vb, "Reached threshold on delta, stopping");
            break;
        }
    }

    u = u * u.constant(scale);
    return u;
}