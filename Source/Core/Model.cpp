#include "Model.h"

#include <random>

namespace QI {

template <typename T>
auto RicianNoise<T>::Add(QI_ARRAY(T) const &s, T const sigma) -> QI_ARRAY(T) {
    QI_ARRAY(T) noise(s.rows());
    // Simple Box Muller transform
    QI_ARRAY(T) U = (QI_ARRAY(T)::Random(s.rows()) * 0.5) + 0.5;
    QI_ARRAY(T) V = (QI_ARRAY(T)::Random(s.rows()) * 0.5) + 0.5;
    QI_ARRAY(T) real     = (sigma / M_SQRT2) * (-2. * U.log()).sqrt() * cos(2. * M_PI * V);
    QI_ARRAY(T) imag     = (sigma / M_SQRT2) * (-2. * V.log()).sqrt() * sin(2. * M_PI * U);
    return ((s + real).abs2() + imag.abs2()).sqrt();
}
template struct RicianNoise<float>;
template struct RicianNoise<double>;

template <typename T>
auto ComplexNoise<T>::Add(QI_ARRAY(std::complex<T>) const &s,
                                             T const sigma) -> QI_ARRAY(std::complex<T>) {
    QI_ARRAY(std::complex<T>) noise(s.rows());
    // Simple Box Muller transform
    QI_ARRAY(T) U = (QI_ARRAY(T)::Random(s.rows()) * 0.5) + 0.5;
    QI_ARRAY(T) V = (QI_ARRAY(T)::Random(s.rows()) * 0.5) + 0.5;
    noise.real()     = (sigma / M_SQRT2) * (-2. * U.log()).sqrt() * cos(2. * M_PI * V);
    noise.imag()     = (sigma / M_SQRT2) * (-2. * V.log()).sqrt() * sin(2. * M_PI * U);
    return s + noise;
}
template struct ComplexNoise<float>;
template struct ComplexNoise<double>;

template <typename T> auto RealNoise<T>::Add(QI_ARRAY(T) const &s, T const sigma) -> QI_ARRAY(T) {
    return s + QI_ARRAY(T)::Random(s.rows()) * sigma;
}
template struct RealNoise<float>;
template struct RealNoise<double>;

} // namespace QI
