/*
 *  RegionContraction.h
 *
 *  Created by Tobias Wood on 17/08/2012.
 *  Copyright (c) 2012-2013 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef DESPOT_RegionContraction_h
#define DESPOT_RegionContraction_h

#include <vector>
#include <random>
#include <iostream>
#include <atomic>
#include <cmath>

#include <Eigen/Dense>

#include "Util.h"

namespace QI {

typedef Eigen::Array<bool, Eigen::Dynamic, 1> ArrayXb;

std::vector<size_t> index_partial_sort(const Eigen::Ref<Eigen::ArrayXd> &x, Eigen::ArrayXd::Index N);
std::vector<size_t> index_partial_sort(const Eigen::Ref<Eigen::ArrayXd> &x, Eigen::ArrayXd::Index N)
{
	eigen_assert(x.size() >= N);
    std::vector<size_t> allIndices(x.size()), indices(N);
    for(size_t i = 0; i < allIndices.size(); i++) {
		allIndices[i] = i;
    }
	partial_sort(allIndices.begin(), allIndices.begin() + N, allIndices.end(),
	             [&x](size_t i1, size_t i2) { return x[i1] < x[i2]; });
    for (Eigen::ArrayXd::Index i = 0; i < N; i++) {
		indices[i] = allIndices[i];
	}
    return indices;
}

enum class RCStatus {
	NotStarted = -1,
	Converged, NoImprovement, IterationLimit, ErrorInvalid, ErrorResidual
};

std::ostream& operator<<(std::ostream &os, const RCStatus &s) {
	switch (s) {
		case RCStatus::NotStarted: os << "Not Started"; break;
		case RCStatus::Converged: os << "Converged"; break;
		case RCStatus::NoImprovement: os << "No Improvement to best residual"; break;
		case RCStatus::IterationLimit: os << "Reached iteration limit"; break;
		case RCStatus::ErrorInvalid: os << "Could not generate valid sample"; break;
		case RCStatus::ErrorResidual: os << "Infinite residual found"; break;
	}
	return os;
}

template <typename Functor_t>
class RegionContraction {
	private:
		Functor_t &m_f;
        std::mt19937_64 m_rng;
        Eigen::ArrayXXd m_startBounds, m_currentBounds;
        Eigen::ArrayXd m_threshes;
		size_t m_nS, m_nR, m_maxContractions, m_contractions;
		double m_expand, m_SoS;
		RCStatus m_status;
		bool m_gaussian, m_debug;
	
	public:
		RegionContraction(Functor_t &f,
                          const Eigen::Ref<Eigen::ArrayXXd> &startBounds, const Eigen::ArrayXd &thresh,
						  const int nS = 5000, const int nR = 50, const int maxContractions = 10,
						  const double expand = 0., const bool gauss = false, const bool debug = false, const int seed = -1) :
				m_f(f), m_startBounds(startBounds), m_currentBounds(startBounds),
				m_nS(nS), m_nR(nR), m_maxContractions(maxContractions),
                m_threshes(thresh), m_expand(expand), m_contractions(0),
                m_status(RCStatus::NotStarted), m_gaussian(gauss), m_debug(debug)
		{
			eigen_assert(f.inputs() == startBounds.rows());
			eigen_assert(startBounds.cols() == 2);
			eigen_assert(thresh.rows() == f.inputs());
			eigen_assert((thresh >= 0.).all() && (thresh <= 1.).all());

			if (seed < 0) {
                m_rng = std::mt19937_64(RandomSeed());
			} else {
                m_rng = std::mt19937_64(seed);
			}
		}
		
        const Eigen::ArrayXXd &startBounds() const { return m_startBounds; }
        void setBounds(const Eigen::Ref<Eigen::ArrayXXd> &b) {
			eigen_assert(m_f.inputs() == b.rows());
			eigen_assert(b.cols() == 2);
			m_startBounds = b;
		}
		const ArrayXb &thresholds() const { return m_threshes; }
        void setThresholds(const Eigen::Ref<Eigen::ArrayXd> &t) {
			eigen_assert(t.rows() == m_f.inputs());
			eigen_assert((t >= 0.).all() && (t <= 1.).all());
			m_threshes = t;
		}
        size_t   contractions() const { return m_contractions; }
        RCStatus       status() const { return m_status; }
        const Eigen::ArrayXXd &currentBounds() const { return m_currentBounds; }
		const double   SoS() const { return m_SoS; }
        Eigen::ArrayXd startWidth() const { return m_startBounds.col(1) - m_startBounds.col(0); }
        Eigen::ArrayXd width() const { return m_currentBounds.col(1) - m_currentBounds.col(0); }
        Eigen::ArrayXd midPoint() const { return (m_currentBounds.rowwise().sum() / 2.); }
		
        void optimise(Eigen::Ref<Eigen::ArrayXd> params) {
            static std::atomic<bool> finiteWarning(false);
            static std::atomic<bool> constraintWarning(false);
            static std::atomic<bool> boundsWarning(false);
            std::mutex warn_mtx;

			eigen_assert(m_f.inputs() == params.size());
			int nP = static_cast<int>(params.size());
            Eigen::ArrayXXd samples(m_f.inputs(), m_nS);
            Eigen::ArrayXXd retained(m_f.inputs(), m_nR);
            Eigen::ArrayXd residuals(m_nS);
            Eigen::ArrayXd retainedRes(m_nR);
            Eigen::ArrayXd gauss_mu(m_f.inputs()), gauss_sigma(m_f.inputs());
            std::vector<size_t> indices(m_nR);
			m_currentBounds = m_startBounds;
			if ((m_startBounds != m_startBounds).any() ||
                (m_startBounds >= std::numeric_limits<double>::infinity()).any() ||
			    (m_startBounds.col(1) < m_startBounds.col(0)).any()) {
				warn_mtx.lock();
				if (!boundsWarning) {
					boundsWarning = true;
                    std::cerr << "Warning: Starting boundaries do not make sense." << std::endl;
                    std::cerr << "Bounds were: " << m_startBounds.transpose() << std::endl;
                    std::cerr << "This warning will only be printed once." << std::endl;
				}
				warn_mtx.unlock();
				params.setZero();
				m_status = RCStatus::ErrorInvalid;

				return;
			}

			if (m_debug) {
                std::cout << std::endl;
                std::cout << "START REGION CONTRACTION" << std::endl;
                std::cout << "Start Boundaries: " << std::endl << m_startBounds.transpose() << std::endl;
			}
			
            std::uniform_real_distribution<double> uniform(0., 1.);
			m_status = RCStatus::IterationLimit;
			for (m_contractions = 0; m_contractions < m_maxContractions; m_contractions++) {
				size_t startSample = 0;
				/*if (m_contractions > 0) {
					// Keep the retained samples to prevent boundary contracting too fast
					for (size_t s = 0; s < m_nR; s++) {
						samples.col(s) = retained.col(s);
						residuals.col(s) = retainedRes.col(s);
					}
					startSample = m_nR;
				}*/
				for (size_t s = startSample; s < m_nS; s++) {
                    Eigen::ArrayXd tempSample(nP);
					size_t nTries = 0;
					do {
						if (!m_gaussian || (m_contractions == 0)) {
							for (int p = 0; p < nP; p++) {
								tempSample(p) = uniform(m_rng);
							}
							tempSample *= width();
							tempSample += m_currentBounds.col(0);
						} else {
							for (int p = 0; p < nP; p++) {
                                if (std::isfinite(gauss_sigma(p))) {
                                    std::normal_distribution<double> gauss(gauss_mu(p), gauss_sigma(p));
									do {
										tempSample(p) = gauss(m_rng);
									} while ((tempSample(p) < m_currentBounds(p, 0)) || (tempSample(p) > m_currentBounds(p, 1)));
								} else {
									tempSample(p) = gauss_mu(p);
								}
							}
						}
						nTries++;
						if (nTries > 100) {
							warn_mtx.lock();
							if (!constraintWarning) {
								constraintWarning = true;
                                std::cerr << "Warning: Cannot fulfill sample constraints after " << std::to_string(nTries) << " attempts, giving up." << std::endl
                                          << "Last attempt was: " << tempSample.transpose() << std::endl
                                          << "This warning will only be printed once." << std::endl;
							}
							warn_mtx.unlock();
							params.setZero();
							m_status = RCStatus::ErrorInvalid;
							return;
						}
					} while (!m_f.constraint(tempSample));
					
                    residuals[s] = m_f(tempSample);
                    if (!std::isfinite(residuals[s])) {
						warn_mtx.lock();
						if (!finiteWarning) {
							finiteWarning = true;
							std::cout << "Warning: Non-finite residual found!" << std::endl
                                      << "Result may be meaningless. This warning will only be printed once." << std::endl
                                      << "Parameters were " << tempSample.transpose() << std::endl;
						}
						warn_mtx.unlock();
						params = retained.col(0);
						m_status = RCStatus::ErrorResidual;
						return;
					}
					samples.col(s) = tempSample;
				}
                indices = index_partial_sort(residuals, m_nR);
                Eigen::ArrayXd previousBest = retained.col(0);
				for (size_t i = 0; i < m_nR; i++) {
					retained.col(i) = samples.col(indices[i]);
                    retainedRes(i) = residuals(indices[i]);
				}
				// Find the min and max for each parameter in the top nR samples
				m_currentBounds.col(0) = retained.rowwise().minCoeff();
				m_currentBounds.col(1) = retained.rowwise().maxCoeff();
				if (m_gaussian) {
					gauss_mu = retained.rowwise().mean();
					gauss_sigma = ((retained.colwise() - gauss_mu).square().rowwise().sum() / (m_f.inputs() - 1)).sqrt();
				}
				if (m_debug) {
                    std::cout << "CONTRACTION:    " << m_contractions << std::endl
                              << "Retained best: " << retainedRes.minCoeff() << " Worst: " << retainedRes.maxCoeff() << std::endl
                              << "All best:      " << residuals.minCoeff() << " Worst: " << residuals.maxCoeff() << std::endl
                              << "Current width%: " << (width() / startWidth()).transpose() << std::endl;
					//cout << "Thresh        : " << m_threshes.transpose() << std::endl;
					//cout << "Width < Thresh: " << (width() <= (m_threshes * startWidth())).transpose() << std::endl;
					//cout << "Converged:      " << (width() <= (m_threshes * startWidth())).all() << std::endl;
					if (m_gaussian) {
                        std::cout << "Gaussian mu:    " << gauss_mu.transpose() << std::endl
                                  << "Gaussian sigma%:"<<  (gauss_sigma / startWidth()).transpose() << std::endl;
					}
				}
				// Terminate if all the desired parameters have converged
				if ((width() <= (m_threshes * startWidth())).all()) {
					m_status = RCStatus::Converged;
					m_contractions++; // Just to give an accurate contraction count.
					break;
				} else if ((previousBest == retained.col(0)).all()) {
					m_status = RCStatus::NoImprovement;
					m_contractions++; // Just to give an accurate contraction count.
					break;
				}
				
				if (m_expand != 0) {
					// Expand the boundaries back out in case we just missed a minima,
					// but don't go past initial boundaries
                    Eigen::ArrayXd tempW = width(); // Because altering .col(0) will change width
					m_currentBounds.col(0) = (m_currentBounds.col(0) - tempW * m_expand).max(m_startBounds.col(0));
					m_currentBounds.col(1) = (m_currentBounds.col(1) + tempW * m_expand).min(m_startBounds.col(1));
					if (m_debug) {
                        std::cout << "Width expanded to: " << width().transpose() << std::endl;
					}
				}
			}

			if (m_gaussian) {
				params = gauss_mu;
			} else {
				// Return the best evaluated solution so far
				params = retained.col(0);
			}
            m_SoS = retainedRes(0);
			if (m_debug) {
                std::cout << "Finished, contractions = " << m_contractions << std::endl;
			}
		}
};

} // End namespace QI

#endif
