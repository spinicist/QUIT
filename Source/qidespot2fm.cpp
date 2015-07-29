/*
 *  despot2fm.cpp
 *
 *  Created by Tobias Wood on 2015/06/03.
 *  Copyright (c) 2015 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <time.h>
#include <getopt.h>
#include <iostream>
#include <atomic>
#include <Eigen/Dense>
#include <unsupported/Eigen/LevenbergMarquardt>
#include <unsupported/Eigen/NumericalDiff>

#include "Util.h"
#include "Filters/ApplyAlgorithmSliceBySliceFilter.h"
#include "Model.h"
#include "Sequence.h"

#include "itkAmoebaOptimizer.h"
#include "itkLBFGSBOptimizer.h"

using namespace std;
using namespace Eigen;

/* The code below is really quite hairy. It relies on template specialisations
 * to ensure the correct behaviour when fitting to complex or magnitude data.
 * The central issue is the DifferenceVector functions, because we have to take
 * the .abs() in a different place with complex or magnitude data.
 *
 * Everything else then becomes tedious C++ to ensure the right version of
 * these functions is called.
 */

template<typename T> ArrayXd DifferenceVector(Ref<const ArrayXcd> a1, Ref<const ArrayXcd> a2, T dummy);
template<typename T>
ArrayXd DifferenceVector(Ref<const ArrayXcd> a1, const Array<T, Dynamic, 1> &a2) {
    //cout << __PRETTY_FUNCTION__ << endl;
    return a1.abs() - a2.abs();
}

template<typename T>
ArrayXd DifferenceVector(Ref<const ArrayXcd> a1, const Array<complex<T>, Dynamic, 1> &a2) {
    //cout << __PRETTY_FUNCTION__ << endl;
    return (a1 - a2).abs();
}

template<typename T>
class FMFunctor : public DenseFunctor<double> {
public:
    typedef Array<T, Eigen::Dynamic, 1> TArray;

	const shared_ptr<SequenceBase> m_sequence;
	shared_ptr<SCD> m_model;
    TArray m_data;
	const double m_T1, m_B1;

    FMFunctor(const shared_ptr<SCD> m, const shared_ptr<SequenceBase> s, const TArray &d, const double T1, const double B1) :
		DenseFunctor<double>(3, s->size()),
		m_model(m), m_sequence(s), m_data(d),
		m_T1(T1), m_B1(B1)
	{
		assert(static_cast<size_t>(m_data.rows()) == values());
	}

	const bool constraint(const VectorXd &params) const {
		Array4d fullparams;
		fullparams << params(0), m_T1, params(1), params(2);
		return m_model->ValidParameters(fullparams);
	}

    int operator()(const Ref<VectorXd> &params, Ref<ArrayXd> diffs) const {
        //cout << __PRETTY_FUNCTION__ << endl;
        eigen_assert(diffs.size() == values());
        ArrayXd fullparams(5);
        fullparams << params(0), m_T1, params(1), params(2), m_B1;
        ArrayXcd s = m_sequence->signal(m_model, fullparams);
        diffs = DifferenceVector(s, m_data);
        return 0;
    }
};

template<typename T>
class FixT2 : public DenseFunctor<double> {
public:
    typedef Array<T, Eigen::Dynamic, 1> TArray;

	const shared_ptr<SequenceBase> m_sequence;
	shared_ptr<SCD> m_model;
    TArray m_data;
	const double m_T1, m_B1;
	double m_T2;

    FixT2(const shared_ptr<SCD> m, const shared_ptr<SequenceBase> s, const TArray &d, const double T1, const double T2, const double B1) :
		DenseFunctor<double>(2, s->size()),
		m_model(m), m_sequence(s), m_data(d),
		m_T1(T1), m_T2(T2), m_B1(B1)
	{
		assert(static_cast<size_t>(m_data.rows()) == values());
	}

	void setT2(double T2) { m_T2 = T2; }
	int operator()(const Ref<VectorXd> &params, Ref<ArrayXd> diffs) const {
        //cout << __PRETTY_FUNCTION__ << endl;
		eigen_assert(diffs.size() == values());

		ArrayXd fullparams(5);
		fullparams << params(0), m_T1, m_T2, params(1), m_B1;
		ArrayXcd s = m_sequence->signal(m_model, fullparams);
        diffs = DifferenceVector(s, m_data);
		return 0;
	}
};

template<typename T>
class FMAlgo : public Algorithm<T> {
protected:
	const shared_ptr<SCD> m_model = make_shared<SCD>();
	shared_ptr<SSFPSimple> m_sequence;

public:
    typedef typename Algorithm<T>::TArray TArray;
    typedef typename Algorithm<T>::TInput TInput;
    typedef typename Algorithm<T>::TIterations TIterations;

    void setSequence(shared_ptr<SSFPSimple> s) { m_sequence = s; }

	size_t numInputs() const override  { return m_sequence->count(); }
	size_t numConsts() const override  { return 2; }
	size_t numOutputs() const override { return 3; }
	size_t dataSize() const override   { return m_sequence->size(); }

	virtual TArray defaultConsts() {
		// T1 & B1
		TArray def = TArray::Ones(2);
		return def;
	}
};

template<typename T>
class LMAlgo : public FMAlgo<T> {
protected:
    static const int m_iterations = 100;
public:
    typedef typename FMAlgo<T>::TArray TArray;
    typedef typename FMAlgo<T>::TInput TInput;
    typedef typename FMAlgo<T>::TIterations TIterations;

    void f0guess(double &lo, double &hi, double &step, const TInput &data) const;
    virtual void apply(const TInput &data, const TArray &inputs, TArray &outputs, TArray &resids, TIterations &its) const override;
};

template<typename T>
void LMAlgo<T>::f0guess(double &lo, double &hi, double &step, const TInput &data) const {
    double bw = 1. / (4. * this->m_sequence->TR());
    lo = -bw + 1.;
    hi = bw + 2.;
    step = bw;
    if (this->m_sequence->isSymmetric()) {
        lo = 1.;
    }
    //cout << __PRETTY_FUNCTION__ << endl;
    //cout << lo << "/" << hi << "/" << step << endl;
}

template<>
void LMAlgo<complex<double>>::f0guess(double &lo, double &hi, double &step, const TInput &data) const {
    double bw = 1. / (4. * this->m_sequence->TR());
    lo = -bw + 1.;
    hi = bw + 2.;
    step = bw;
    if (this->m_sequence->isSymmetric()) {
        lo = 1.;
    }
    //cout << endl << __PRETTY_FUNCTION__ << endl;
    //cout << lo << "/" << a << "/" << hi << "/" << step << endl;
}

template<typename T>
void LMAlgo<T>::apply(const TInput &data, const TArray &inputs, TArray &outputs, TArray &resids, TIterations &its) const
{
    //cout << __PRETTY_FUNCTION__ << endl;
    const double T1 = inputs[0];
    if (isfinite(T1) && (T1 > 0.001)) {
        const double B1 = inputs[1];
        double bestF = numeric_limits<double>::infinity();
        for (int j = 0; j < 2; j++) {
            const double T2guess = (0.05 + j * 0.2) * T1; // From a Yarnykh paper T2/T1 = 0.045 in brain at 3T. Try the longer value for CSF
            double lo, hi, step;
            this->f0guess(lo, hi, step, data);
            //cout << lo << "/" << hi << "/" << step << endl;
            its = 0;
            for (float f0guess = lo; f0guess < hi; f0guess += step) {
                // First fix T2 and fit
                FixT2<T> fixT2(this->m_model, this->m_sequence, data, T1, T2guess, B1);
                NumericalDiff<FixT2<T>> fixT2Diff(fixT2);
                LevenbergMarquardt<NumericalDiff<FixT2<T>>> fixT2LM(fixT2Diff);
                fixT2LM.setMaxfev(this->m_iterations * (this->m_sequence->size() + 1));
                VectorXd g(2); g << data.abs().maxCoeff() * 10.0, f0guess;
                //cout << "T1 " << T1 << " T2 " << T2guess << " g " << g.transpose() << endl;
                fixT2LM.minimize(g);

                // Now fit everything together
                FMFunctor<T> full(this->m_model, this->m_sequence, data, T1, B1);
                NumericalDiff<FMFunctor<T>> fullDiff(full);
                LevenbergMarquardt<NumericalDiff<FMFunctor<T>>> fullLM(fullDiff);
                VectorXd fullP(3); fullP << g[0], T2guess, g[1]; // Now include T2
                //cout << "Before " << fullP.transpose() << endl;
                fullLM.minimize(fullP);

                double F = fullLM.fnorm();
                //cout << "After  " << fullP.transpose() << " F " << F << endl;
                if (F < bestF) {
                    outputs = fullP;
                    bestF = F;
                }
                its += fixT2LM.iterations() + fullLM.iterations();
            }
        }
        outputs[1] = clamp(outputs[1], 0.001, T1);
        VectorXd pfull(5); pfull << outputs[0], T1, outputs[1], outputs[2], B1; // Now include EVERYTHING to get a residual
        ArrayXcd theory = this->m_sequence->signal(this->m_model, pfull);
        resids = DifferenceVector(theory, data);
    } else {
        // No point in processing -ve T1
        outputs.setZero();
        resids.setZero();
        its = 0;
    }
}

namespace itk {
template<typename T>
class FMCostFunction : public SingleValuedCostFunction
{
public:
    /** Standard class typedefs. */
    typedef FMCostFunction            Self;
    typedef SingleValuedCostFunction  Superclass;
    typedef SmartPointer<Self>        Pointer;
    typedef SmartPointer<const Self>  ConstPointer;

    typedef Eigen::Array<T, Dynamic, 1> TArray;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(CostFunction, SingleValuedCostfunction);

    TArray m_data;
    double m_T1, m_B1;
    shared_ptr<SequenceBase> m_sequence;
    shared_ptr<SCD> m_model;

    unsigned int GetNumberOfParameters(void) const { return 3; } // itk::CostFunction
    ArrayXd residuals(const ParametersType &p) const {
        ArrayXd fullparams(5);
        fullparams << p(0), m_T1, p(1), p(2), m_B1;
        ArrayXcd s = m_sequence->signal(m_model, fullparams);
        //cout << "fullp " << fullparams.transpose() << endl;
        return DifferenceVector(s, m_data);
    }
    MeasureType GetValue(const ParametersType &p1) const {
        double r = residuals(p1).square().sum();
        //cout << "r " << r << endl;
        return r;
    }
    void GetDerivative(const ParametersType &p, DerivativeType &df) const {
        using std::sqrt;
        using std::abs;
        double eps = sqrt(numeric_limits<double>::epsilon());
        df.SetSize(p.Size());
        ParametersType x = p;
        for (int i = 0; i < p.Size(); ++i) {
            double h = eps * abs(x[i]);
            if (h == 0.) {
                h = eps;
            }
            double _x = x[i];
            x[i] += h;
            double v2 = GetValue(x);
            x[i] -= 2*h;
            double v1 = GetValue(x);
            df[i] = (v2 - v1)/(2*h);
            //cout << "h " << h << " x " << x << " v1 " << v1 << " v2 " << v2 << endl;
            x[i] = _x;
        }
        //cout << "df " << df << endl;
    }

    protected:
    FMCostFunction(){};
    ~FMCostFunction(){};

    private:
    FMCostFunction(const Self &); //purposely not implemented
    void operator = (const Self &); //purposely not implemented
};
}

template<typename T>
class LBFGSBAlgo : public FMAlgo<T> {
public:
    using typename FMAlgo<T>::TArray;
    using typename FMAlgo<T>::TInput;
    using typename FMAlgo<T>::TIterations;
    typedef typename itk::LBFGSBOptimizer TOptimizer;

    virtual void apply(const TInput &indata, const TArray &consts, TArray &outputs, TArray &resids, TIterations &its) const override {
        double T1 = consts[0];
        double B1 = consts[1];
        if (isfinite(T1) && (T1 > 0.001)) {
            // Improve scaling by dividing the PD down to something sensible.
            // This gets scaled back up at the end.
            const auto data = indata / indata.abs().maxCoeff();
            TOptimizer::Pointer optimizer = TOptimizer::New();
            // Set properties pertinent to convergence
            optimizer->SetCostFunctionConvergenceFactor(1.e3);
            optimizer->SetProjectedGradientTolerance(1e-10);
            optimizer->SetMaximumNumberOfIterations(100);
            optimizer->SetMaximumNumberOfEvaluations(9999);
            optimizer->SetMaximumNumberOfCorrections(10);
            // Instantiate the cost function
            auto cost = itk::FMCostFunction<T>::New();
            cost->m_data = data;
            cost->m_T1 = T1;
            cost->m_B1 = B1;
            cost->m_sequence = this->m_sequence;
            cost->m_model = this->m_model;
            optimizer->SetCostFunction(cost.GetPointer());

            TOptimizer::BoundSelectionType select(3);
            TOptimizer::BoundValueType lower(3), upper(3);
            double f0_lo, f0_hi, f0_step;
            if (this->m_sequence->isSymmetric()) {
                select[0] = 1; select[1] = 2; select[2] = 2; // Lower bounds for PD, upper and lower for f0 & T2
                lower[0] = 0.001; lower[1] = this->m_sequence->TR() * 2.0; lower[2] = 0.001;
                upper[0] = 0; upper[1] = T1; upper[2] = 0.6/this->m_sequence->TR(); // Allow for a bit of fuzz on upper boundary
                f0_lo = 1.0;
                f0_step = 0.4 / this->m_sequence->TR();
                f0_hi = f0_lo + 1. + f0_step;
            } else {
                select[0] = 1; select[1] = 2; select[2] = 2; // Lower bounds for PD, upper and lower for f0 & T2
                lower[0] = 0.001; lower[1] = this->m_sequence->TR() * 2.0; lower[2] = -0.5/this->m_sequence->TR();
                upper[0] = 0; upper[1] = T1; upper[2] = 0.5/this->m_sequence->TR();
                f0_lo = -(0.4 / this->m_sequence->TR());
                f0_hi =  (0.4 / this->m_sequence->TR()) + 1.;
                f0_step = 4./15. / this->m_sequence->TR(); // 2/3 * 0.4
            }
            optimizer->SetLowerBound(lower);
            optimizer->SetUpperBound(upper);
            optimizer->SetBoundSelection(select);
            double best = numeric_limits<double>::infinity();
            TOptimizer::ParametersType bestP;
            its = 0;
            for (double f0 = f0_lo; f0 < f0_hi; f0 += f0_step) {
                TOptimizer::ParametersType p(3);
                p[0] = 10.; p[1] = 0.1 * T1; p[2] = f0; // Yarnykh gives T2 = 0.045 * T1 in brain, but best to overestimate for CSF
                optimizer->SetInitialPosition(p);
                optimizer->StartOptimization();
                p = optimizer->GetCurrentPosition();
                double r = cost->GetValue(p);
                if (r < best) {
                    best = r;
                    bestP = p;
                }
                its += optimizer->GetCurrentIteration();
            }
            outputs[0] = bestP[0] * indata.abs().maxCoeff();
            outputs[1] = bestP[1];
            outputs[2] = bestP[2];
            resids = cost->residuals(bestP) * indata.abs().maxCoeff();
        } else {
            outputs.setZero();
            resids.setZero();
            its = 0;
        }
    }
};


const string usage {
"Usage is: despot2-fm [options] T1_map ssfp_file\n\
\
Options:\n\
    --help, -h        : Print this message\n\
    --verbose, -v     : Print slice processing times\n\
    --no-prompt, -n   : Suppress input prompts\n\
    --mask, -m file   : Mask input with specified file\n\
    --out, -o path    : Add a prefix to the output filenames\n\
    --B1, -b file     : B1 Map file (ratio)\n\
    --algo, -a l      : Use 2-step LM algorithm\n\
               b      : Use LBFGSB algorithm (default)\n\
    --complex, -x     : Fit to complex data\n\
    --start, -s N     : Start processing from slice N\n\
    --stop, -p  N     : Stop processing at slice N\n\
    --flip, -F        : Data order is phase, then flip-angle (default opposite)\n\
    --finite          : Use finite pulse length correction\n\
    --resids, -r      : Write out per flip-angle residuals\n\
    --threads, -T N   : Use N threads (default=hardware limit)\n"
};
/* --complex, -x     : Fit to complex data\n\ */

struct option long_opts[] = {
    {"help", no_argument, 0, 'h'},
    {"verbose", no_argument, 0, 'v'},
    {"no-prompt", no_argument, 0, 'n'},
    {"mask", required_argument, 0, 'm'},
    {"out", required_argument, 0, 'o'},
    {"B1", required_argument, 0, 'b'},
    {"algo", required_argument, 0, 'a'},
    {"complex", no_argument, 0, 'x'},
    {"start", required_argument, 0, 's'},
    {"stop", required_argument, 0, 'p'},
    {"flip", required_argument, 0, 'F'},
    {"threads", required_argument, 0, 'T'},
    {"finite", no_argument, 0, 'f'},
    {"resids", no_argument, 0, 'r'},
    {0, 0, 0, 0}
};
const char* short_opts = "hvnm:o:b:a:xs:p:FT:frd:";
int indexptr = 0;
char c;

template<typename T>
int run_main(int argc, char **argv) {
    typedef itk::Image<T, 4> TSeries;
    typedef itk::VectorImage<T, 3> TVectorImage;
    typedef itk::ImageFileReader<TSeries> TReader;
    typedef itk::ImageToVectorFilter<TSeries> TToVector;
    typedef itk::ReorderVectorFilter<TVectorImage> TReorder;
    typedef itk::ApplyAlgorithmSliceBySliceFilter<FMAlgo<T>, T, float, 3> TApply;

    int start_slice = 0, stop_slice = 0;
    int verbose = false, prompt = true, all_residuals = false,
        fitFinite = false, flipData = false, use_LBFGSB = true;
    string outPrefix;
    QI::ReadImageF::Pointer mask = ITK_NULLPTR, B1 = ITK_NULLPTR;

    optind = 1;
    while ((c = getopt_long(argc, argv, short_opts, long_opts, &indexptr)) != -1) {
        switch (c) {
        case 'x': case 'h': break; //Already handled in main
        case 'v': verbose = true; break;
        case 'n': prompt = false; break;
        case 'a':
        switch (*optarg) {
            case 'l': use_LBFGSB = false; if (verbose) cout << "LM algorithm selected." << endl; break;
            case 'b': use_LBFGSB = true; if (verbose) cout << "LBFGSB algorithm selected." << endl; break;
            default: throw(runtime_error(string("Unknown algorithm type ") + optarg)); break;
        } break;
        case 'm':
            if (verbose) cout << "Reading mask file " << optarg << endl;
            mask = QI::ReadImageF::New();
            mask->SetFileName(optarg);
            break;
        case 'o':
            outPrefix = optarg;
            if (verbose) cout << "Output prefix will be: " << outPrefix << endl;
            break;
        case 'b':
            if (verbose) cout << "Reading B1 file: " << optarg << endl;
            B1 = QI::ReadImageF::New();
            B1->SetFileName(optarg);
            break;
        case 's': start_slice = atoi(optarg); break;
        case 'p': stop_slice = atoi(optarg); break;
        case 'F': flipData = true; break;
        case 'f': fitFinite = true; break;
        case 'T': itk::MultiThreader::SetGlobalMaximumNumberOfThreads(atoi(optarg)); break;
        case 'r': all_residuals = true; break;
        case 0: break; // Just a flag
        default:
            cout << "Unhandled option " << string(1, c) << endl;
            return EXIT_FAILURE;
        }
    }
    if ((argc - optind) != 2) {
        cout << usage << endl;
        cout << "Wrong number of arguments. Need a T1 map and one SSFP file." << endl;
        return EXIT_FAILURE;
    }

    shared_ptr<SSFPSimple> ssfpSequence;
    if (fitFinite) {
        cout << "Using finite pulse model." << endl;
        ssfpSequence = make_shared<SSFPFinite>(prompt);
    } else {
        ssfpSequence = make_shared<SSFPSimple>(prompt);
    }
    if (verbose) cout << *ssfpSequence << endl;

    if (verbose) cout << "Reading T1 Map from: " << argv[optind] << endl;
    auto T1 = QI::ReadImageF::New();
    T1->SetFileName(argv[optind++]);
    if (verbose) cout << "Opening SSFP file: " << argv[optind] << endl;
    auto ssfpFile = TReader::New();
    auto ssfpData = TToVector::New();
    auto ssfpFlip = TReorder::New();
    ssfpFile->SetFileName(argv[optind++]);
    ssfpData->SetInput(ssfpFile->GetOutput());
    ssfpFlip->SetInput(ssfpData->GetOutput());
    if (flipData) {
        ssfpFlip->SetStride(ssfpSequence->phases());
    }
    auto apply = TApply::New();
    shared_ptr<FMAlgo<T>> algo;
    if (use_LBFGSB) {
        itk::MultiThreader::SetGlobalMaximumNumberOfThreads(1);
        algo = make_shared<LBFGSBAlgo<T>>();
    } else {
        algo = make_shared<LMAlgo<T>>();
    }
    algo->setSequence(ssfpSequence);
    apply->SetAlgorithm(algo);
    apply->SetInput(0, ssfpFlip->GetOutput());
    apply->SetConst(0, T1->GetOutput());
    apply->SetSlices(start_slice, stop_slice);
    if (B1) {
        apply->SetConst(1, B1->GetOutput());
    }
    if (mask) {
        apply->SetMask(mask->GetOutput());
    }
    time_t startTime;
    if (verbose) {
        startTime = QI::printStartTime();
        auto monitor = QI::SliceMonitor<TApply>::New();
        apply->AddObserver(itk::IterationEvent(), monitor);
    }
    apply->Update();
    if (verbose) {
        QI::printElapsedTime(startTime);
        cout << "Writing output files. Prefix is " << outPrefix << endl;
    }
    outPrefix = outPrefix + "FM_";
    QI::writeResult(apply->GetOutput(0), outPrefix + "PD.nii");
    QI::writeResult(apply->GetOutput(1), outPrefix + "T2.nii");
    QI::writeResult(apply->GetOutput(2), outPrefix + "f0.nii");
    QI::writeResiduals(apply->GetResidOutput(), outPrefix, all_residuals);

    return EXIT_SUCCESS;
}

//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv) {
	Eigen::initParallel();

    // Check for complex, do everything else inside templated function
    bool use_complex = false;
	while ((c = getopt_long(argc, argv, short_opts, long_opts, &indexptr)) != -1) {
		switch (c) {
            case 'x': use_complex = true; break;
            case 'h':
                cout << usage << endl;
                return EXIT_SUCCESS;
            case '?': // getopt will print an error message
                return EXIT_FAILURE;
			default: break;
		}
	}

    if (use_complex) {
        return run_main<complex<double>>(argc, argv);
    } else {
        return run_main<double>(argc, argv);
    }
}
