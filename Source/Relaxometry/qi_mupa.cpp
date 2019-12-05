/*
 *  qi_vfa_prep.cpp
 *
 *  Copyright (c) 2019 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <type_traits>
#include <unsupported/Eigen/MatrixFunctions>

// #define QI_DEBUG_BUILD 1

#include "Args.h"
#include "ImageIO.h"
#include "Macro.h"
#include "Model.h"
#include "ModelFitFilter.h"
#include "SequenceBase.h"
#include "SimulateModel.h"
#include "Util.h"

struct MUPASequence : QI::SequenceBase {
    double                   TR, Tramp, FA;
    int                      SPS;
    std::vector<std::string> prep_type;
    std::vector<double>      prep_time;
    QI_SEQUENCE_DECLARE(MUPA);
    Eigen::Index size() const override { return prep_type.size(); };
};
void from_json(const json &j, MUPASequence &s) {
    j.at("TR").get_to(s.TR);
    j.at("Tramp").get_to(s.Tramp);
    s.FA = j.at("FA").get<double>() * M_PI / 180.0;
    j.at("SPS").get_to(s.SPS);
    j.at("prep_type").get_to(s.prep_type);
    j.at("prep_time").get_to(s.prep_time);
}

void to_json(json &j, const MUPASequence &s) {
    j = json{{"TR", s.TR},
             {"Tramp", s.Tramp},
             {"FA", s.FA * 180 / M_PI},
             {"SPS", s.SPS},
             {"prep_type", s.prep_type},
             {"prep_time", s.prep_time}};
}

struct MUPAModel {
    using DataType      = double;
    using ParameterType = double;

    static constexpr int NV = 3; // Number of varying parameters
    static constexpr int ND = 0; // Number of derived parameters
    static constexpr int NF = 0; // Number of fixed parameters
    static constexpr int NI = 1; // Number of inputs

    using VaryingArray = QI_ARRAYN(ParameterType, NV); // Type for the varying parameter array
    using FixedArray   = QI_ARRAYN(ParameterType, NF); // Type for the fixed parameter array

    // Sequence paramter structs
    MUPASequence &sequence;

    // Fitting start point and bounds
    // The PD will be scaled by the fitting function to keep parameters roughly the same magnitude
    VaryingArray const start{30., 1., 0.1};
    VaryingArray const bounds_lo{1, 0.01, 0.01};
    VaryingArray const bounds_hi{150, 10.0, 10.0};

    std::array<std::string, NV> const varying_names{"PD", "T1", "T2"};
    std::array<std::string, NF> const fixed_names{};
    // If fixed parameters not supplied, use these default values
    FixedArray const fixed_defaults{};

    template <typename Derived>
    auto signal(Eigen::ArrayBase<Derived> const &v, FixedArray const &) const
        -> QI_ARRAY(typename Derived::Scalar) {
        using T    = typename Derived::Scalar;
        using Vec3 = Eigen::Vector<T, 3>;
        // using Vec4  = Eigen::Vector<T, 4>;
        using Mat33 = Eigen::Matrix<T, 3, 3>;
        using Mat44 = Eigen::Matrix<T, 4, 4>;
        T const &PD = v[0];
        T const &R1 = 1. / v[1];
        T const &R2 = 1. / v[2];

        auto const Relax = [&PD, &R1, &R2](double const t) {
            Mat44 R;
            R << -R2, 0, 0, 0,      //
                0, -R2, 0, 0,       //
                0, 0, -R1, PD * R1, //
                0, 0, 0, 0;
            Mat44 eRt = (R * t).exp();
            return eRt;
        };

        auto const RF = [](double const a, double const ph) {
            auto const ca = cos(a);
            auto const sa = sin(a);
            // I got the definition of the rotation matrix around the wrong axis,
            // so rotate the axis
            auto const      ux = cos(ph - M_PI_2);
            auto const      uy = sin(ph - M_PI_2);
            Eigen::Matrix4d A;
            A << ca + ux * ux * (1 - ca), ux * uy * (1 - ca), -uy * sa, 0., //
                ux * uy * (1 - ca), ca + uy * uy * (1 - ca), ux * sa, 0.,   //
                uy * sa, -ux * sa, ca, 0.,                                  //
                0., 0., 0., 1.;
            return A;
        };

        auto const  A     = RF(sequence.FA, 0);
        auto const  R     = Relax(sequence.TR);
        auto const  S     = Eigen::DiagonalMatrix<double, 4, 4>({0, 0, 1., 1.}).toDenseMatrix();
        Mat44 const ramp  = Relax(sequence.Tramp);
        Mat44 const RUFIS = S * R * A;
        Mat44 const seg   = RUFIS.pow(sequence.SPS);
        // First calculate the system matrix
        Mat44 X = Mat44::Identity();
        for (int is = 0; is < sequence.size(); is++) {
            double const tprep = sequence.prep_time[is];
            if (sequence.prep_type[is] == "inversion") {
                auto const I  = RF(M_PI, 0);
                auto const TI = Relax(tprep);
                X             = S * TI * I * X;
            } else if (sequence.prep_type[is] == "t2-prep") {
                auto const TD = RF(M_PI_2, 0);
                auto const TE = Relax(tprep);
                auto const TU = RF(-M_PI_2, 0);
                X             = S * TU * TE * TD * X;
            } else if (sequence.prep_type[is] == "delay") {
                auto const D = Relax(tprep);
                X            = D * X;
            }
            X = ramp * seg * ramp * X;
        }

        // Solve for steady-state and re-augment
        Mat33           Xr    = (X - Mat44::Identity()).template topLeftCorner<3, 3>();
        Eigen::Vector3d b     = -X.template topRightCorner<3, 1>();
        Eigen::Vector3d m_ss  = Xr.partialPivLu().solve(b);
        Eigen::Vector4d m_aug = Eigen::Vector4d::Ones();
        m_aug.head(3)         = m_ss;
        QI_DBVEC(m_ss);
        // Now loop through the segments and record the signal for each
        Eigen::ArrayXd sig(sequence.size());
        for (int is = 0; is < sequence.size(); is++) {
            double const tprep = sequence.prep_time[is];
            if (sequence.prep_type[is] == "inversion") {
                auto const I  = RF(M_PI, 0);
                auto const TI = Relax(tprep);
                m_aug         = S * TI * I * m_aug;
            } else if (sequence.prep_type[is] == "t2-prep") {
                auto const TD = RF(M_PI / 2, 0);
                auto const TE = Relax(tprep);
                auto const TU = RF(-M_PI / 2, 0);
                m_aug         = S * TU * TE * TD * m_aug;
            } else if (sequence.prep_type[is] == "delay") {
                auto const D = Relax(tprep);
                m_aug        = D * m_aug;
            }
            m_aug = ramp * m_aug;
            // Geometric Sum formula to get average
            Mat44 const LHS = (Mat44::Identity() - RUFIS);
            QI_DBMAT(LHS);
            Vec3 const RHS = ((Mat44::Identity() - seg) * m_aug).template head<3>() -
                             (sequence.SPS * LHS.template topRightCorner<3, 1>());
            QI_DBVEC(RHS);
            Vec3 const m_gm =
                LHS.template topLeftCorner<3, 3>().partialPivLu().solve(RHS) / sequence.SPS;
            QI_DBVEC(m_gm);
            sig[is] = m_gm[2] * sin(sequence.FA);
            QI_DB(sig[is]);
            m_aug = ramp * seg * m_aug;
        }
        QI_DB(PD);
        QI_DB(1 / R1);
        QI_DB(1 / R2);
        QI_DBVEC(sig);
        return sig;
    }
};

template <> struct QI::NoiseFromModelType<MUPAModel> : QI::RealNoise {};

struct MUPACost {
    MUPAModel const &     model;
    MUPAModel::FixedArray fixed;
    QI_ARRAY(double) const data;

    template <typename T> bool operator()(T const *const vin, T *rin) const {
        Eigen::Map<QI_ARRAYN(T, MUPAModel::NV) const> const varying(vin);
        Eigen::Map<QI_ARRAY(T)>                             residuals(rin, data.rows());
        residuals = data - model.signal(varying, fixed);
        QI_DBVEC(residuals);
        return true;
    }
};

struct MUPAFit {
    // Boilerplate information required by ModelFitFilter
    static const bool Blocked = false; // = input is in blocks and outputs have multiple entries
    static const bool Indexed = false; // = the voxel index will be passed to the fit
    using RMSErrorType        = double;
    using FlagType            = int; // Almost always the number of iterations

    using ModelType = MUPAModel;
    ModelType model;

    // Have to tell the ModelFitFilter how many volumes we expect in each input
    int input_size(const int) const { return model.sequence.size(); }

    // This has to match the function signature that will be called in ModelFitFilter (which depends
    // on Blocked/Indexed. The return type is a simple struct indicating success, and on failure
    // also the reason for failure
    QI::FitReturnType
    fit(std::vector<Eigen::ArrayXd> const &inputs,    // Input: signal data
        Eigen::ArrayXd const &             fixed,     // Input: Fixed parameters
        ModelType::VaryingArray &          varying,   // Output: Varying parameters
        RMSErrorType &                     rmse,      // Output: root-mean-square error
        std::vector<Eigen::ArrayXd> &      residuals, // Optional output: point residuals
        FlagType &                         iterations /* Usually iterations */) const {
        // First scale down the raw data so that PD will be roughly the same magnitude as other
        // parameters This is important for numerical stability in the optimiser

        QI_DBVEC(inputs[0]);

        double scale = inputs[0].maxCoeff();
        QI_DB(scale);
        if (scale < std::numeric_limits<double>::epsilon()) {
            varying = ModelType::VaryingArray::Zero();
            rmse    = 0.0;
            return {false, "Maximum data value was zero or less"};
        }
        Eigen::ArrayXd const data = inputs[0] / scale;

        // Setup Ceres
        ceres::Problem problem;
        using VFADiff =
            ceres::NumericDiffCostFunction<MUPACost, ceres::CENTRAL, ceres::DYNAMIC, ModelType::NV>;
        auto *vfa_prep_cost = new VFADiff(
            new MUPACost{model, fixed, data}, ceres::TAKE_OWNERSHIP, model.sequence.size());
        ceres::LossFunction *loss = new ceres::HuberLoss(1.0); // Don't know if this helps
        // This is where the parameters and cost functions actually get added to Ceres
        problem.AddResidualBlock(vfa_prep_cost, loss, varying.data());

        // Set up parameter bounds
        for (int i = 0; i < ModelType::NV; i++) {
            problem.SetParameterLowerBound(varying.data(), i, model.bounds_lo[i]);
            problem.SetParameterUpperBound(varying.data(), i, model.bounds_hi[i]);
        }

        ceres::Solver::Options options;
        ceres::Solver::Summary summary;
        options.max_num_iterations  = 30;
        options.function_tolerance  = 1e-5;
        options.gradient_tolerance  = 1e-6;
        options.parameter_tolerance = 1e-4;
        options.logging_type        = ceres::SILENT;

        varying = model.start;
        ceres::Solve(options, &problem, &summary);
        if (!summary.IsSolutionUsable()) {
            return {false, summary.FullReport()};
        }

        Eigen::ArrayXd const residual = (data - model.signal(varying, fixed)) * scale;
        rmse                          = sqrt(residual.square().mean());
        if (residuals.size() > 0) {
            residuals[0] = residual;
        }
        varying[0] *= scale; // Multiply signals/proton density back up
        QI_DBVEC(varying);
        iterations = summary.iterations.size();
        return {true, ""};
    }
};

/*
 * Main
 */
int mupa_main(int argc, char **argv) {
    Eigen::initParallel();
    args::ArgumentParser          parser("Calculates PD/T1/T2 from MUPA data "
                                "data.\nhttp://github.com/spinicist/QUIT");
    args::Positional<std::string> input_path(parser, "INPUT", "Input MUPA file");

    QI_COMMON_ARGS;

    QI::ParseArgs(parser, argc, argv, verbose, threads);

    QI::CheckPos(input_path);

    QI::Log(verbose, "Reading sequence parameters");
    json doc = json_file ? QI::ReadJSON(json_file.Get()) : QI::ReadJSON(std::cin);

    MUPASequence sequence(doc["MUPA"]);
    MUPAModel    model{sequence};
    if (simulate) {
        QI::SimulateModel<MUPAModel, false>(
            doc, model, {}, {input_path.Get()}, verbose, simulate.Get());
    } else {
        MUPAFit fit{model};
        auto fit_filter = QI::ModelFitFilter<MUPAFit>::New(&fit, verbose, resids, subregion.Get());
        fit_filter->ReadInputs({input_path.Get()}, {}, mask.Get());
        fit_filter->Update();
        fit_filter->WriteOutputs(prefix.Get() + "MUPA_");
    }
    QI::Log(verbose, "Finished.");
    return EXIT_SUCCESS;
}
