#include "prep_sequence.h"
#include "Log.h"

PrepSequence::PrepSequence(double         TR_,
                           double         Tramp_,
                           double         Trf_,
                           double         Dprep_,
                           double         Dseg_,
                           int            SPS_,
                           int            spoilers_,
                           Eigen::ArrayXd FA_,
                           Eigen::ArrayXd FAprep_,
                           Eigen::ArrayXd Tprep_,
                           Eigen::ArrayXd fprep_) :
    TR{TR_}, Tramp{Tramp_}, Trf{Trf_}, Dprep{Dprep_}, Dseg{Dseg_}, SPS{SPS_}, spoilers{spoilers_},
    FA{FA_}, FAprep{FAprep_}, Tprep{Tprep_}, fprep{fprep_} {
    if (FA.size() != FAprep.size()) {
        QI::Fail("Number of prep flip-angles {} does not match number of flip-angles {}",
                 FAprep.size(),
                 FA.size());
    }
    if (FA.size() != Tprep.size()) {
        QI::Fail("Number of prep pulse times {} does not match number of flip-angles {}",
                 Tprep.size(),
                 FA.size());
    }
    if (FA.size() != fprep.size()) {
        QI::Fail("Number of prep frequencies {} does not match number of flip-angles {}",
                 fprep.size(),
                 FA.size());
    }
}

auto PrepSequence::preps() const -> Eigen::Index {
    return FAprep.size();
}
auto PrepSequence::size() const -> Eigen::Index {
    return SPS * FAprep.size();
}

void from_json(const json &j, PrepSequence &s) {
    QI::GetJSON(j, "TR", s.TR);
    QI::GetJSON(j, "Trf", s.Trf);
    QI::GetJSON(j, "Tramp", s.Tramp);
    QI::GetJSON(j, "Dprep", s.Dprep);
    QI::GetJSON(j, "Dseg", s.Dseg);
    QI::GetJSON(j, "SPS", s.SPS);
    QI::GetJSON(j, "spoilers", s.spoilers);
    s.FA     = QI::ArrayFromJSON(j, "FA", (M_PI / 180.0));
    s.FAprep = QI::ArrayFromJSON(j, "FAprep", (M_PI / 180.0));
    s.Tprep  = QI::ArrayFromJSON(j, "Tprep", 1.0);
    if (j.contains("fprep")) {
        s.fprep = QI::ArrayFromJSON(j, "fprep", 1.0);
    } else {
        s.fprep = Eigen::ArrayXd::Zero(s.FAprep.size());
    }
    if (s.FA.size() != s.FAprep.size()) {
        QI::Fail("Number of prep flip-angles {} does not match number of flip-angles {}",
                 s.FAprep.size(),
                 s.FA.size());
    }
    if (s.FA.size() != s.Tprep.size()) {
        QI::Fail("Number of prep pulse times {} does not match number of flip-angles {}",
                 s.Tprep.size(),
                 s.FA.size());
    }
    if (s.FA.size() != s.fprep.size()) {
        QI::Fail("Number of prep frequencies {} does not match number of flip-angles {}",
                 s.fprep.size(),
                 s.FA.size());
    }
}