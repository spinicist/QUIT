#include "prep_sequence.h"
#include "Log.h"

PrepSequence::PrepSequence(double         TR_,
                           double         Tramp_,
                           double         Trf_,
                           int            SPS_,
                           int            spoilers_,
                           Eigen::ArrayXd FA_,
                           Eigen::ArrayXd FAprep_,
                           Eigen::ArrayXd Tprep_,
                           Eigen::ArrayXd Tpreseg_,
                           Eigen::ArrayXd Tpostseg_,
                           Eigen::ArrayXd fprep_) :
    TR{TR_}, Tramp{Tramp_}, Trf{Trf_}, SPS{SPS_}, spoilers{spoilers_}, FA{FA_ * M_PI / 180.},
    FAprep{FAprep_ * M_PI / 180.}, Tprep{Tprep_}, Tpreseg{Tpreseg_}, Tpostseg{Tpostseg_},
    fprep{fprep_} {
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
    QI::GetJSON(j, "SPS", s.SPS);
    QI::GetJSON(j, "spoilers", s.spoilers);
    s.FA       = QI::ArrayFromJSON(j, "FA", (M_PI / 180.0));
    s.FAprep   = QI::ArrayFromJSON(j, "FAprep", (M_PI / 180.0));
    s.Tprep    = QI::ArrayFromJSON(j, "Tprep", 1.);
    s.Tpreseg  = QI::ArrayFromJSON(j, "Tpreseg", 1.);
    s.Tpostseg = QI::ArrayFromJSON(j, "Tpostseg", 1.);
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
    if (s.FA.size() != s.Tpreseg.size()) {
        QI::Fail("Number of pre-segment times {} does not match number of flip-angles {}",
                 s.Tpreseg.size(),
                 s.FA.size());
    }
    if (s.FA.size() != s.Tpostseg.size()) {
        QI::Fail("Number of post-segment pulse times {} does not match number of flip-angles {}",
                 s.Tpostseg.size(),
                 s.FA.size());
    }
    if (s.FA.size() != s.fprep.size()) {
        QI::Fail("Number of prep frequencies {} does not match number of flip-angles {}",
                 s.fprep.size(),
                 s.FA.size());
    }
}

void to_json(json &j, const PrepSequence &s) {
    j = json{{"TR", s.TR},
             {"Trf", s.Trf},
             {"Tramp", s.Tramp},
             {"Tpreseg", s.Tpreseg},
             {"Tpostseg", s.Tpostseg},
             {"SPS", s.SPS},
             {"spoilers", s.spoilers},
             {"FA", s.FA * 180 / M_PI},
             {"FAprep", s.FAprep * 180 / M_PI},
             {"Tprep", s.Tprep},
             {"fprep", s.fprep}};
}