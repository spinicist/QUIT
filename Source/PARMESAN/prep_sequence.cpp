#include "prep_sequence.h"
#include "Log.h"

auto PrepSequence::preps() const -> Eigen::Index
{ return FAprep.size(); }
auto PrepSequence::size() const -> Eigen::Index
{ return SPS * FAprep.size(); }

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
    if (s.FA.rows() != static_cast<long>(s.FAprep.size())) {
        QI::Fail("Number of prep flip-angles {} does not match number of flip-angles {}",
                 s.FAprep.size(),
                 s.FA.rows());
    }
    if (s.FA.rows() != static_cast<long>(s.Tprep.size())) {
        QI::Fail("Number of prep pulse times {} does not match number of flip-angles {}",
                 s.Tprep.size(),
                 s.FA.rows());
    }
    if (s.FA.rows() != static_cast<long>(s.fprep.size())) {
        QI::Fail("Number of prep frequencies {} does not match number of flip-angles {}",
                 s.fprep.size(),
                 s.FA.rows());
    }
}