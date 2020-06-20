#include "Log.h"
#include "ss_sequence.h"

Eigen::Index SSSequence::size() const {
    return FA.rows();
}

void from_json(const json &j, SSSequence &s) {
    QI::GetJSON(j, "TR", s.TR);
    QI::GetJSON(j, "Trf", s.Trf);
    QI::GetJSON(j, "Tramp", s.Tramp);
    QI::GetJSON(j, "spokes_per_seg", s.spokes_per_seg);
    s.FA = QI::ArrayFromJSON(j, "FA", M_PI / 180.0);
    QI::GetJSON(j, "prep_p1", s.prep_p1);
    QI::GetJSON(j, "prep_p2", s.prep_p2);
    QI::GetJSON(j, "prep_Trf", s.prep_Trf);
    s.prep_FA = QI::ArrayFromJSON(j, "prep_FA", M_PI / 180.0, s.FA.rows());
    s.prep_df = QI::ArrayFromJSON(j, "prep_df", 1., s.FA.rows());
}