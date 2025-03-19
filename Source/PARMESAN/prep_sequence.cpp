#include "prep_sequence.h"
#include "Log.h"

void from_json(const json &j, PrepZTESequence &s) {
    QI::GetJSON(j, "TR", s.TR);
    QI::GetJSON(j, "Trf", s.Trf);
    QI::GetJSON(j, "Tprep", s.Tprep);
    QI::GetJSON(j, "Tramp", s.Tramp);
    QI::GetJSON(j, "SPS", s.SPS);
    QI::GetJSON(j, "spoilers", s.spoilers);
    s.FA     = QI::ArrayFromJSON(j, "FA", M_PI / 180.0);
    s.FAprep = QI::ArrayFromJSON(j, "FAprep", M_PI / 180.0);
    if (j.contains("basis")) {
        s.basis = QI::MatrixFromJSON(j, "basis", 1.0, -1, s.SPS * s.FAprep.size());
    }
    if (s.FA.rows() != static_cast<long>(s.FAprep.size())) {
        QI::Fail("Number preps {} does not match number of flip-angles {}",
                 s.FAprep.size(),
                 s.FA.rows());
    }
}