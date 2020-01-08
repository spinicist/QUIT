#include "Log.h"
#include "mupa_sequence.h"

void from_json(const json &j, MUPASequence &s) {
    j.at("TR").get_to(s.TR);
    j.at("Tramp").get_to(s.Tramp);
    s.FA  = QI::ArrayFromJSON(j, "FA", M_PI / 180.0);
    s.Trf = QI::ArrayFromJSON(j, "Trf", 1.e-6);
    j.at("SPS").get_to(s.SPS);
    j.at("prep_pulses").get_to(s.prep_pulses);
    j.at("prep").get_to(s.prep);
    if (s.FA.rows() != static_cast<long>(s.prep.size())) {
        QI::Fail(
            "Number preps {} does not match number of flip-angles {}", s.prep.size(), s.FA.rows());
    }
}