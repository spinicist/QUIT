#include "mt_sequence.h"
#include "Log.h"

void from_json(json const &j, MTSequence &s) {
    j.at("TR").get_to(s.TR);
    j.at("Trf").get_to(s.Trf);
    j.at("Tramp").get_to(s.Tramp);
    j.at("Tspoil").get_to(s.Tspoil);
    j.at("SPS").get_to(s.SPS);
    j.at("MT_pulse").get_to(s.MT_pulse);
    j.at("MT_pulsewidth").get_to(s.MT_pulsewidth);
    s.RUFIS_FA   = QI::ArrayFromJSON(j, "RUFIS_FA", M_PI / 180.0);
    s.MT_FA      = QI::ArrayFromJSON(j, "MT_FA", M_PI / 180.0);
    s.MT_offsets = QI::ArrayFromJSON(j, "MT_offsets");
    if (s.RUFIS_FA.rows() != s.MT_FA.rows()) {
        QI::Fail("Number of RUFIS FAs {} does not match number of MT FAs {}",
                 s.RUFIS_FA.rows(),
                 s.MT_FA.rows());
    }
}

void to_json(json &j, MTSequence const &s) {
    j = json{{"TR", s.TR},
             {"Trf", s.Trf},
             {"Tramp", s.Tramp},
             {"Tspoil", s.Tspoil},
             {"SPS", s.SPS},
             {"RUFIS_FA", s.RUFIS_FA * 180. / M_PI},
             {"MT_FA", s.MT_FA * 180. / M_PI},
             {"MT_offsets", s.MT_offsets},
             {"MT_pulsewidth", s.MT_pulsewidth},
             {"MT_pulse", s.MT_pulse}};
}