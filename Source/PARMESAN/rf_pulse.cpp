#include "JSON.h"
#include "Log.h"
#include "rf_pulse.h"

void from_json(const json &j, RFPulse &p) {
    p.B1x      = QI::ArrayFromJSON(j, "B1x", 1.);
    p.B1y      = QI::ArrayFromJSON(j, "B1y", 1.);
    p.timestep = QI::ArrayFromJSON(j, "timestep", 1.);

    if (p.B1x.rows() != p.B1y.rows()) {
        QI::Fail("B1x and B1y array lengths did not match {} vs {}", p.B1x.rows(), p.B1y.rows());
    }
}

void to_json(json &j, RFPulse const &p) {
    j = json{{"B1x", p.B1x}, {"B1y", p.B1y}, {"timestep", p.timestep}};
}

void from_json(const json &j, MTPulse &p) {
    p.B1x = QI::ArrayFromJSON(j, "B1x", 1.);
    p.B1y = QI::ArrayFromJSON(j, "B1y", 1.);
    p.FA  = j.at("FA").get<double>() * M_PI / 180;
    j.at("width").get_to(p.width);

    if (p.B1x.rows() != p.B1y.rows()) {
        QI::Fail("B1x and B1y array lengths did not match {} vs {}", p.B1x.rows(), p.B1y.rows());
    }
}

void to_json(json &j, MTPulse const &p) {
    j = json{{"B1x", p.B1x}, {"B1y", p.B1y}, {"FA", p.FA * 180 / M_PI}, {"width", p.width}};
}

void from_json(const json &j, PrepPulse &p) {
    p.FAeff = j.at("FAeff").get<double>() * M_PI / 180.0;
    j.at("int_b1_sq").get_to(p.int_b1_sq);
    j.at("T_long").get_to(p.T_long);
    j.at("T_trans").get_to(p.T_trans);
}

void to_json(json &j, PrepPulse const &p) {
    j = json{{"FAeff", p.FAeff * 180 / M_PI},
             {"int_b1_sq", p.int_b1_sq},
             {"T_long", p.T_long},
             {"T_trans", p.T_trans}};
}