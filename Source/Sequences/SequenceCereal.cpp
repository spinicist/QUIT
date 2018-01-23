
#include "SequenceCereal.h"
#include "SPGRSequence.h"
#include "SSFPSequence.h"

namespace cereal {

void save(cereal::JSONOutputArchive &ar, std::shared_ptr<QI::SequenceBase> const &s) {
    if (s->name() == "SPGR") { ar(cereal::make_nvp("SPGR", *std::static_pointer_cast<QI::SPGRSequence>(s))); }
    else if (s->name() == "SPGREcho") { ar(cereal::make_nvp("SPGREcho", *std::static_pointer_cast<QI::SPGREchoSequence>(s))); }
    else if (s->name() == "SPGRFiniteEcho") { ar(cereal::make_nvp("SPGRFinite", *std::static_pointer_cast<QI::SPGRFiniteSequence>(s))); }
    else { QI_FAIL("Unimplemented read for sequence type: " << s->name()); }
}

void load(cereal::JSONInputArchive &ar, std::shared_ptr<QI::SequenceBase> &sb) {
    std::string seq_type = ar.getNodeName();
    if (seq_type == "SPGR") { QI::SPGRSequence s; ar(s); sb = std::make_shared<QI::SPGRSequence>(s); }
    else if (seq_type == "SPGREcho") { QI::SPGREchoSequence s; ar(s); sb = std::make_shared<QI::SPGREchoSequence>(s); }
    else if (seq_type == "SPGRFiniteEcho") { QI::SPGRFiniteSequence s; ar(s); sb = std::make_shared<QI::SPGRFiniteSequence>(s); }
    else { QI_FAIL("Unimplemented load for sequence type: " << seq_type); }
}

} // End namespace cereal