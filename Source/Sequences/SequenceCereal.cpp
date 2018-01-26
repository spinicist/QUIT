
#include "SequenceCereal.h"
#include "SPGRSequence.h"
#include "SSFPSequence.h"
#include "MultiEchoSequence.h"

namespace cereal {

void save(cereal::JSONOutputArchive &ar, std::shared_ptr<QI::SequenceBase> const &s) {
    #define QI_SAVE( NAME ) \
        (s->name() == #NAME ) { ar(cereal::make_nvp(s->name(), *std::static_pointer_cast< QI::NAME ## Sequence >(s))); }
    if QI_SAVE( SPGR )
    else if QI_SAVE( SPGREcho )
    else if QI_SAVE( SPGRFinite )
    else if QI_SAVE( SSFP )
    else if QI_SAVE( SSFPEcho )
    else if QI_SAVE( SSFPFinite )
    else if QI_SAVE( SSFPGS )
    else if QI_SAVE( MultiEcho )
    else { QI_FAIL("Unimplemented save for sequence type: " << s->name()); }
    #undef QI_SAVE
}

void load(cereal::JSONInputArchive &ar, std::shared_ptr<QI::SequenceBase> &sb) {
    std::string seq_type = ar.getNodeName();
    #define QI_LOAD( NAME ) \
        (seq_type == #NAME ) { QI::NAME ## Sequence s; ar(s); sb = std::make_shared< QI::NAME ## Sequence >(s); }
    if QI_LOAD( SPGR )
    else if QI_LOAD( SPGREcho )
    else if QI_LOAD( SPGRFinite )
    else if QI_LOAD( SSFP )
    else if QI_LOAD( SSFPEcho )
    else if QI_LOAD( SSFPFinite )
    else if QI_LOAD( SSFPGS )
    else if QI_LOAD( MultiEcho )
    else { QI_FAIL("Unimplemented load for sequence type: " << seq_type); }
    #undef QI_LOAD
}

} // End namespace cereal