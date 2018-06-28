
#include "SequenceCereal.h"
#include "SequenceCereal.hpp"
#include "SPGRSequence.h"
#include "MPRAGESequence.h"
#include "MTSatSequence.h"
#include "SSFPSequence.h"
#include "MultiEchoSequence.h"
#include "CASLSequence.h"
#include "SequenceGroup.h"

namespace cereal {

void save(cereal::JSONOutputArchive &ar, std::shared_ptr<QI::SequenceBase> const &s) {
    #define QI_SAVE( NAME ) \
        (s->name() == #NAME ) { ar(cereal::make_nvp(s->name(), *std::static_pointer_cast< QI::NAME ## Sequence >(s))); }
    if QI_SAVE( SPGR )
    else if QI_SAVE( SPGREcho )
    else if QI_SAVE( SPGRFinite )
    else if QI_SAVE( MPRAGE )
    else if QI_SAVE( SSFP )
    else if QI_SAVE( SSFPEcho )
    else if QI_SAVE( SSFPFinite )
    else if QI_SAVE( SSFPGS )
    else if QI_SAVE( SSFPMT )
    else if QI_SAVE( MTSat )
    else if QI_SAVE( MultiEcho )
    else if QI_SAVE( MultiEchoFlex )
    else if QI_SAVE( CASL )
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
    else if QI_LOAD( MPRAGE )
    else if QI_LOAD( SSFP )
    else if QI_LOAD( SSFPEcho )
    else if QI_LOAD( SSFPFinite )
    else if QI_LOAD( SSFPGS )
    else if QI_LOAD( SSFPMT )
    else if QI_LOAD( MTSat )
    else if QI_LOAD( MultiEcho )
    else if QI_LOAD( MultiEchoFlex )
    else if QI_LOAD( CASL )
    else { QI_FAIL("Unimplemented load for sequence type: " << seq_type); }
    #undef QI_LOAD
}

} // End namespace cereal

namespace QI {
    #define QI_READSEQ( Seq ) \
        template auto ReadSequence< Seq >(cereal::JSONInputArchive &in_archive, bool verbose, std::ostream &os = std::cout) -> Seq;\
        template auto ReadSequence< Seq >(std::istream &is, bool verbose, std::ostream &os = std::cout) -> Seq;

    QI_READSEQ( SPGRSequence )
    QI_READSEQ( SPGREchoSequence )
    QI_READSEQ( SPGRFiniteSequence )
    QI_READSEQ( MPRAGESequence )
    QI_READSEQ( SSFPSequence )
    QI_READSEQ( SSFPEchoSequence )
    QI_READSEQ( SSFPFiniteSequence )
    QI_READSEQ( SSFPGSSequence )
    QI_READSEQ( SSFPEllipseSequence )
    QI_READSEQ( SSFPMTSequence )
    QI_READSEQ( MultiEchoSequence )
    QI_READSEQ( MultiEchoFlexSequence )
    QI_READSEQ( MP2RAGESequence )
    QI_READSEQ( CASLSequence )
    QI_READSEQ( SequenceGroup )
    #undef QI_READSEQ
}