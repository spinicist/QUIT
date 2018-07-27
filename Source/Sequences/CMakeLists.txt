add_library( qi_sequences
                SequenceBase.cpp
                AFISequence.cpp
                CASLSequence.cpp
                MPRAGESequence.cpp
                MTSatSequence.cpp
                MultiEchoSequence.cpp
                SPGRSequence.cpp
                SSFPSequence.cpp
                SequenceGroup.cpp
           )
target_link_libraries( qi_sequences qi_core )
target_include_directories( qi_sequences PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )
set_target_properties( qi_sequences PROPERTIES VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}
                                             SOVERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH} )