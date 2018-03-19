add_library( qi_sequences
                SequenceBase.cpp
                SPGRSequence.cpp SSFPSequence.cpp AFISequence.cpp
                MPRAGESequence.cpp MultiEchoSequence.cpp CASLSequence.cpp
                SequenceGroup.cpp SequenceCereal.cpp )
target_link_libraries( qi_sequences qi_models qi_core )
target_include_directories( qi_sequences PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )
install( TARGETS qi_sequences
         ARCHIVE DESTINATION lib
         LIBRARY DESTINATION lib
         RUNTIME DESTINATION bin )