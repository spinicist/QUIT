add_library( qi_sequences STATIC
                SequenceBase.cpp SequenceGroup.cpp
                SPGR.cpp SSFP.cpp AFI.cpp
                MPRAGESequence.cpp SpinEcho.cpp )
target_link_libraries( qi_sequences qi_models qi_core )
target_include_directories( qi_sequences PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )
