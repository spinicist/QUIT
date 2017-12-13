add_library( qi_signals STATIC
             Common.cpp SignalEquations.cpp Lineshape.cpp
             SPGR.cpp SSFP.cpp SSFP_MC.cpp MPRAGE.cpp )
target_include_directories( qi_signals PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )