add_library( qi_signals
             Common.cpp SignalEquations.cpp Lineshape.cpp
             SPGR.cpp SSFP.cpp SSFP_MC.cpp MPRAGE.cpp )
target_include_directories( qi_signals PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )
set_target_properties( qi_signals PROPERTIES VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}
                                           SOVERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH} )