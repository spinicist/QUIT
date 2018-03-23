add_library( qi_signals
             Common.cpp SignalEquations.cpp Lineshape.cpp
             SPGR.cpp SSFP.cpp SSFP_MC.cpp MPRAGE.cpp )
target_include_directories( qi_signals PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )
install( TARGETS qi_signals
         ARCHIVE DESTINATION lib
         LIBRARY DESTINATION lib
         RUNTIME DESTINATION bin )
set_target_properties( qi_signals PROPERTIES VERSION 2.0.0
                                           SOVERSION 2.0.0 )