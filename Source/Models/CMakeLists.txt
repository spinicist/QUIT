add_library( qi_models
             Model.cpp DESPOT_2C.cpp DESPOT_3C.cpp qMT.cpp )
target_link_libraries( qi_models qi_signals qi_core )
target_include_directories( qi_models PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )
install (TARGETS qi_models
         ARCHIVE DESTINATION lib
         LIBRARY DESTINATION lib
         RUNTIME DESTINATION bin)
