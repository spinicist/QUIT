add_library( qi_models
             Model.cpp DESPOT_2C.cpp DESPOT_3C.cpp qMT.cpp )
target_link_libraries( qi_models qi_signals qi_core )
target_include_directories( qi_models PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )
set_target_properties( qi_models PROPERTIES VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}
                                          SOVERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH} )
