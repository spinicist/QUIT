add_library( qi_models
             ModelBase.cpp
             OnePoolModel.cpp TwoPoolModel.cpp ThreePoolModel.cpp
             RamaniModel.cpp )
target_link_libraries( qi_models qi_signals qi_core )
target_include_directories( qi_models PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )
set_target_properties( qi_models PROPERTIES VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}
                                          SOVERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH} )
