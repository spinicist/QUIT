include_directories(${CMAKE_CURRENT_BINARY_DIR}) # For version file
set( VERSION_FILE_NAME "VersionFile")
add_custom_target( qi_version ALL DEPENDS version_depend_dummy )
add_custom_command(OUTPUT version_depend_dummy
                   COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_SOURCE_DIR}/CMake/GetVersion.cmake)
set_source_files_properties(${CMAKE_CURRENT_BINARY_DIR}/${VERSION_FILE_NAME} PROPERTIES GENERATED TRUE HEADER_FILE_ONLY TRUE )

add_library( qi_core
             Macro.h Args.h IO.h EigenCereal.h ImageTypes.h
             Util.cpp ThreadPool.cpp
             GoldenSection.cpp Masking.cpp
             Kernels.cpp Fit.cpp Spline.cpp )
add_dependencies( qi_core qi_version )
target_include_directories( qi_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )
target_link_libraries( qi_core PRIVATE ${ITK_LIBRARIES} )
install( TARGETS  qi_core
         ARCHIVE DESTINATION lib
         LIBRARY DESTINATION lib
         RUNTIME DESTINATION bin )
set_target_properties( qi_core PROPERTIES VERSION 2.0.0
                                        SOVERSION 2.0.0 )