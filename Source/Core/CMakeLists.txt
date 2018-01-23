include_directories(${CMAKE_CURRENT_BINARY_DIR}) # For version file
set( VERSION_FILE_NAME "VersionFile")
add_custom_target( qi_version ALL DEPENDS version_depend_dummy )
add_custom_command(OUTPUT version_depend_dummy
                   COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_SOURCE_DIR}/CMake/GetVersion.cmake)
set_source_files_properties(${CMAKE_CURRENT_BINARY_DIR}/${VERSION_FILE_NAME} PROPERTIES GENERATED TRUE HEADER_FILE_ONLY TRUE )

set(CORE "qi_core")
add_library( ${CORE} STATIC
             Types.h Macro.h Args.h IO.h EigenCereal.h
             Util.cpp ThreadPool.cpp
             GoldenSection.cpp Masking.cpp
             Kernels.cpp Fit.cpp Spline.cpp )
add_dependencies( ${CORE} qi_version )
target_link_libraries( ${CORE} qi_filters )
target_include_directories( ${CORE} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> )

set( PROGRAMS qinewimage qisignal qidiff )
#set( PROGRAMS qitest )

foreach( PROGRAM ${PROGRAMS} )
    add_executable( ${PROGRAM} ${PROGRAM}.cpp )
    target_link_libraries( ${PROGRAM} qi_sequences qi_core ${ITK_LIBRARIES} )
endforeach( PROGRAM )

install( TARGETS ${PROGRAMS} RUNTIME DESTINATION bin )