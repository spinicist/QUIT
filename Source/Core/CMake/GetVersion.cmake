# Set up the version file
set( VERSION_FILE_NAME "VersionFile")

find_package(Git)
set( BUILD_VERSION "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}" )
if( GIT_FOUND )
    execute_process(COMMAND ${GIT_EXECUTABLE} describe --tags --dirty RESULT_VARIABLE GIT_RES OUTPUT_VARIABLE GIT_DESCRIPTION )
    if( NOT ${GIT_RES} EQUAL 0 )
        message( WARNING "Git failed (not a repo, or no tags). Defaulting to static build number." )
    else( NOT ${GIT_RES} EQUAL 0 )
        string( REPLACE "\n" "" BUILD_VERSION ${GIT_DESCRIPTION} )
    endif( NOT ${GIT_RES} EQUAL 0 )
endif( GIT_FOUND )
message( STATUS "Version: ${BUILD_VERSION}" )

set( VERSION_STRING "const static std::string Version = \"${BUILD_VERSION}\"\;\n")

file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/${VERSION_FILE_NAME}.tmp ${VERSION_STRING} )
execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_BINARY_DIR}/${VERSION_FILE_NAME}.tmp ${CMAKE_CURRENT_BINARY_DIR}/${VERSION_FILE_NAME} )