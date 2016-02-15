# Set up the version file
FIND_PACKAGE(Git)
SET( BUILD_VERSION 1.1 )
IF( GIT_FOUND )
    EXECUTE_PROCESS(COMMAND ${GIT_EXECUTABLE} describe --tags RESULT_VARIABLE GIT_RES OUTPUT_VARIABLE GIT_DESCRIPTION )
    MESSAGE( STATUS ${GIT_RES} ${GIT_DESCRIPTION} )
    IF( NOT ${GIT_RES} EQUAL 0 )
        MESSAGE( WARNING "Git failed (not a repo, or no tags). Defaulting to static build number." )
    ENDIF()
    STRING( REPLACE "\n" "" BUILD_VERSION ${GIT_DESCRIPTION} )
    MESSAGE( STATUS ${BUILD_VERSION} ${GIT_DESCRIPTION} )
ENDIF( GIT_FOUND )
MESSAGE( STATUS "Version: ${BUILD_VERSION}" )

SET( VERSION_STRING "//version_string.hpp - written by cmake. changes will be lost!\n"
             "const std::string Version = \"${BUILD_VERSION}\"\;\n")

FILE(WRITE GetVersion.h.txt ${VERSION_STRING} )
# copy the file to the final header only if the version changes
# reduces needless rebuilds
execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different
                        GetVersion.h.txt ${CMAKE_CURRENT_BINARY_DIR}/BuildVersion.h)