option( BUILD_RELAX "Build the relaxometry (DESPOT etc.) programs" ON )
if( ${BUILD_RELAX} )
    file(GLOB SOURCES *.cpp)
    target_sources(qi PRIVATE ${SOURCES})
    target_compile_definitions(qi PRIVATE "-DBUILD_RELAX")
endif()