file(GLOB SOURCES *.cpp)
target_sources(qi PRIVATE ${SOURCES})
target_include_directories(qi PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
