add_executable(qi qi_main.cpp)
# set(ASAN $<$<CONFIG:RelWithDebInfo>:-fsanitize=address>)
target_compile_options(
    qi PRIVATE
    "-Werror" "-Wall" "-Wextra" "-Wshadow"
    "-ftemplate-backtrace-limit=30"
    ${ASAN}
)
target_link_libraries(qi ${ITK_LIBRARIES} ${CERES_LIBRARIES} fmtlib ${ASAN})
install( TARGETS qi RUNTIME DESTINATION bin )

add_subdirectory( Core )
add_subdirectory( ImageIO )
add_subdirectory( Sequences )
add_subdirectory( CoreProgs )
add_subdirectory( B1 )
add_subdirectory( MT )
add_subdirectory( Perfusion )
add_subdirectory( Relaxometry )
add_subdirectory( RUFIS )
add_subdirectory( Stats )
add_subdirectory( Susceptibility )
add_subdirectory( Utils )
