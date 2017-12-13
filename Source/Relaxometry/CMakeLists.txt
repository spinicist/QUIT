set( PROGRAMS
    qiafi qidream
    qidespot1 qidespot2 qidespot1hifi qidespot2fm qimcdespot 
    qimultiecho qimp2rage )

foreach(PROGRAM ${PROGRAMS})
    add_executable(${PROGRAM} ${PROGRAM}.cpp)
    target_link_libraries(${PROGRAM} qi_sequences qi_core ${ITK_LIBRARIES} ${CERES_LIBRARIES})
endforeach(PROGRAM)

install( TARGETS ${PROGRAMS} RUNTIME DESTINATION bin )