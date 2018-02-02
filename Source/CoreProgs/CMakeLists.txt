set( PROGRAMS qinewimage qisignal qidiff qitest )

foreach( PROGRAM ${PROGRAMS} )
    add_executable( ${PROGRAM} ${PROGRAM}.cpp )
    target_link_libraries( ${PROGRAM} qi_imageio qi_core ${ITK_LIBRARIES} )
endforeach( PROGRAM )
target_link_libraries( qisignal qi_filters qi_sequences )
target_link_libraries( qitest   qi_sequences )

install( TARGETS ${PROGRAMS} RUNTIME DESTINATION bin )