set( IMAGEIO "qi_imageio")
add_library( ${IMAGEIO}
             ImageRead.cpp ImageWrite.cpp
             VectorImageRead.cpp VectorImageWrite.cpp )
target_link_libraries( ${IMAGEIO} qi_filters qi_core )
target_include_directories( ${IMAGEIO} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )
install (TARGETS ${IMAGEIO}
         ARCHIVE DESTINATION lib
         LIBRARY DESTINATION lib
         RUNTIME DESTINATION bin)
