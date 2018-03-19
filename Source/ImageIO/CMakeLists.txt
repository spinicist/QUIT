add_library( qi_imageio
             ImageRead.cpp ImageWrite.cpp
             VectorImageRead.cpp VectorImageWrite.cpp )
target_link_libraries( qi_imageio qi_filters qi_core )
target_include_directories( qi_imageio PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )
install( TARGETS qi_imageio
         ARCHIVE DESTINATION lib
         LIBRARY DESTINATION lib
         RUNTIME DESTINATION bin )
