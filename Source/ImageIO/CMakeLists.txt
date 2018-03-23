add_library( qi_imageio
             ImageRead.cpp ImageWrite.cpp
             VectorImageRead.cpp VectorImageWrite.cpp )
target_link_libraries( qi_imageio PRIVATE qi_filters qi_core ${ITK_LIBRARIES} )
target_include_directories( qi_imageio PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )
install( TARGETS qi_imageio
         ARCHIVE DESTINATION lib
         LIBRARY DESTINATION lib
         RUNTIME DESTINATION bin )
set_target_properties( qi_imageio PROPERTIES VERSION 2.0.0
                                           SOVERSION 2.0.0 )