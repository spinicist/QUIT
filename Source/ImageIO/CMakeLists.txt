add_library( qi_imageio
             ImageRead.cpp ImageWrite.cpp
             VectorImageRead.cpp VectorImageWrite.cpp )
target_link_libraries( qi_imageio PRIVATE qi_core ${ITK_LIBRARIES} )
target_include_directories( qi_imageio PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )
set_target_properties( qi_imageio PROPERTIES VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}
                                           SOVERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH} )