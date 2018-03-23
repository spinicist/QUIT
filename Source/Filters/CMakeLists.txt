add_library( qi_filters
             ImageToVectorFilter.h VectorToImageFilter.h
             ApplyAlgorithmFilter.h ApplyTypes.h
             VolumeFilters.cpp VectorVolumeFilters.cpp )
target_link_libraries( qi_filters PRIVATE qi_core ${ITK_LIBRARIES} )
target_include_directories( qi_filters PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )
install( TARGETS qi_filters
         ARCHIVE DESTINATION lib
         LIBRARY DESTINATION lib
         RUNTIME DESTINATION bin )
set_target_properties( qi_filters PROPERTIES VERSION 2.0.0
                                           SOVERSION 2.0.0 )
