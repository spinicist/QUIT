set( FILTERS "qi_filters")
add_library( ${FILTERS}
             ImageToVectorFilter.h
             VectorToImageFilter.h
             ApplyAlgorithmFilter.h
             ReorderVectorFilter.h
             ReorderImageFilter.h
             ApplyTypes.h
             VolumeFilters.cpp VectorVolumeFilters.cpp )
target_link_libraries( ${FILTERS} qi_core ${ITK_LIBRARIES} )
target_include_directories( ${FILTERS} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )
install (TARGETS ${FILTERS}
         ARCHIVE DESTINATION lib
         LIBRARY DESTINATION lib
         RUNTIME DESTINATION bin)