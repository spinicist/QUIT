set( FILTERS "qi_filters")
add_library( ${FILTERS} STATIC
             ImageToVectorFilter.h
             VectorToImageFilter.h
             ApplyAlgorithmFilter.h
             ReorderVectorFilter.h
             ReorderImageFilter.h
             ApplyTypes.h
             VolumeFilters.cpp VectorVolumeFilters.cpp )
target_link_libraries( ${FILTERS} qi_core ${ITK_LIBRARIES} )
target_include_directories( ${FILTERS} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )