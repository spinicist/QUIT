set( header_files
        ${CMAKE_CURRENT_SOURCE_DIR}/ImageToVectorFilter.h
        ${CMAKE_CURRENT_SOURCE_DIR}/ImageToVectorFilter.hxx
        ${CMAKE_CURRENT_SOURCE_DIR}/VectorToImageFilter.h
        ${CMAKE_CURRENT_SOURCE_DIR}/VectorToImageFilter.hxx
        ${CMAKE_CURRENT_SOURCE_DIR}/ApplyAlgorithmFilter.h
        ${CMAKE_CURRENT_SOURCE_DIR}/ApplyAlgorithmFilter.hxx
        ${CMAKE_CURRENT_SOURCE_DIR}/ReorderVectorFilter.h
        ${CMAKE_CURRENT_SOURCE_DIR}/ReorderVectorFilter.hxx
        ${CMAKE_CURRENT_SOURCE_DIR}/ReorderImageFilter.h
        ${CMAKE_CURRENT_SOURCE_DIR}/ReorderImageFilter.hxx )


add_library( qi_filters INTERFACE )
target_sources( qi_filters INTERFACE ${header_files} )
target_include_directories( qi_filters INTERFACE
                                ${CMAKE_CURRENT_SOURCE_DIR} )