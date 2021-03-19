add_executable(qi qi_main.cpp)
set_target_properties(qi PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON
    CXX_EXTENSIONS OFF)
target_link_libraries(qi PRIVATE
    taywee::args
    nlohmann_json nlohmann_json::nlohmann_json
    fmt::fmt
    ITKCommon ITKStatistics ITKTransform ITKSpatialObjects ITKLabelMap
    ITKIOImageBase ITKIONIFTI
    ITKTransformFactory ITKIOTransformBase ITKIOTransformInsightLegacy
    ceres
    Eigen3::Eigen
)
install( TARGETS qi RUNTIME DESTINATION bin )

add_subdirectory( Core )
add_subdirectory( ImageIO )
add_subdirectory( Sequences )
add_subdirectory( CoreProgs )
add_subdirectory( B1 )
add_subdirectory( MT )
add_subdirectory( Perfusion )
add_subdirectory( Relaxometry )
add_subdirectory( PARMESAN )
add_subdirectory( Stats )
add_subdirectory( Susceptibility )
add_subdirectory( Utils )
