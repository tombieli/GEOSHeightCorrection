cmake_minimum_required(VERSION 3.0.0)
project(GEOSHeightCorrection VERSION 0.1.0)

add_executable(geosheightcorrection main.cpp correction.cpp)
target_link_libraries(geosheightcorrection boost_program_options gdal)