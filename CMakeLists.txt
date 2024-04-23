cmake_minimum_required(VERSION 3.9.0)
project(GEOSHeightCorrection VERSION 1.0)

find_package(boost_program_options 1.70 REQUIRED)
find_package(dlib REQUIRED)
find_package(OpenMP)

execute_process(
    COMMAND gdal-config --cflags 
    OUTPUT_VARIABLE GDAL_FLAGS
)

add_compile_options( ${GDAL_FLAGS} )

set(CMAKE_CXX_STANDARD 17)

add_executable(geosheightcorrection main.cpp correction.cpp)
target_link_libraries(geosheightcorrection boost_program_options gdal dlib blas)

add_executable(geostablegenerator tableGenerator.cpp correction.cpp)
set(TABLE_GENERATOR_LIBS boost_program_options gdal dlib blas)

if(OpenMP_CXX_FOUND)
    list(APPEND TABLE_GENERATOR_LIBS ${OpenMP_CXX_LIBRARIES})
    set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()


target_link_libraries(geostablegenerator ${TABLE_GENERATOR_LIBS})

enable_testing()

add_subdirectory(tests)