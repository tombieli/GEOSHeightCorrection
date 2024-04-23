enable_testing()

find_package(boost_unit_test_framework 1.70 REQUIRED)

add_executable(tests main_test.cpp cloud_simulation.cpp fixtures.cpp correction_tests.cpp pixel_correction_tests.cpp ../correction.cpp)
include_directories( ${CMAKE_CURRENT_LIST_DIR}/.. )
# target_compile_features(tests PRIVATE cxx_std_17)
target_link_libraries(tests boost_unit_test_framework gdal dlib blas)

add_test(NAME tests
         COMMAND tests)
