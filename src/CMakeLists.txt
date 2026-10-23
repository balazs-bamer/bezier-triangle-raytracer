cmake_minimum_required(VERSION 3.16)
project(BezierTrangleReference)

enable_testing()
find_package(GTest REQUIRED)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type (Debug or Release)" FORCE)
endif()

add_compile_options(-std=c++17)
set(CMAKE_CXX_FLAGS_DEBUG "-ggdb -D_GLIBCXX_DEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "-O2")

SET(PreprocessBezierSrc mesh.cpp bezierTriangle.cpp bezierMesh.cpp bezierLens.cpp hostUtil.cpp)

INCLUDE_DIRECTORIES ( "/opt/alk/eigen-5.0.1"
                      "../repos/stl_reader" )
ADD_LIBRARY (PreprocessBezierLib SHARED ${PreprocessBezierSrc})

add_executable(googleTest googleTest.cpp)
target_link_libraries(googleTest GTest::GTest GTest::Main)
add_test(google-test googleTest)

add_executable(manual-test manualTest.cpp)
target_link_libraries(manual-test PreprocessBezierLib)

add_executable(performance-test performanceTest.cpp)
target_link_libraries(performance-test PreprocessBezierLib)
