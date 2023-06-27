cmake_minimum_required(VERSION 3.26)

project(gem VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
set(CMAKE_CXX_EXTENSIONS FALSE)

add_library(gem INTERFACE)

target_sources(gem
    PRIVATE    
    ${CMAKE_CURRENT_SOURCE_DIR}/include/vector.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/matrix.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/quaternion.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/transform.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/ray.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/line.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/plane.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/range.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/common/defines.h
)

target_include_directories(gem
  INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/include)

set(GEM_UNIT_TEST_ENABLE OFF CACHE BOOL "Builds unit tests using gtest")
if (GEM_UNIT_TEST_ENABLE)
    enable_testing()

    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/unit_tests)
endif()