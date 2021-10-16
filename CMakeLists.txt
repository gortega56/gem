cmake_minimum_required(VERSION 3.20)

project(gem VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
set(CMAKE_CXX_EXTENSIONS FALSE)

add_library(gem INTERFACE)

set(SOURCE_FILES ${CMAKE_CURRENT_SOURCE_DIR}/include/vector.h
                 ${CMAKE_CURRENT_SOURCE_DIR}/include/matrix.h
                 ${CMAKE_CURRENT_SOURCE_DIR}/include/quaternion.h
                 ${CMAKE_CURRENT_SOURCE_DIR}/include/transform.h
                 ${CMAKE_CURRENT_SOURCE_DIR}/include/ray.h
                 ${CMAKE_CURRENT_SOURCE_DIR}/include/line.h
                 ${CMAKE_CURRENT_SOURCE_DIR}/include/plane.h
                 ${CMAKE_CURRENT_SOURCE_DIR}/include/range.h
                 ${CMAKE_CURRENT_SOURCE_DIR}/include/common/defines.h
)

target_sources(${PROJECT_NAME} INTERFACE "$<BUILD_INTERFACE:${SOURCE_FILES}>")

target_include_directories(${PROJECT_NAME} SYSTEM INTERFACE 
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include">
    "$<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_INCLUDEDIR}>")

add_custom_target(game_engine_math SOURCES ${SOURCE_FILES})

enable_testing()

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/unit_tests)