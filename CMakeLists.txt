cmake_minimum_required(VERSION 3.10)

project(Game-Development-Mathematics VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(GDM_BUILD_TESTS "Build tests" OFF)

file(GLOB_RECURSE GDM_SOURCES "src/*.cpp" "src/*.hpp")

list(REMOVE_ITEM GDM_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/precompiled_headers.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/src/precompiled_headers.hpp")

add_subdirectory(src)

if(GDM_BUILD_TESTS)
    add_subdirectory(test)
endif()
