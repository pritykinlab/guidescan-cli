cmake_minimum_required(VERSION 2.8.11)

option(LINK_STATICALLY "Statically link the executable against c++11 standard libary" OFF)

project(guidescan-cli)

set(PROJECT_VERSION VERSION 2.0)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(GUIDESCAN_VERSION "2.0.0" CACHE STRING "The version number of GuideScan")
if(NOT GUIDESCAN_VERSION STREQUAL "2.0.0")
    configure_file (
        "${PROJECT_SOURCE_DIR}/configure/version.hpp.in"
        "${PROJECT_SOURCE_DIR}/include/version.hpp"
    )
endif()

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -Wfatal-errors -g -pg")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall -Wextra")

# set(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -pg")

add_subdirectory(sdsl)
add_subdirectory(src bin)
add_subdirectory(test test_bin)
