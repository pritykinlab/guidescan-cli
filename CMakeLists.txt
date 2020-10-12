cmake_minimum_required(VERSION 2.8.11)

option(LINK_STATICALLY "Statically link the executable against c++11 standard libary" OFF)

project(guidescan-cli)

set(PROJECT_VERSION VERSION 2.0)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -Wfatal-errors -g -pg")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall -Wextra")

# set(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -pg")

add_subdirectory(sdsl)
add_subdirectory(src bin)
add_subdirectory(test test_bin)
