cmake_minimum_required(VERSION 3.10)

project(guidescan-cli VERSION 2.0)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_CXX_FLAGS_DEBUG_INIT "-Wall -Wextra -Wfatal-errors -g -pg")
set(CMAKE_CXX_FLAGS_RELEASE_INIT "-Wall -Wextra -O3")

add_subdirectory(sdsl)
add_subdirectory(src bin)
