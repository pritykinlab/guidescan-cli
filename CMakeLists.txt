cmake_minimum_required(VERSION 3.10)

project(guidescan-cli VERSION 2.0)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_subdirectory(sdsl)
add_subdirectory(src bin)
