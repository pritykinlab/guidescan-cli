cmake_minimum_required(VERSION 2.8.11)

project(guidescan-cli)
set(PROJECT_VERSION VERSION 2.0)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_CXX_FLAGS_DEBUG_INIT "-Wall -Wextra -Wfatal-errors -g -pg")
set(CMAKE_CXX_FLAGS_RELEASE_INIT "-Wall -Wextra -O3 ")

# Static linking trick comes from here:
# https://stackoverflow.com/questions/35116327/when-g-static-link-pthread-cause-segmentation-fault-why
# pretty dirty but that is how it goes.
IF(CMAKE_BUILD_TYPE MATCHES RELEASE)
  set(CMAKE_EXE_LINKER_FLAGS "-static")
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lrt -Wl,--whole-archive -lpthread -Wl,--no-whole-archive")
ENDIF(CMAKE_BUILD_TYPE MATCHES RELEASE)

add_subdirectory(sdsl)
add_subdirectory(src bin)
