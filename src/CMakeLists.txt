add_executable(guidescan guidescan.cxx
  genomics/seq_io.cxx
  genomics/kmer.cxx
  genomics/structures.cxx
  genomics/sequences.cxx )

target_include_directories(guidescan PUBLIC 
  "${CMAKE_SOURCE_DIR}/include"
  "${PROJECT_BINARY_DIR}/sdsl/include"
  "${PROJECT_BINARY_DIR}/sdsl/external/libdivsufsort/include")

# What if pthread isn't found? Find alternatives...
target_link_libraries(guidescan PUBLIC sdsl divsufsort divsufsort64 pthread)

# Static linking trick comes from here:
# https://stackoverflow.com/questions/35116327/when-g-static-link-pthread-cause-segmentation-fault-why
# pretty dirty but that is how it goes.
if(LINK_STATICALLY)
  set(CMAKE_EXE_LINKER_FLAGS "-static")
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lrt -Wl,--whole-archive -lpthread -Wl,--no-whole-archive")
endif()
