add_executable(guidescan guidescan.cxx
  genomics/seq_io.cxx
  genomics/kmer.cxx
  genomics/structures.cxx
  genomics/sequences.cxx )

target_include_directories(guidescan PUBLIC 
  "${CMAKE_SOURCE_DIR}/include"
  "${CMAKE_SOURCE_DIR}/include/genomics"
  "${PROJECT_BINARY_DIR}/sdsl/include"
  "${PROJECT_BINARY_DIR}/sdsl/external/libdivsufsort/include")

target_link_libraries(guidescan PUBLIC sdsl divsufsort divsufsort64 pthread)
