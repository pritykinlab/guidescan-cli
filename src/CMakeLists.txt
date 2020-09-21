add_executable(guidescan seq_io.cxx kmer.cxx guidescan.cxx)

target_include_directories(guidescan PUBLIC 
  "${CMAKE_SOURCE_DIR}/include"
  "${PROJECT_BINARY_DIR}/sdsl/include"
  "${PROJECT_BINARY_DIR}/sdsl/external/libdivsufsort/include")

target_link_libraries(guidescan PUBLIC sdsl divsufsort divsufsort64 pthread)
