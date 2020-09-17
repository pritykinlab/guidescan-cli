add_executable(main seq_io.cxx main.cxx)

target_include_directories(main PUBLIC 
  "${CMAKE_SOURCE_DIR}/include"
  "${PROJECT_BINARY_DIR}/sdsl/include"
  "${PROJECT_BINARY_DIR}/sdsl/external/libdivsufsort/include")

target_link_libraries(main PUBLIC sdsl divsufsort divsufsort64)
