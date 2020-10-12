add_executable(csa_sada_test csa_sada_test.cxx)

target_include_directories(csa_sada_test PUBLIC 
  "${CMAKE_SOURCE_DIR}/include"
  "${PROJECT_BINARY_DIR}/sdsl/include"
  "${PROJECT_BINARY_DIR}/sdsl/external/libdivsufsort/include")

# What if pthread isn't found? Find alternatives...
target_link_libraries(csa_sada_test PUBLIC sdsl divsufsort divsufsort64)
