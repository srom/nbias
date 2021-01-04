cmake_minimum_required(VERSION 3.0.0)

project(nbias)

# Compile with C++17
set(CMAKE_CXX_STANDARD 17)

# External libraries (typically single header libraries)
include_directories(external)

# Boost library
find_package(Boost 1.75 COMPONENTS program_options iostreams filesystem unit_test_framework REQUIRED)

# Executables
file(GLOB tri_bias_SOURCES "src/tri_bias/tri_bias.cpp")
add_executable(tri_bias ${tri_bias_SOURCES})
target_link_libraries(tri_bias ${Boost_LIBRARIES})

# Test Executables
file(GLOB test_fasta_SOURCES "test/fasta/*.cpp")
add_executable(test_fasta ${test_fasta_SOURCES})
target_link_libraries(test_fasta ${Boost_LIBRARIES})

file(GLOB test_assembly_SOURCES "test/assembly/*.cpp")
add_executable(test_assembly ${test_assembly_SOURCES})
target_link_libraries(test_assembly ${Boost_LIBRARIES})

file(GLOB test_tri_count_SOURCES "test/count/*.cpp")
add_executable(test_tri_count ${test_tri_count_SOURCES})
target_link_libraries(test_tri_count ${Boost_LIBRARIES})

file(GLOB test_dna_SOURCES "test/dna/*.cpp")
add_executable(test_dna ${test_dna_SOURCES})
target_link_libraries(test_dna ${Boost_LIBRARIES})
