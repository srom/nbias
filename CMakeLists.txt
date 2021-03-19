cmake_minimum_required(VERSION 3.0.0)

project(nbias CXX)

enable_testing()

# Compile with C++17
set(CMAKE_CXX_STANDARD 17)

# External libraries (typically single header libraries)
include_directories(external)

# Boost library
find_package(Boost 1.75 COMPONENTS program_options iostreams filesystem unit_test_framework REQUIRED)

# XTensor library
set(XTENSOR_USE_XSIMD 1)
find_package(xtensor REQUIRED)

# Executables
file(GLOB compute_distribution_SOURCES "src/compute_distribution.cpp")
add_executable(compute_distribution ${compute_distribution_SOURCES})
target_link_libraries(compute_distribution ${Boost_LIBRARIES})

file(GLOB compute_distance_SOURCES "src/compute_distance.cpp")
add_executable(compute_distance ${compute_distance_SOURCES})
target_link_libraries(compute_distance ${Boost_LIBRARIES} xtensor)
target_compile_options(compute_distance PRIVATE -march=native)

file(GLOB domain_probability_SOURCES "src/domain_probability.cpp")
add_executable(domain_probability ${domain_probability_SOURCES})
target_link_libraries(domain_probability ${Boost_LIBRARIES} xtensor)
target_compile_options(domain_probability PRIVATE -march=native)

file(GLOB intervention_SOURCES "src/intervention.cpp")
add_executable(intervention ${intervention_SOURCES})
target_link_libraries(intervention ${Boost_LIBRARIES} xtensor)
target_compile_options(intervention PRIVATE -march=native)

# Test Executables
file(GLOB test_fasta_SOURCES "test/fasta/*.cpp")
add_executable(test_fasta ${test_fasta_SOURCES})
target_link_libraries(test_fasta ${Boost_LIBRARIES})
add_test(NAME fasta COMMAND test_fasta)

file(GLOB test_assembly_SOURCES "test/assembly/*.cpp")
add_executable(test_assembly ${test_assembly_SOURCES})
target_link_libraries(test_assembly ${Boost_LIBRARIES})
add_test(NAME assembly COMMAND test_assembly)

file(GLOB test_tri_count_SOURCES "test/count/*.cpp")
add_executable(test_tri_count ${test_tri_count_SOURCES})
target_link_libraries(test_tri_count ${Boost_LIBRARIES})
add_test(NAME tri_count COMMAND test_tri_count)

file(GLOB test_dna_SOURCES "test/dna/*.cpp")
add_executable(test_dna ${test_dna_SOURCES})
target_link_libraries(test_dna ${Boost_LIBRARIES})
add_test(NAME dna COMMAND test_dna)

file(GLOB test_distance_SOURCES "test/distance/*.cpp")
add_executable(test_distance ${test_distance_SOURCES})
target_link_libraries(test_distance ${Boost_LIBRARIES} xtensor)
target_compile_options(test_distance PRIVATE -march=native)
add_test(NAME distance COMMAND test_distance)

file(GLOB test_distribution_SOURCES "test/distribution/*.cpp")
add_executable(test_distribution ${test_distribution_SOURCES})
target_link_libraries(test_distribution ${Boost_LIBRARIES} xtensor)
target_compile_options(test_distribution PRIVATE -march=native)
add_test(NAME distribution COMMAND test_distribution)

file(GLOB test_probability_SOURCES "test/probability/*.cpp")
add_executable(test_probability ${test_probability_SOURCES})
target_link_libraries(test_probability ${Boost_LIBRARIES} xtensor)
target_compile_options(test_probability PRIVATE -march=native)
add_test(NAME probability COMMAND test_probability)

file(GLOB test_domain_SOURCES "test/domain/*.cpp")
add_executable(test_domain ${test_domain_SOURCES})
target_link_libraries(test_domain ${Boost_LIBRARIES} xtensor)
target_compile_options(test_domain PRIVATE -march=native)
add_test(NAME domain COMMAND test_domain)

file(GLOB test_codons_SOURCES "test/codons/*.cpp")
add_executable(test_codons ${test_codons_SOURCES})
target_link_libraries(test_codons ${Boost_LIBRARIES})
add_test(NAME codons COMMAND test_codons)

file(GLOB test_context_SOURCES "test/task/*.cpp")
add_executable(test_context ${test_context_SOURCES})
target_link_libraries(test_context ${Boost_LIBRARIES})
add_test(NAME context COMMAND test_context)
