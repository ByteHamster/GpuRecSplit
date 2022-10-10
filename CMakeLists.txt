cmake_minimum_required (VERSION 3.16)
cmake_policy(SET CMP0048 NEW)
project(GPURecSplit LANGUAGES CXX)

# ---------------------------- Library Setup ----------------------------

# Standard RecSplit
add_library(RecSplit INTERFACE)
target_compile_features(RecSplit INTERFACE cxx_std_17)
target_include_directories(RecSplit INTERFACE src)
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
	target_compile_options(RecSplit INTERFACE -Wall -Wextra -march=native)
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
	target_compile_options(RecSplit INTERFACE /constexpr:steps10000000)
endif()

# GPURecSplit
include (CheckLanguage)
check_language(CUDA)
if (CMAKE_CUDA_COMPILER)
	enable_language(CUDA)
endif()
if (CMAKE_CUDA_COMPILER)
	add_library(GPURecSplit INTERFACE)
	target_include_directories(GPURecSplit INTERFACE src)
	target_compile_features(GPURecSplit INTERFACE cxx_std_17)
	target_compile_options(GPURecSplit INTERFACE --expt-relaxed-constexpr)
	set_target_properties(GPURecSplit PROPERTIES CUDA_SEPERABLE_COMPILATION ON)
	set_target_properties(GPURecSplit PROPERTIES CUDA_ARCHITECTURES 61-virtual)
	if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
		target_compile_options(GPURecSplit INTERFACE -Wall -Wextra -march=native)
		target_link_libraries(GPURecSplit INTERFACE pthread)
	elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
		target_compile_options(GPURecSplit INTERFACE -constexpr:steps10000000 -arch:AVX2 -DINSTRSET=8)
		target_include_directories(GPURecSplit INTERFACE src/support)
	endif()
endif()

# SIMDRecSplit
add_library(SIMDRecSplit INTERFACE)
target_include_directories(SIMDRecSplit INTERFACE src)
target_compile_features(SIMDRecSplit INTERFACE cxx_std_17)
target_compile_options(SIMDRecSplit INTERFACE -DSIMD)
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
	target_compile_options(SIMDRecSplit INTERFACE -Wall -Wextra -march=native)
	target_link_libraries(SIMDRecSplit INTERFACE pthread)
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
	target_compile_options(SIMDRecSplit INTERFACE /constexpr:steps10000000 /arch:AVX2 /DINSTRSET=8)
endif()

# ---------------------------- Benchmarks ----------------------------
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
	# Only added if this is the main cmake file

	# Standard RecSplit
	add_executable(recsplit_dump128 benchmark/function/recsplit_dump128.cpp)
	target_link_libraries(recsplit_dump128 PUBLIC RecSplit)

	add_executable(recsplit_interactive test/interactive.cpp)
	target_link_libraries(recsplit_interactive PUBLIC RecSplit)

	add_executable(recsplit_correctness test/correctness.cpp)
	target_link_libraries(recsplit_correctness PUBLIC RecSplit)

	add_executable(recsplit_construction benchmark/function/benchmark_construction.cpp)
	target_link_libraries(recsplit_construction PUBLIC RecSplit)

	# GPURecSplit
	if (CMAKE_CUDA_COMPILER)
		add_executable(gpurecsplit_dump128 benchmark/function/recsplit_dump128.cu)
		target_link_libraries(gpurecsplit_dump128 PUBLIC GPURecSplit)

		add_executable(gpurecsplit_interactive test/interactive.cu)
		target_link_libraries(gpurecsplit_interactive PUBLIC GPURecSplit)

		add_executable(gpurecsplit_correctness test/correctness.cu)
		target_link_libraries(gpurecsplit_correctness PUBLIC GPURecSplit)

		add_executable(gpurecsplit_construction benchmark/function/benchmark_construction.cu)
		target_link_libraries(gpurecsplit_construction PUBLIC GPURecSplit)

		add_executable(gpu_simd_equivalence test/GPUandSIMDequivalence.cu)
		target_link_libraries(gpurecsplit_construction PUBLIC GPURecSplit SIMDRecSplit)
	endif()

	# SIMDRecSplit
	add_executable(simdrecsplit_dump128 benchmark/function/recsplit_dump128.cpp)
	target_link_libraries(simdrecsplit_dump128 PUBLIC SIMDRecSplit)

	add_executable(simdrecsplit_interactive test/interactive.cpp)
	target_link_libraries(simdrecsplit_interactive PUBLIC SIMDRecSplit)

	add_executable(simdrecsplit_correctness test/correctness.cpp)
	target_link_libraries(simdrecsplit_correctness PUBLIC SIMDRecSplit)

	add_executable(simdrecsplit_construction benchmark/function/benchmark_construction.cpp)
	target_link_libraries(simdrecsplit_construction PUBLIC SIMDRecSplit)
endif()
