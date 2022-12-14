cmake_minimum_required (VERSION 3.16)
cmake_policy(SET CMP0048 NEW)
project(GPURecSplit LANGUAGES CXX)

# ---------------------------- Dependencies ----------------------------
add_library(vectorclass INTERFACE)
target_include_directories(vectorclass SYSTEM INTERFACE extlib/vectorclass)

add_library(gcem INTERFACE)
target_include_directories(gcem SYSTEM INTERFACE extlib/gcem/include)

add_library(Sux INTERFACE)
target_include_directories(Sux SYSTEM INTERFACE extlib/sux)

add_subdirectory(extlib/simpleRibbon)

# ---------------------------- Library Setup ----------------------------

# Standard RecSplit
add_library(RecSplit INTERFACE)
target_compile_features(RecSplit INTERFACE cxx_std_20)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
	target_include_directories(RecSplit INTERFACE include)
else()
	target_include_directories(RecSplit SYSTEM INTERFACE include)
endif()
target_link_libraries(RecSplit INTERFACE gcem)
target_compile_options(RecSplit INTERFACE $<$<COMPILE_LANGUAGE:CXX>:-march=native>)

# GPURecSplit
include (CheckLanguage)
check_language(CUDA)
if (CMAKE_CUDA_COMPILER)
	enable_language(CUDA)

	add_library(GPURecSplit INTERFACE)
	if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
		target_include_directories(GPURecSplit INTERFACE include)
	else()
		target_include_directories(GPURecSplit SYSTEM INTERFACE include)
	endif()
	target_compile_features(GPURecSplit INTERFACE cxx_std_20)
	target_compile_options(GPURecSplit INTERFACE --expt-relaxed-constexpr)
	target_compile_options(GPURecSplit INTERFACE -Xcompiler="-march=native")
	target_link_libraries(GPURecSplit INTERFACE gcem pthread)
endif()

# SIMDRecSplit
add_library(SIMDRecSplit INTERFACE)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
	target_include_directories(SIMDRecSplit INTERFACE include)
else()
	target_include_directories(SIMDRecSplit SYSTEM INTERFACE include)
endif()
target_compile_features(SIMDRecSplit INTERFACE cxx_std_20)
target_compile_options(SIMDRecSplit INTERFACE -DSIMD)
target_link_libraries(SIMDRecSplit INTERFACE gcem vectorclass pthread)
target_compile_options(SIMDRecSplit INTERFACE $<$<COMPILE_LANGUAGE:CXX>:-march=native>)

# ShockHash
add_library(ShockHash INTERFACE)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
	target_include_directories(ShockHash INTERFACE include)
else()
	target_include_directories(ShockHash SYSTEM INTERFACE include)
endif()
target_compile_features(ShockHash INTERFACE cxx_std_20)
target_link_libraries(ShockHash INTERFACE Sux SimpleRibbon)

# ---------------------------- Benchmarks ----------------------------
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
	add_library(BenchmarkUtils INTERFACE)
	target_include_directories(BenchmarkUtils INTERFACE benchmark)
	target_include_directories(BenchmarkUtils INTERFACE test)
	target_include_directories(BenchmarkUtils INTERFACE extlib/util/include)

	set(TLX_INSTALL_INCLUDE_DIR tlx CACHE PATH "Workaround for TLX breaking the first cmake call")
	add_subdirectory(extlib/tlx EXCLUDE_FROM_ALL)
	target_link_libraries(BenchmarkUtils INTERFACE tlx)

	# Warnings if this is the main project
	target_compile_options(RecSplit INTERFACE $<$<COMPILE_LANGUAGE:CXX>:-Wall -Wextra -frecord-gcc-switches>)
	target_compile_options(SIMDRecSplit INTERFACE $<$<COMPILE_LANGUAGE:CXX>:-Wall -Wextra -frecord-gcc-switches>)
	if (CMAKE_CUDA_COMPILER)
		target_compile_options(GPURecSplit INTERFACE -Xcompiler="-Wall"  -Xcompiler="-Wextra" -Xcompiler="-frecord-gcc-switches")
	endif()

	# Standard RecSplit
	add_executable(recsplit_dump128 benchmark/recsplit_dump128.cpp)
	target_link_libraries(recsplit_dump128 PUBLIC RecSplit BenchmarkUtils)

	add_executable(recsplit_interactive test/interactive.cpp)
	target_link_libraries(recsplit_interactive PUBLIC RecSplit BenchmarkUtils)

	add_executable(recsplit_correctness test/correctness.cpp)
	target_link_libraries(recsplit_correctness PUBLIC RecSplit BenchmarkUtils)

	add_executable(recsplit_construction benchmark/benchmark_construction.cpp)
	target_link_libraries(recsplit_construction PUBLIC RecSplit BenchmarkUtils ShockHash)

	# GPURecSplit
	if (CMAKE_CUDA_COMPILER)
		function(prepareCudaTarget target)
			# Need to be set on the top-level target, not the GPURecSplit target. Otherwise cmake complains
			set_target_properties(${target} PROPERTIES CUDA_SEPERABLE_COMPILATION ON)
			set_target_properties(${target} PROPERTIES CUDA_ARCHITECTURES "50;61;72")
		endfunction()

		add_executable(gpurecsplit_dump128 benchmark/recsplit_dump128.cu)
		target_link_libraries(gpurecsplit_dump128 PUBLIC GPURecSplit BenchmarkUtils)
		prepareCudaTarget(gpurecsplit_dump128)

		add_executable(gpurecsplit_interactive test/interactive.cu)
		target_link_libraries(gpurecsplit_interactive PUBLIC GPURecSplit BenchmarkUtils)
		prepareCudaTarget(gpurecsplit_interactive)

		add_executable(gpurecsplit_correctness test/correctness.cu)
		target_link_libraries(gpurecsplit_correctness PUBLIC GPURecSplit BenchmarkUtils)
		prepareCudaTarget(gpurecsplit_correctness)

		add_executable(gpurecsplit_construction benchmark/benchmark_construction.cu)
		target_link_libraries(gpurecsplit_construction PUBLIC GPURecSplit BenchmarkUtils)
		prepareCudaTarget(gpurecsplit_construction)

		add_executable(gpu_simd_equivalence test/GPUandSIMDequivalence.cu)
		target_link_libraries(gpu_simd_equivalence PUBLIC GPURecSplit SIMDRecSplit BenchmarkUtils)
		prepareCudaTarget(gpu_simd_equivalence)
	endif()

	# SIMDRecSplit
	add_executable(simdrecsplit_dump128 benchmark/recsplit_dump128.cpp)
	target_link_libraries(simdrecsplit_dump128 PUBLIC SIMDRecSplit BenchmarkUtils)

	add_executable(simdrecsplit_interactive test/interactive.cpp)
	target_link_libraries(simdrecsplit_interactive PUBLIC SIMDRecSplit BenchmarkUtils)

	add_executable(simdrecsplit_correctness test/correctness.cpp)
	target_link_libraries(simdrecsplit_correctness PUBLIC SIMDRecSplit BenchmarkUtils)

	add_executable(simdrecsplit_construction benchmark/benchmark_construction.cpp)
	target_link_libraries(simdrecsplit_construction PUBLIC SIMDRecSplit BenchmarkUtils)
endif()
