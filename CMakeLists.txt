cmake_minimum_required (VERSION 3.25...4.0)
cmake_policy(SET CMP0048 NEW)
project(GPURecSplit LANGUAGES CXX)

if(TARGET RecSplit)
	return()
endif()

if (NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE "Release")
endif ()

# ---------------------------- Dependencies ----------------------------
if(NOT TARGET vectorclass)
	add_library(vectorclass INTERFACE)
	target_include_directories(vectorclass SYSTEM INTERFACE extlib/vectorclass)
endif()

if(NOT TARGET gcem)
	add_library(gcem INTERFACE)
	target_include_directories(gcem SYSTEM INTERFACE extlib/gcem/include)
endif()

if(NOT TARGET Sux)
	add_library(Sux INTERFACE)
	target_include_directories(Sux SYSTEM INTERFACE extlib/sux)
endif()

set(TLX_INSTALL_INCLUDE_DIR tlx CACHE PATH "Workaround for TLX breaking the first cmake call")
add_subdirectory(extlib/tlx EXCLUDE_FROM_ALL)

if(NOT TARGET ips2ra)
	set(IPS2RA_DISABLE_PARALLEL ON CACHE PATH "ips2ra's FindTBB greps a file that does not exist in recent TBB versions")
	add_subdirectory(extlib/ips2ra EXCLUDE_FROM_ALL)
	find_package(TBB)
	target_compile_options(ips2ra INTERFACE -D_REENTRANT)
	target_link_libraries(ips2ra INTERFACE pthread atomic TBB::tbb)
endif()

if(NOT TARGET Ips2raRecSplitSorter)
	add_library(Ips2raRecSplitSorter SHARED src/Sorter.cpp)
	target_compile_features(Ips2raRecSplitSorter PRIVATE cxx_std_20)
	target_include_directories(Ips2raRecSplitSorter PRIVATE include)
	target_link_libraries(Ips2raRecSplitSorter PUBLIC ips2ra tlx Sux)
	target_compile_options(Ips2raRecSplitSorter PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-march=native>)
endif()

# ---------------------------- Library Setup ----------------------------

# Standard RecSplit
add_library(RecSplit INTERFACE)
target_compile_features(RecSplit INTERFACE cxx_std_20)
target_include_directories(RecSplit INTERFACE include)
target_link_libraries(RecSplit INTERFACE gcem pthread Ips2raRecSplitSorter Sux)
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
	target_compile_options(GPURecSplit INTERFACE $<$<COMPILE_LANGUAGE:CUDA>: --expt-relaxed-constexpr -Xcompiler="-march=native" -diag-suppress 177 >)
	target_link_libraries(GPURecSplit INTERFACE gcem pthread Ips2raRecSplitSorter Sux)
else()
	message(STATUS "No CUDA compiler found, skipping GPURecSplit target")
endif()

# SIMDRecSplit
include(${CMAKE_CURRENT_SOURCE_DIR}/extlib/cmake-findsse/FindSSE.cmake)
FindSSE()
if(SSE4_2_FOUND)
	add_library(SIMDRecSplit INTERFACE)
	if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
		target_include_directories(SIMDRecSplit INTERFACE include)
	else()
		target_include_directories(SIMDRecSplit SYSTEM INTERFACE include)
	endif()
	target_compile_features(SIMDRecSplit INTERFACE cxx_std_20)
	target_compile_options(SIMDRecSplit INTERFACE -DSIMD)
	target_link_libraries(SIMDRecSplit INTERFACE gcem vectorclass pthread Ips2raRecSplitSorter Sux)
	target_compile_options(SIMDRecSplit INTERFACE $<$<COMPILE_LANGUAGE:CXX>:-march=native>)
else()
	message(STATUS "No SIMD support found, skipping SIMDRecSplit target")
endif(SSE4_2_FOUND)

# ---------------------------- Benchmarks ----------------------------
if(PROJECT_IS_TOP_LEVEL)
	target_compile_options(RecSplit INTERFACE $<$<COMPILE_LANGUAGE:CXX>:-Wall -Wextra -Wpedantic -Werror -frecord-gcc-switches>)

	add_library(BenchmarkUtils INTERFACE)
	target_include_directories(BenchmarkUtils INTERFACE benchmark)
	target_include_directories(BenchmarkUtils INTERFACE test)
	target_link_libraries(BenchmarkUtils INTERFACE tlx)

	add_subdirectory(extlib/util EXCLUDE_FROM_ALL)
	target_link_libraries(BenchmarkUtils INTERFACE ByteHamster::Util)

	# Warnings if this is the main project
	target_compile_options(RecSplit INTERFACE $<$<COMPILE_LANGUAGE:CXX>:-Wall -Wextra -frecord-gcc-switches>)
	if(SSE4_2_FOUND)
		target_compile_options(SIMDRecSplit INTERFACE $<$<COMPILE_LANGUAGE:CXX>:-Wall -Wextra -Wpedantic -Werror -frecord-gcc-switches>)
	endif()
	if (CMAKE_CUDA_COMPILER)
		target_compile_options(GPURecSplit INTERFACE -Xcompiler="-Wall" -Xcompiler="-Wextra" -Xcompiler="-Werror" -Xcompiler="-frecord-gcc-switches")
	endif()

	# Standard RecSplit
	add_executable(recsplit_dump128 benchmark/recsplit_dump128.cpp)
	target_link_libraries(recsplit_dump128 PUBLIC RecSplit BenchmarkUtils)

	add_executable(recsplit_interactive test/interactive.cpp)
	target_link_libraries(recsplit_interactive PUBLIC RecSplit BenchmarkUtils)

	add_executable(recsplit_correctness test/correctness.cpp)
	target_link_libraries(recsplit_correctness PUBLIC RecSplit BenchmarkUtils)

	add_executable(recsplit_construction benchmark/benchmark_construction.cpp)
	target_link_libraries(recsplit_construction PUBLIC RecSplit BenchmarkUtils)

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
	if(SSE4_2_FOUND)
		add_executable(simdrecsplit_dump128 benchmark/recsplit_dump128.cpp)
		target_link_libraries(simdrecsplit_dump128 PUBLIC SIMDRecSplit BenchmarkUtils)

		add_executable(simdrecsplit_interactive test/interactive.cpp)
		target_link_libraries(simdrecsplit_interactive PUBLIC SIMDRecSplit BenchmarkUtils)

		add_executable(simdrecsplit_correctness test/correctness.cpp)
		target_link_libraries(simdrecsplit_correctness PUBLIC SIMDRecSplit BenchmarkUtils)

		add_executable(simdrecsplit_construction benchmark/benchmark_construction.cpp)
		target_link_libraries(simdrecsplit_construction PUBLIC SIMDRecSplit BenchmarkUtils)
	endif()
endif()
