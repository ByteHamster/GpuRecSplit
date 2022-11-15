#define GPU
#include "benchmark_construction.hpp"

int main(int argc, const char* const* argv) {
    setenv ("CUDA_DEVICE_MAX_CONNECTIONS", "32", 1);
    return constructAll(argc, argv);
}
