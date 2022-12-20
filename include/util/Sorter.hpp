#include "Hash128.h"

namespace bez::sorter {
    /**
     * Reasons for moving ips2ra to its own compilation unit:
     * - nvcc somehow doesn't like ips2ra
     * - Compilation is faster
     */
    void sortParallel_hash128_t(bez::function::hash128_t *objects, size_t n, size_t numThreads);
}
