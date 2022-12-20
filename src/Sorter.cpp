#include <util/Sorter.hpp>
#include <ips2ra.hpp>

void bez::sorter::sortParallel_hash128_t(bez::function::hash128_t* objects, size_t n, size_t numThreads) {
    ips2ra::parallel::sort(objects, objects + n, [] (bez::function::hash128_t t) { return t.first; }, numThreads);
}
