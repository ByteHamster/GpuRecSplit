#include <util/Sorter.hpp>
#include <ips2ra.hpp>
#include <ips4o.hpp>

void bez::sorter::sortParallel_hash128_t(bez::function::hash128_t* objects, size_t n, size_t numThreads) {
    //ips2ra::parallel::sort(objects, objects + n, [] (bez::function::hash128_t t) { return t.first; }, numThreads);
    ips4o::parallel::sort(objects, objects + n, [] (bez::function::hash128_t x, bez::function::hash128_t y) { return x.first < y.first; }, numThreads);
}
