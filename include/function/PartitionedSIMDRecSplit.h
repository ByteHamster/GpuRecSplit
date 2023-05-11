#pragma once
#include <function/SIMDRecSplit.hpp>
#include <thread>
#include <Function.h>

namespace bez::function {
/**
 * Multiple SIMDRecSplit perfect hash functions, one for each construction thread.
 * If you use one single thread, this is slower than simply using SIMDRecSplit directly.
 * This also has some query overhead.
 */
template<size_t leafSize>
class PartitionedSIMDRecSplit {
    private:
        size_t numThreads;
        static constexpr size_t HASH_FUNCTION_CHILD_ASSIGNMENT = 43;
        std::vector<SIMDRecSplit<leafSize>*> children;
        std::vector<uint64_t> childOffsets;

    public:
        PartitionedSIMDRecSplit(const vector<string> &keys, const size_t bucket_size, size_t numThreads)
                : numThreads(numThreads) {
            std::vector<std::vector<hash128_t>> childInput;
            childInput.resize(numThreads);
            const size_t N = keys.size();
            for (auto &singleChildInput: childInput) {
                singleChildInput.reserve(N / numThreads);
            }
            for (const std::string &key : keys) {
                hash128_t hash = first_hash(key.c_str(), key.size());
                // RecSplit internally uses hash.first for bucket assignment, but not for anything else.
                // So we should be fine using remixed hash.first for partition assignment.
                size_t child = ::util::fastrange64(remix(hash.first), numThreads);
                childInput[child].push_back(hash);
            }
            buildChildren(childInput, bucket_size);
        }

        PartitionedSIMDRecSplit(const vector<hash128_t> &keys, const size_t bucket_size, size_t numThreads)
                : numThreads(numThreads) {
            std::vector<std::vector<hash128_t>> childInput;
            childInput.resize(numThreads);
            const size_t N = keys.size();
            for (auto &singleChildInput: childInput) {
                singleChildInput.reserve(N / numThreads);
            }
            for (const hash128_t &hash : keys) {
                size_t child = ::util::fastrange64(remix(hash.first), numThreads);
                childInput[child].push_back(hash);
            }
            buildChildren(childInput, bucket_size);
        }

    private:
        void buildChildren(std::vector<std::vector<hash128_t>> &childInput, size_t bucket_size) {
            children.resize(numThreads);
            std::vector<std::thread> threads;
            std::atomic<bool> hadException = false;
            uint64_t childOffset = 0;
            for (size_t i = 0; i < numThreads; i++) {
                childOffsets.push_back(childOffset);
                childOffset += childInput[i].size();
                threads.emplace_back([&, i]() {
                    try {
                        children[i] = new SIMDRecSplit<leafSize>(childInput[i], bucket_size, 1);
                    } catch (const std::exception& e) {
                        std::cout<<"Error: "<<e.what()<<std::endl;
                        hadException = true;
                    }
                });
            }
            for (size_t i = 0; i < numThreads; i++) {
                threads[i].join();
            }
            if (hadException) {
                throw std::logic_error("One construction thread experienced a problem. Read output for details.");
            }
        }

    public:
        ~PartitionedSIMDRecSplit() {
            for (auto &child : children) {
                delete child;
            }
        }

        /** Estimate for the space usage of this structure, in bits */
        [[nodiscard]] size_t getBits() const {
            size_t spaceUsage = sizeof(*this) * 8;
            for (auto &child : children) {
                spaceUsage += child->getBits();
            }
            return spaceUsage;
        }

        size_t operator()(const hash128_t &hash) const {
            size_t child = ::util::fastrange64(remix(hash.first), numThreads);
            return children[child]->operator()(hash) + childOffsets[child];
        }

        size_t operator() (std::string &key) const {
            return operator()(first_hash(key.c_str(), key.size()));
        }
};
} // Namespace bez::function
