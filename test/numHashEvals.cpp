#include <iostream>
#include <vector>
#include <cmath>
#include "../include/shockhash/TinyBinaryCuckooHashTable.h"

static constexpr uint64_t rotate(size_t l, uint64_t val, uint32_t x) {
return ((val << x) | (val >> (l - x))) & ((1 << l) - 1);
}

void testRotationFitting(size_t l) {
    size_t numIterations = std::max(2ul, (size_t) 7e7 / (1<<l));
    size_t totalTries = 0;
    size_t hfEvals = 0;
    for (size_t iteration = 0; iteration < numIterations; iteration++) {
        std::vector<shockhash::HashedKey> keys;
        for (size_t i = 0; i < l; i++) {
            keys.emplace_back(std::to_string(i) + " " + std::to_string(iteration));
        }
        while (true) {
            uint64_t taken1 = 0;
            uint64_t taken2 = 0;
            for (size_t i = 0; i < l; i++) {
                size_t pos = keys[i].hash(totalTries, l);
                hfEvals++;
                if (keys[i].mhc % 2 == 0) {
                    taken1 |= (1<<pos);
                } else {
                    taken2 |= (1<<pos);
                }
            }
            bool canBeRotated = false;
            for (size_t r = 0; r < l; r++) {
                if ((taken1 | rotate(l, taken2, r)) == (1<<l) - 1) {
                    canBeRotated = true;
                    break;
                }
            }
            if (canBeRotated) {
                break;
            }
            totalTries++;
        }
    }
    std::cout<<"RESULT"
             <<" method=rotations"
             <<" l="<<l
             <<" hfEvals="<<(double)hfEvals/(double)numIterations
             <<" tries="<<(double)totalTries / (double)numIterations
             <<" iterations="<<numIterations
             <<" spaceEstimate="<<log2((double)totalTries / (double)numIterations * l) / l
             <<std::endl;

}

void testBruteForce(size_t l) {
    size_t numIterations = std::max(2ul, (size_t) 4e7 / (1<<l));
    size_t totalTries = 0;
    size_t hfEvals = 0;
    for (size_t iteration = 0; iteration < numIterations; iteration++) {
        std::vector<shockhash::HashedKey> keys;
        for (size_t i = 0; i < l; i++) {
            keys.emplace_back(std::to_string(i) + " " + std::to_string(iteration));
        }
        while (true) {
            uint64_t taken = 0;
            size_t i = 0;
            for (; i < l; i++) {
                size_t pos = keys[i].hash(totalTries, l);
                hfEvals++;
                if (taken & (1<<pos)) {
                    break;
                }
                taken |= (1<<pos);
            }
            if (i == l) {
                break;
            }
            totalTries++;
        }
    }
    std::cout<<"RESULT"
             <<" method=bruteforce"
             <<" l="<<l
             <<" hfEvals="<<(double)hfEvals/(double)numIterations
             <<" tries="<<(double)totalTries / (double)numIterations
             <<" iterations="<<numIterations
             <<" spaceEstimate="<<log2((double)totalTries / (double)numIterations) / l
             <<std::endl;
}

void testShockHash(size_t l) {
    size_t numIterations = l <= 30 ? 20000 : (l <= 43 ? 4000 : 1000);
    size_t totalTries = 0;
    size_t hfEvals = 0;
    for (size_t iteration = 0; iteration < numIterations; iteration++) {
        shockhash::TinyBinaryCuckooHashTable table(l, l);
        for (size_t i = 0; i < l; i++) {
            table.prepare(shockhash::HashedKey(std::to_string(i) + " " + std::to_string(iteration)));
        }
        while (!table.construct(totalTries)) { // totalTries is the (unique) seed here
            totalTries++;
            hfEvals += 2 * l;
        }
    }
    std::cout<<"RESULT"
            <<" method=cuckoo"
            <<" l="<<l
            <<" hfEvals="<<(double)hfEvals/(double)numIterations
            <<" tries="<<(double)totalTries / (double)numIterations
            <<" iterations="<<numIterations
            <<" spaceEstimate="<<log2((double)totalTries / (double)numIterations) / l + 1
            <<std::endl;
}

int main() {
    for (size_t l = 2; l <= 20; l++) {
        if (l <= 25) {
            testRotationFitting(l);
        }
        if (l <= 22) {
            testBruteForce(l);
        }
        testShockHash(l);
    }
}
