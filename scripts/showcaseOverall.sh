#!/bin/bash

function repeat() {
    repetitions=$1
    shift
    # shellcheck disable=SC2034
    for i in $(seq "$repetitions"); do
      # shellcheck disable=SC2068
      $@
    done
}

repeat 2 ./recsplit_construction     --numObjects 5M --numQueries 50M --leafMethod bruteforce --leafSize 16 --bucketSize 2000 --numThreads 1
repeat 2 ./recsplit_construction     --numObjects 5M --numQueries 50M --leafMethod bruteforce --leafSize 16 --bucketSize 2000 --numThreads 16
repeat 2 ./simdrecsplit_construction --numObjects 5M --numQueries 50M --leafMethod rotations  --leafSize 16 --bucketSize 2000 --numThreads 1
repeat 3 ./simdrecsplit_construction --numObjects 5M --numQueries 50M --leafMethod rotations  --leafSize 16 --bucketSize 2000 --numThreads 16
repeat 5 ./gpurecsplit_construction  --numObjects 5M --numQueries 50M --leafMethod rotations  --leafSize 16 --bucketSize 2000 --numThreads 1

repeat 1 ./recsplit_construction     --numObjects 5M --numQueries 50M --leafMethod bruteforce --leafSize 18 --bucketSize 50 --numThreads 1
repeat 1 ./recsplit_construction     --numObjects 5M --numQueries 50M --leafMethod bruteforce --leafSize 18 --bucketSize 50 --numThreads 16
repeat 2 ./simdrecsplit_construction --numObjects 5M --numQueries 50M --leafMethod rotations  --leafSize 18 --bucketSize 50 --numThreads 1
repeat 3 ./simdrecsplit_construction --numObjects 5M --numQueries 50M --leafMethod rotations  --leafSize 18 --bucketSize 50 --numThreads 16
repeat 5 ./gpurecsplit_construction  --numObjects 5M --numQueries 50M --leafMethod rotations  --leafSize 18 --bucketSize 50 --numThreads 1

repeat 1 ./gpurecsplit_construction  --numObjects 5M --numQueries 50M --leafMethod rotations  --leafSize 24 --bucketSize 2000 --numThreads 1
