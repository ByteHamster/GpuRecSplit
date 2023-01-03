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

repeat  2 ./recsplit_construction     --numObjects 1M --numQueries 50M --leafMethod bruteforce --leafSize 16 --bucketSize 2000
repeat  3 ./simdrecsplit_construction --numObjects 1M --numQueries 50M --leafMethod rotations  --leafSize 16 --bucketSize 2000 --numThreads 1
repeat  5 ./simdrecsplit_construction --numObjects 1M --numQueries 50M --leafMethod rotations  --leafSize 16 --bucketSize 2000 --numThreads 16
repeat 10 ./gpurecsplit_construction  --numObjects 1M --numQueries 50M --leafMethod rotations  --leafSize 16 --bucketSize 2000 --numThreads 1

repeat  1 ./recsplit_construction     --numObjects 1M --numQueries 50M --leafMethod bruteforce --leafSize 18 --bucketSize 50
repeat  3 ./simdrecsplit_construction --numObjects 1M --numQueries 50M --leafMethod rotations  --leafSize 18 --bucketSize 50 --numThreads 1
repeat  5 ./simdrecsplit_construction --numObjects 1M --numQueries 50M --leafMethod rotations  --leafSize 18 --bucketSize 50 --numThreads 16
repeat 10 ./gpurecsplit_construction  --numObjects 1M --numQueries 50M --leafMethod rotations  --leafSize 18 --bucketSize 50 --numThreads 1

repeat  1 ./gpurecsplit_construction  --numObjects 1M --numQueries 50M --leafMethod rotations  --leafSize 24 --bucketSize 2000 --numThreads 1
