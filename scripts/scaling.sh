#!/bin/bash
# shellcheck disable=SC2086
hostname
strings gpurecsplit_construction | grep " -m"
strings simdrecsplit_construction | grep " -m"

for threads in $(seq 1 1 8); do
  params="--numObjects $((20*threads))M --numQueries 0 --leafMethod rotations --leafSize 5 --bucketSize 5 --numThreads $threads"
  [ $threads == 1 ] && ./recsplit_construction $params
  ./gpurecsplit_construction $params
  ./simdrecsplit_construction $params

  params="--numObjects $((20*threads))M --numQueries 0 --leafMethod rotations --leafSize 8 --bucketSize 100 --numThreads $threads"
  [ $threads == 1 ] && ./recsplit_construction $params
  ./gpurecsplit_construction $params
  ./simdrecsplit_construction $params

  params="--numObjects $((20*threads))M --numQueries 0 --leafMethod rotations --leafSize 16 --bucketSize 2000 --numThreads $threads"
  [ $threads == 1 ] && ./recsplit_construction $params
  ./gpurecsplit_construction $params
  ./simdrecsplit_construction $params
done
