#!/bin/bash
# shellcheck disable=SC2086
hostname
strings gpurecsplit_construction | grep " -m"
strings simdrecsplit_construction | grep " -m"

for threads in $(seq 1 1 16); do
  params="--numQueries 0 --leafMethod rotations --leafSize 5 --bucketSize 5 --numThreads $threads"
  ./gpurecsplit_construction $params --numObjects $((2*threads))M
  ./simdrecsplit_construction $params --numObjects $((20*threads))M

  params="--numObjects $((20*threads))M --numQueries 0 --leafMethod rotations --leafSize 8 --bucketSize 100 --numThreads $threads"
  ./gpurecsplit_construction $params
  ./simdrecsplit_construction $params

  params="--numObjects $((20*threads))M --numQueries 0 --leafMethod rotations --leafSize 12 --bucketSize 9 --numThreads $threads"
  ./gpurecsplit_construction $params
  ./simdrecsplit_construction $params

  params="--numObjects $((200*threads))k --numQueries 0 --leafMethod rotations --leafSize 16 --bucketSize 2000 --numThreads $threads"
  ./gpurecsplit_construction $params
  ./simdrecsplit_construction $params
done
