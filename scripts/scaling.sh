#!/bin/bash
# shellcheck disable=SC2086
hostname
strings gpurecsplit_construction | grep " -m"
strings simdrecsplit_construction | grep " -m"

for threads in $(seq 1 1 16); do
  if [[ $threads -eq 0 ]]; then
    threads=1
  fi
  params="--numQueries 0 --leafMethod rotations --numThreads $threads"
  ./simdrecsplit_construction $params --leafSize  5 --bucketSize    5 --numObjects 1200M
  ./simdrecsplit_construction $params --leafSize  8 --bucketSize  100 --numObjects 400M
  ./simdrecsplit_construction $params --leafSize 12 --bucketSize    9 --numObjects 160M
  ./simdrecsplit_construction $params --leafSize 16 --bucketSize 2000 --numObjects 320k
done
