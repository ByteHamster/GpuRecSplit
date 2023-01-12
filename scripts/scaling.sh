#!/bin/bash
# shellcheck disable=SC2086
hostname
strings gpurecsplit_construction | grep " -m"
strings simdrecsplit_construction | grep " -m"

for threads in $(seq 1 1 16); do
  if [[ $threads -eq 0 ]]; then
    threads=1
  fi

  params="--numQueries 0 --leafMethod rotations --leafSize 5 --bucketSize 5 --numThreads $threads"
  #./gpurecsplit_construction $params --numObjects $((2*threads))M
  ./simdrecsplit_construction $params --numObjects $((80*threads))M
  #./simdrecsplit_construction $params --numObjects $((80*threads))M --partitioned

  params="--numQueries 0 --leafMethod rotations --leafSize 8 --bucketSize 100 --numThreads $threads"
  #./gpurecsplit_construction $params --numObjects $((20*threads))M
  ./simdrecsplit_construction $params --numObjects $((30*threads))M
  #./simdrecsplit_construction $params --numObjects $((30*threads))M --partitioned

  params="--numQueries 0 --leafMethod rotations --leafSize 12 --bucketSize 9 --numThreads $threads"
  #./gpurecsplit_construction $params --numObjects $((5*threads))M
  ./simdrecsplit_construction $params --numObjects $((10*threads))M
  #./simdrecsplit_construction $params --numObjects $((10*threads))M --partitioned

  params="--numQueries 0 --leafMethod rotations --leafSize 16 --bucketSize 2000 --numThreads $threads"
  #./gpurecsplit_construction $params --numObjects $((500*threads))k
  ./simdrecsplit_construction $params --numObjects $((20*threads))k
  #./simdrecsplit_construction $params --numObjects $((20*threads))k --partitioned
done
