#!/bin/bash
hostname
strings gpurecsplit_construction | grep " -m"
strings simdrecsplit_construction | grep " -m"

# shellcheck disable=SC2086
for leafSize in $(seq 4 1 14); do
  params="--numObjects 20M --numQueries 0 --leafMethod rotations --leafSize $leafSize --bucketSize 2000"
  ./gpurecsplit_construction $params --numThreads 1
  ./simdrecsplit_construction $params --numThreads 1
  ./gpurecsplit_construction $params --numThreads 8
  ./simdrecsplit_construction $params --numThreads 8
done
