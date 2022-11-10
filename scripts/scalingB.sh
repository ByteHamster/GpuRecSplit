#!/bin/bash
hostname
strings gpurecsplit_construction | grep " -m"
strings simdrecsplit_construction | grep " -m"

# shellcheck disable=SC2086
for bucketSize in $(seq 50 50 500); do
  params="--numObjects 50M --numQueries 0 --leafMethod rotations --leafSize 10 --bucketSize $bucketSize"
  ./gpurecsplit_construction $params --numThreads 1
  ./simdrecsplit_construction $params --numThreads 1
  ./gpurecsplit_construction $params --numThreads 8
  ./simdrecsplit_construction $params --numThreads 8
done
