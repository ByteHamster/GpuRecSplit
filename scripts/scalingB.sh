#!/bin/bash
hostname
strings gpurecsplit_construction | grep " -m"
strings simdrecsplit_construction | grep " -m"

# shellcheck disable=SC2086
for bucketSize in $(seq 100 100 2000); do
  params="--numObjects 20M --numQueries 0 --leafMethod rotations --leafSize 8 --bucketSize $bucketSize"
  ./gpurecsplit_construction $params --numThreads 1
  ./simdrecsplit_construction $params --numThreads 1
  ./gpurecsplit_construction $params --numThreads 8
  ./simdrecsplit_construction $params --numThreads 8
done
