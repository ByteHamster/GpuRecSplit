#!/bin/bash
hostname
strings recsplit_construction | grep " -m"
strings gpurecsplit_construction | grep " -m"
strings simdrecsplit_construction | grep " -m"

# shellcheck disable=SC2086
for leafSize in $(seq 5 1 16); do
    for bucketSize in 5 50 500 2000; do
        params="--numObjects 5M --numQueries 0 --leafSize $leafSize --bucketSize $bucketSize --leafMethod rotations"
        ./recsplit_construction $params
        ./gpurecsplit_construction $params
        ./simdrecsplit_construction $params
    done
done
