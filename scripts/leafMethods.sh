#!/bin/bash
hostname
strings recsplit_construction | grep " -m"

# shellcheck disable=SC2086
for leafSize in $(seq 2 1 16); do
    for bucketSize in 100 500 1000 2000; do
        params="--numObjects 5M --numQueries 0 --leafSize $leafSize --bucketSize $bucketSize"
        ./recsplit_construction $params --leafMethod bruteforce
        ./recsplit_construction $params --leafMethod rotations
    done
done

for leafSize in $(seq 2 2 34); do
    for bucketSize in 100 500 1000 2000; do
        params="--numObjects 5M --numQueries 0 --leafSize $leafSize --bucketSize $bucketSize"
        ./recsplit_construction $params --leafMethod cuckoo
    done
done
