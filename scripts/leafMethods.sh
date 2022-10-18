#!/bin/bash
hostname
strings recsplit_construction | grep " -m"

# shellcheck disable=SC2086
for leafSize in $(seq 2 1 34); do
    for bucketSize in 50 100 500 1000 1500 2000; do
        params="--numObjects 5M --numQueries 0 --leafSize $leafSize --bucketSize $bucketSize"

        if [[ "$leafSize" -lt 18 ]]; then
          ./recsplit_construction $params --leafMethod bruteforce
          ./recsplit_construction $params --leafMethod rotations
        fi
        ./recsplit_construction $params --leafMethod cuckoo
    done
done
