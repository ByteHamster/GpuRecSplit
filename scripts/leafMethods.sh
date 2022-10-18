#!/bin/bash
hostname
strings recsplit_construction | grep " -m"

# shellcheck disable=SC2086
for leafSize in $(seq 6 2 30); do
    for bucketSize in 100 500 1000 2000; do
        params="--numObjects 1k --numQueries 0 --leafSize $leafSize --bucketSize $bucketSize"

        if [[ "$leafSize" -lt 16 ]]; then
          ./recsplit_construction $params --leafMethod bruteforce
          ./recsplit_construction $params --leafMethod rotations
        fi
        ./recsplit_construction $params --leafMethod cuckoo
    done
done
