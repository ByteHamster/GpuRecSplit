#!/bin/bash

./recsplit_construction     --numObjects 1M --numQueries 50M --leafMethod bruteforce --leafSize 16 --bucketSize 2000
./simdrecsplit_construction --numObjects 1M --numQueries 50M --leafMethod rotations  --leafSize 16 --bucketSize 2000 --numThreads 1
./simdrecsplit_construction --numObjects 1M --numQueries 50M --leafMethod rotations  --leafSize 16 --bucketSize 2000 --numThreads 16
./gpurecsplit_construction  --numObjects 1M --numQueries 50M --leafMethod rotations  --leafSize 16 --bucketSize 2000 --numThreads 4

./recsplit_construction     --numObjects 1M --numQueries 50M --leafMethod bruteforce --leafSize 18 --bucketSize 50
./simdrecsplit_construction --numObjects 1M --numQueries 50M --leafMethod rotations  --leafSize 18 --bucketSize 50 --numThreads 1
./simdrecsplit_construction --numObjects 1M --numQueries 50M --leafMethod rotations  --leafSize 18 --bucketSize 50 --numThreads 16
./gpurecsplit_construction  --numObjects 1M --numQueries 50M --leafMethod rotations  --leafSize 18 --bucketSize 50 --numThreads 4

./gpurecsplit_construction  --numObjects 1M --numQueries 50M --leafMethod rotations  --leafSize 21 --bucketSize 2000 --numThreads 4
