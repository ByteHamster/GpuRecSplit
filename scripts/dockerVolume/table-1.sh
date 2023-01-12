#!/bin/bash

./recsplit_construction     --numObjects 1M --numQueries 50M --leafMethod bruteforce --leafSize 16 --bucketSize 2000 | tee table-1.txt
./simdrecsplit_construction --numObjects 1M --numQueries 50M --leafMethod rotations  --leafSize 16 --bucketSize 2000 --numThreads 1 | tee --append table-1.txt
./simdrecsplit_construction --numObjects 1M --numQueries 50M --leafMethod rotations  --leafSize 16 --bucketSize 2000 --numThreads 16 | tee --append table-1.txt

./recsplit_construction     --numObjects 1M --numQueries 50M --leafMethod bruteforce --leafSize 18 --bucketSize 50 | tee --append table-1.txt
./simdrecsplit_construction --numObjects 1M --numQueries 50M --leafMethod rotations  --leafSize 18 --bucketSize 50 --numThreads 1 | tee --append table-1.txt
./simdrecsplit_construction --numObjects 1M --numQueries 50M --leafMethod rotations  --leafSize 18 --bucketSize 50 --numThreads 16 | tee --append table-1.txt

# Build plot
cp table-1.txt /opt/dockerVolume
cd /opt/dockerVolume
/opt/sqlplot-tools/build/src/sqlplot-tools table-1.tex
pdflatex table-1.tex
pdflatex table-1.tex

