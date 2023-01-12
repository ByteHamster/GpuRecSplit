#!/bin/bash
# shellcheck disable=SC2086
hostname
strings simdrecsplit_construction | grep " -m" | tee figure-7.txt

for threads in $(seq 1 2 16); do
  if [[ $threads -eq 0 ]]; then
    threads=1
  fi

  params="--numQueries 0 --leafMethod rotations --leafSize 5 --bucketSize 5 --numThreads $threads"
  ./simdrecsplit_construction $params --numObjects $((80*threads))M | tee --append figure-7.txt

  params="--numQueries 0 --leafMethod rotations --leafSize 8 --bucketSize 100 --numThreads $threads"
  ./simdrecsplit_construction $params --numObjects $((30*threads))M | tee --append figure-7.txt

  params="--numQueries 0 --leafMethod rotations --leafSize 12 --bucketSize 9 --numThreads $threads"
  ./simdrecsplit_construction $params --numObjects $((10*threads))M | tee --append figure-7.txt

  params="--numQueries 0 --leafMethod rotations --leafSize 16 --bucketSize 2000 --numThreads $threads"
  ./simdrecsplit_construction $params --numObjects $((20*threads))k | tee --append figure-7.txt
done

# Build plot
cp figure-7.txt /opt/dockerVolume
cd /opt/dockerVolume
/opt/sqlplot-tools/build/src/sqlplot-tools figure-7.tex
pdflatex figure-7.tex
pdflatex figure-7.tex

