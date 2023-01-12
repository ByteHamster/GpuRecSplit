#!/bin/bash
hostname
strings recsplit_construction | grep " -m" | tee figure-5.txt

# shellcheck disable=SC2086
for leafSize in $(seq 2 1 16); do
    for bucketSize in 100 500 1000 2000; do
        params="--numObjects 1M --numQueries 0 --leafSize $leafSize --bucketSize $bucketSize"
        ./recsplit_construction $params --leafMethod bruteforce | tee --append figure-5.txt
        ./recsplit_construction $params --leafMethod rotations | tee --append figure-5.txt
    done
done

for leafSize in $(seq 2 2 34); do
    for bucketSize in 100 500 1000 2000; do
        params="--numObjects 1M --numQueries 0 --leafSize $leafSize --bucketSize $bucketSize"
        ./recsplit_construction $params --leafMethod cuckoo | tee --append figure-5.txt
    done
done

# Build plot
cp figure-5.txt /opt/dockerVolume
cd /opt/dockerVolume
/opt/sqlplot-tools/build/src/sqlplot-tools figure-5.tex
pdflatex figure-5.tex
pdflatex figure-5.tex

