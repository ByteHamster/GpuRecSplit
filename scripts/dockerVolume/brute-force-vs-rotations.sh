#!/bin/bash
hostname
strings recsplit_construction | grep " -m" | tee /opt/dockerVolume/brute-force-vs-rotations.txt

# shellcheck disable=SC2086
for leafSize in $(seq 4 2 14); do
    for bucketSize in 100 500 2000; do
        params="--numObjects 1M --numQueries 0 --leafSize $leafSize --bucketSize $bucketSize"
        ./recsplit_construction $params --leafMethod bruteforce | tee --append /opt/dockerVolume/brute-force-vs-rotations.txt
        ./recsplit_construction $params --leafMethod rotations | tee --append /opt/dockerVolume/brute-force-vs-rotations.txt
    done
done

# Build plot
cd /opt/dockerVolume
/opt/sqlplot-tools/build/src/sqlplot-tools brute-force-vs-rotations.tex
rm -f brute-force-vs-rotations.pdf
pdflatex brute-force-vs-rotations.tex
pdflatex brute-force-vs-rotations.tex
rm -f *.out *.log *.aux

