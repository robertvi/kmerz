#!/bin/sh

set -eu

g++ -g -o kmerz \
./src/kmerz.cpp \
./src/config.cpp \
./src/main.cpp

g++ -g -o test_kmerz \
./src/kmerz.cpp \
./src/config.cpp \
./src/test_kmerz.cpp

./test_kmerz && echo okay
