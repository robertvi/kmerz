#!/bin/sh

#
# build the main binary and the test binary
# then run the test binary and report "okay" if all test were passed
#
# you can use the Makefile if you only want to build the main binary (kmerz)
#

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
