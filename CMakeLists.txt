cmake_minimum_required (VERSION 2.6)

project (kmerz)

add_executable(kmerz src/main.cpp src/kmerz.cpp src/config.cpp)
add_executable(test_kmerz src/test_kmerz.cpp src/kmerz.cpp src/config.cpp)
