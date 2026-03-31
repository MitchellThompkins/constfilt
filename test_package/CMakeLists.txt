cmake_minimum_required(VERSION 3.13)
project(test_package CXX)

find_package(constfilt CONFIG REQUIRED)

add_executable(test_package test_package.cpp)
target_link_libraries(test_package PRIVATE constfilt::constfilt)
