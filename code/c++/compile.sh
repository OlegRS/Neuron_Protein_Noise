#!/bin/bash

mkdir build
cd build
echo "Running from $(pwd)"

cmake ..

cmake --build . --target optimised_moments
