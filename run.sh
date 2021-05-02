# !/bin/bash

cd ./build/tests/
make 
timeout 0.005 ./PlaygroundTest
