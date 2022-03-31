#!/bin/bash

# This script should be runned from the root of the repository

#######################################
### building gauss-js documentation ###
#######################################

cd gaussjs/ && npm i

./node_modules/.bin/jsdoc --configure ./jsdoc.json --verbose

cd ..

cp -r ./gaussjs/docs ./docs/gaussjs

#######################################
### building gauss documentation ###
#######################################

doxygen gaussdoxy.conf

cd ./docs/gauss/latex && make

cd ../../..

cp ./docs/gauss/latex/refman.pdf ./docs/gauss/docs.pdf

rm -rf ./docs/gauss/latex
