#!/bin/bash

#Change the name of the file where are the constants
myname="$2"
sed -i "s/^  this->name =.*/  this->name = \"${myname}\";/" src/constants.cpp

#Compilation process
cd build/
cmake ../ -DMY_COMPILE_NAME=$1
make
cd ../

#Excution process
time ./build/$1
rm ./build/$1
