#!/bin/bash
## need to use gcc44/g++44 of CentOS 5 to compile.

# for Boost
mkdir pathdir/
ln -s /usr/bin/gcc44 pathdir/gcc
export PATH="$(pwd)/pathdir:$PATH"

echo "using gcc :  : /usr/bin/g++44 ; " > ~/user-config.jam

# for NLopt and SAINT
export CXX=g++44
export CC=gcc44

make -j
