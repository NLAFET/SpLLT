#! /bin/sh

git clone https://github.com/ralna/spral.git
cd spral
./autogen.sh
cd ..
mkdir spral-build
cd spral-build
../spral/configure FC=gfotran-6 CC=gcc-6 CXX=g++-6 --disable-gpu --disable-openmp
