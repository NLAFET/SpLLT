#! /bin/sh

git clone https://github.com/ralna/spral.git
cd spral
./autogen.sh
cd ..
mkdir spral-build
cd spral-build
../spral/configure F77=${FC} --disable-gpu --disable-openmp
