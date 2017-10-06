#! /bin/sh

git clone https://github.com/ralna/spral.git
cd spral
./autogen.sh
cd ..
mkdir spral-build
cd spral-build
~/spral/configure --disable-gpu --disable-openmp
