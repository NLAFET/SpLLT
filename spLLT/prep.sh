#!/bin/bash

CACHE=/numerical/matrices/symmetric/UFL

COLLECTION=`echo $1 | awk -F / '{ print $1 }'`
PRBLM=`echo $1 | awk -F / '{ print $2 }'`

rm matrix.rb
tar -zxvpf $CACHE/$COLLECTION/$PRBLM.tar.gz > /dev/null
cp $PRBLM/$PRBLM.rb matrix.rb
rm -r $PRBLM
