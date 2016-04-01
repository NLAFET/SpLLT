#!/bin/bash

trace=$1

# remove arrows: 18 and 19
# remove number of submitted/ready tasks: 13
sed -i '/^18\|^19\|^13/d' $trace
 
