#!/bin/bash
mklvars=$(find /opt/intel/* -name mklvars.sh)
source $mklvars intel64
echo $MKLROOT
