#!/bin/sh
#
# chip-analyse.sh
#
# Wrapper for RnaChipIntegrator.py which only runs ChIP-seq analysis
#
if [ -z "$1" ] ; then
    echo Usage: `basename $0` \[OPTIONS\] \<gene-data\> \<chip-data\>
    exit
fi
python `dirname $0`/RnaChipIntegrator.py --chip $@
##
#

