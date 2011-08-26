#!/bin/sh
#
# rna-analyse.sh
#
# Wrapper for RnaChipIntegrator.py which only runs RNA-seq analysis
#
if [ -z "$1" ] ; then
    echo Usage: `basename $0` \[OPTIONS\] \<transcript-data\> \<chip-data\>
    exit
fi
python `dirname $0`/RnaChipIntegrator.py --rna $@
##
#

