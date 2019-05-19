#!/bin/bash

cd /Volumes/Yolanda/CRF_Screen/Ref

awk -F "\t" '{if \
( (NR==1) || (($10 == 9606) && ($11 == 9606)) || (($10 == 10090) && ($11 == 10090)) || (($10 == 10116) && ($11 == 10116)) )\
{print $0} }' \
BIOGRID-ALL-3.5.172_simplify.tab.txt \
> BIOGRID.use.txt