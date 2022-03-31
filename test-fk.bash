#!/bin/bash -x
cd fk_src/fk/
# to obtain Greens function hk_15/dist.grn.[0-8]
fk.pl -Mhk/15/k -N512/0.1 10 30 50 70 
# use -S0 to obtain Greens functions for explosive source:
# hk_15/dist.grn.[abc] (the .c file is empty)
fk.pl -Mhk/15/k -N512/0.1 -S0 10 30 50 70 

cd -
# here is a command to obtain a distance list for fk.pl
# data file wildcards need to be modified; sometimes for simplicity
# distances may be rounded to the closest 0s and 5s'.
#dist_list=`saclst dist f data/*.z | sort -g -k 2 | awk '{print $2}' | tr '\n' ' '` 
#fk.pl -Mcus/15/k -N1024/0.05 ${dist_list}
