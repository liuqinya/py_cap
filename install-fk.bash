#!/bin/bash -x
mkdir -p fk_src
cd fk_src/
# download the fk code from Lupei's website (last tested on Apr 2022)
wget  https://www.eas.slu.edu/People/LZhu/downloads/fk3.3.tar -O fk3.3.tar
tar -xvf fk3.3.tar
cd fk/
# download the patch from SeisMan's website
wget http://www.physics.utoronto.ca/~liuqy/temp/fk-20220331.patch -O fk-20220331.patch
patch < fk-20220331.patch
# compile the code
make
make clean
# test if installation is successful
fk.pl
# test if $SACHOME is set:
#if [ -z ${SACHOME} ]; then 
#   echo '****SACHOME variable needs to be set first in .bashrc or .bash_profile ****'
#else 
#   echo 'SACHOME variable has already been set'
#fi
