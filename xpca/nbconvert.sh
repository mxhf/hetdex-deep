#/bin/bash

jupyter nbconvert -to python $1 
n=`echo $1 | sed 's/.ipynb/.py/g'`
cat $n | sed 's/COMMANDLINE = False/COMMANDLINE = True/g' > j
mv j $n
