#!/bin/bash

# Script created by Oriol Borrajo on 20 November 2023
# https://github.com/uriborrajo/HETGEN1000/

## ./exabayes.sh {PATH} *.phylip config.nex
cd "$1"
echo "ENTERING: $1" 
mpirun exabayes -np 4 -R 1 -C 4 -f "$2" -m DNA -c "$3" -n run1 -s 1234 -M 1 #exabayes run1
# Total walltime elapsed: 43:30:4.35 (hh:mm:ss) -M 3-
# Total CPU time elapsed: 1044:01:44.38 (hh:mm:ss) -M 3-
mpirun exabayes -np 4 -R 1 -C 4 -f "$2" -m DNA -c "$3" -n run2 -s 1234 -M 1 #exabayes run2
# Total walltime elapsed: 39:11:50.19 (hh:mm:ss) -M 1-
# Total CPU time elapsed: 940:44:4.51 (hh:mm:ss) -M 1-
mpirun exabayes -np 4 -R 1 -C 4 -f "$2" -m DNA -c "$3" -n run3 -s 1234 -M 1 #exabayes run3
# Total walltime elapsed: 27:06:23.76 (hh:mm:ss) -M 1-
# Total CPU time elapsed: 650:33:30.29 (hh:mm:ss) -M 1-
mpirun exabayes -np 4 -R 1 -C 4 -f "$2" -m DNA -c "$3" -n run4 -s 1234 -M 1 #exabayes run4
# Total walltime elapsed: 25:51:55.34 (hh:mm:ss) -M 1-
# Total CPU time elapsed: 620:46:8.13 (hh:mm:ss) -M 1-
# mpirun exabayes -np 16 -R 4 -C 4 -f "$2" -m DNA -c "$3" -n run1 -s 1234 -M 1 

## ./exabayes.sh {PATH} *.phylip config.nex aln.part
# cd "$1"
# echo "ENTERING: $1" 
# mpirun exabayes -np 4 -R 1 -C 4 -f "$2" -m DNA -c "$3" -q "$4" -n run1 -s 1234 -M 1 #exabayes run1 w/ partition file
# mpirun exabayes -np 4 -R 1 -C 4 -f "$2" -m DNA -c "$3" -q "$4" -n run2 -s 1234 -M 1 #exabayes run2 w/ partition file
# mpirun exabayes -np 4 -R 1 -C 4 -f "$2" -m DNA -c "$3" -q "$4" -n run3 -s 1234 -M 1 #exabayes run3 w/ partition file
# mpirun exabayes -np 4 -R 1 -C 4 -f "$2" -m DNA -c "$3" -q "$4" -n run4 -s 1234 -M 1 #exabayes run4 w/ partition file
# mpirun exabayes -np 16 -R 4 -C 4 -f "$2" -m DNA -c "$3" -q "$4" -n run1 -s 1234 -M 1
