#!/bin/bash
#$ -S /bin/bash
#$ -wd /home/uctpln0/FruitDemand/code/fortran/source
#$ -o /home/uctpln0/FruitDemand/code/fortran/output/sleep1.out
#$ -N A27
#$ -j y
#$ -R y
#$ -cwd
#$ -V
#$ -t 1-400
#$ -l tmem=1.8G,h_vmem=1.8G
#$ -l hostname="burns*|zeppo*|fry*|larry*|cheech*|hale*"
##$ -l h_rt=240:0:0
##$ -l hostname="fry*|burns*|larry*|zeppo*"
##$ -l hostname=fry*\|burns-*\|zeppo*\|larry*
##$ -pe smp 2
##  -terse  -R y -t 1-$np -V $res_list 

cd /home/uctpln0/FruitDemand/code/fortran/source
date
./sleep.exe
date
