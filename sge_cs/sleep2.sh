#!/bin/bash
#$ -S /bin/bash
#$ -wd /home/uctpln0/FruitDemand/code/fortran/sge_cs
#$ -o /home/uctpln0/FruitDemand/code/fortran/output/sleep2.out
#$ -N plotA27A
#$ -j y
#$ -R y
#$ -cwd
#$ -V
#$ -t 1-101
#$ -l tmem=2.4G,h_vmem=2.4G
##$ -l h_rt=240:0:0
##$ -l hostname="fry*|burns*|larry*|zeppo*"
##$ -l hostname=fry*\|burns-*\|zeppo*\|larry*
##$ -pe orte 401
##  -terse  -R y -t 1-$np -V $res_list 

cd /home/uctpln0/FruitDemand/code/fortran/source
date
./sleep.exe
date
