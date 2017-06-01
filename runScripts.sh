#!/bin/bash


scriptDir=/scratch/d/danfldhs/shunjan4/group_project/pipelineScripts/normal  #output directory for created scripts

#run scripts created by makeScripts.sh


for script in $scriptDir/*  #for all scripts in this directory
do
echo "$script" #test which scripts are running
chmod +x $script  #change execution permission for script
qsub $script  #qsub the script

done

