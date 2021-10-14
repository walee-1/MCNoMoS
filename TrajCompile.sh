#!/bin/bash


#SBATCH --job-name = TrajCompile
#SBATCH --output = CompileOut.txt

#SBATCH --ntasks = 1
#SBATCH --time = 00:00:30

EXEDIR="/users/waleed.khalid/"
PROGRAMFOLDER="TrajFiles_v4_NoTZ"
echo ${EXEDIR}${PROGRAMFOLDER}

#if [ $1 == '' ]
#	then
#		echo "no argument given to bash!"
#		exit 1
#	else PROGRAMFOLDER=$1
#fi
#
#
srun hostname

srun make all

