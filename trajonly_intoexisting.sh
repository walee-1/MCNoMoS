#!/bin/bash
#
MYPATH="/home/dmoser/MagSimResults/"
CONFIGPATH="/home/dmoser/Dropbox/PhD/Magfield/MagfieldDaniel/"

if [ $1 == '' ] 
	then 
		echo "BASH: Give Folder Name!"
		exit 1
	else MYNAME=$1
fi

if [ $2 == '' ]
	then 
		echo "BASH: Give config name!"
		exit 1
	else CONFIG=$2
fi

FULLPATH=${MYPATH}${MYNAME}"/"
echo 
echo "BASH: "$FULLPATH

if [ -d "$FULLPATH" ] ;
	then
  		echo "BASH: PATH exists! Write in this Path"
		mkdir $FULLPATH"trajectories"
		./Maintraj $FULLPATH $CONFIG
	else
		mkdir $FULLPATH
		echo "BASH: New PATH made"
		
		cp ${CONFIGPATH}${CONFIG} $FULLPATH
		mkdir $FULLPATH"trajectories"
		./Maintraj $FULLPATH $CONFIG
fi
#
exit
