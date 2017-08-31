#!/bin/bash
#
MYPATH="/Users/Moser/Desktop/Dropbox/PhD/Magfield/MagfieldDaniel/"
if [ $1 == '' ] 
	then 
		echo "BASH: Give Folder Name!"
		exit 1
	else MYNAME=$1
fi
#
if [ $2 == '' ]
	then 
		echo "BASH: Give config name!"
		exit 1
	else CONFIG=$2
fi
#
FULLPATH=${MYPATH}${MYNAME}"/"
echo 
echo "BASH: "$FULLPATH
#
if [ -d "$FULLPATH" ] ;
	then
  		echo "BASH: PATH exists!"
	else
		mkdir $FULLPATH
		echo "BASH: New PATH made"
		
		cp ${MYPATH}${CONFIG} $FULLPATH
		./mainmagJ $FULLPATH $CONFIG
fi
#
exit
