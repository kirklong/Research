#!/bin/bash
#source /usr/bin/xpaget
#xpaget xpans
echo "starting ds9..."
IMGFILE=$1
REGIONFILE=$2
if [ `xpaaccess ds9` = no ]; then
  ds9&

i=1
while [ "$i" -le 30 ]
do
sleep 2
if [ `xpaaccess ds9` = yes ];
then
  echo "access established -- converting region file"
  break
fi

i=`expr $i + 1`
done
fi
xpaset -p ds9 fits $IMGFILE #THIS COMMAND DOESN"T WORK FROM SHELL???
xpaset -p ds9 scale log
xpaset -p ds9 regions load $REGIONFILE
xpaget ds9 regions -format ds9 -system physical > regions.reg
xpaset -p ds9 exit
#'xpaget ds9 regions -format ds9 -system physical > regions.reg'
#xpaset -p ds9 exit
