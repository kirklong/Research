#!/bin/bash
#this script takes inlist and finds relevant files
#on Daryl's computer for LC processing later

#code block below takes file separated by , (with a space after comma) and stores each col as
#variables obsid and src respectively
while IFS=$'\ ,' read -r obsid src
do
  filePath="/Users/djm/Public/cxc/smc/${obsid}/primary/${obsid}_src_${src}.times"
  #echo $filePath
  scp klong@claptrap:$filePath ./test
done < "obsid_src_test.txt"

#cd /Users/djm/Public/cxc/smc
#pwd
