#!/bin/bash
. /etc/scripts/mx_sls.sh

reso="3.5 3.6 3.8 4.0 4.2 4.5 4.8 5.0 5.5"
#reso="4.0 4.2 4.4 4.6 4.8 5.0 5.5 6.0"
emins="1.2 1.3 1.4 1.5"

for i in {60..120..20}
 do
   for res in $reso
    do
      for val in $emins
       do
         sbatch -J tr-$i-$res-$val b1.sh $i $res $val 
      done
   done
done

