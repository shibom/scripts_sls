#!/bin/bash
. /etc/scripts/mx_sls.sh

hklfile=../24xturns.ahkl

reso="3.5 3.6 3.8 4.0 4.2 4.5 4.8 5.0 5.5"
#reso="4.0 4.2 4.4 4.6 4.8 5.0 5.5 6.0"
emins="1.2 1.3 1.4 1.5"
ntrials="5000 10000"

if [! -e shelx_grid]
  then
      mkdir shelx_grid_1
      cd shelx_grid_1
else
    num_folder=`ls shelx_g*`
    mkdir shelx_grid_$num_folder
    cd shelx_grid_$num_folder
fi

echo "#!/bin/bash" >> b1.sh
echo "python ../shelx_batch.py $1 $2 $3 $4 $5" >> b1.sh
chmod +x b1.sh

for i in {60..120..20}
 do
   for res in $reso
    do
      for val in $emins
       do
         for try in $ntrials
          do
            sbatch -J tr-$i-$res-$val ./b1.sh $hklfile $i $res $val $try
          done
      done
   done
done
