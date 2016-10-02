#!/bin/bash                                                                                                                

in_dataset="CCVars.idx"
out_dataset="CCVars.gidx"
filename="${in_dataset%.*}"

echo "<?xml version='1.0' encoding='utf-8'?>" > $out_dataset
echo "<datasets name=\"$filename\">" >> $out_dataset

for d in */ ; do
    if [[ $d == "t"* ]]
    then
        for t in $(ls $d/l0/$filename/) ; do

            strtime=`echo "${d:1:5}" | bc`;

            echo "\t<dataset name=\"${d:0:5}\" log_time=\"$strtime\" url=\"file://$PWD/$d/l0/$in_dataset\" />" >> $out_dataset

            echo "adding timestep $d/l0/$in_dataset to ./$out_dataset";
        done
    fi
done

echo "</datasets>" >> $out_dataset

