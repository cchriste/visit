#!/bin/bash                                                                                                                

dataset="CCVars.idx"
filename="${dataset%.*}"

if [ -d $filename ]; then
  #rm -R $filename
  echo "Please remove existing directory $filename if this is an old generated directory"
  exit 1
fi

echo "creating dir $filename"
mkdir -p $filename

last_idx=""

for d in */ ; do
    if [[ $d == "t"* ]]
    then
        last_idx=$d/l0/
        #echo "found $last_idx"                                                                                            
    fi
done

echo "copying $last_idx/$dataset as global idx"
cp $last_idx/$dataset ./
cp $last_idx/$filename"_OFFSET" ./
cp $last_idx/$filename"_SIZE" ./

for d in */ ; do
    if [[ $d == "t"* ]]
    then
        for t in $(ls $d/l0/$filename/) ; do
            echo "adding timestep $d/l0/$filename/$t to ./$filename/";
            ln -s ../$d/l0/$filename/$t $filename/$t
        done
    fi
done

