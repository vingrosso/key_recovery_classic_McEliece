#!/bin/bash

REPS=10
ns=(512 3488 4608 6688 6960 8192)
ts=(20  64   96   128  119  128 )
ms=(9   12   13   13   13   13  )

rm -f time_*.csv
rm -f average_time_*.txt
for i in ${!ns[*]}; do
    echo "${ns[$i]}"
    for rep in $(seq 1 $REPS); do
	date
	echo "$rep"
	python3.8 -O attack.py ${ns[$i]} ${ts[$i]} ${ms[$i]} >> time_${ns[$i]}.csv 2> /dev/null;
    done 
    awk '{ for(i=1;i<=NF;i++) total[i]+=$i ; }
    	 END { for(i=1;i<=NF;i++) printf "%f ",total[i]/NR ;}' time_${ns[$i]}.csv > average_time_${ns[$i]}.txt;
    cat average_time_${ns[$i]}.txt;
    echo ""
done
