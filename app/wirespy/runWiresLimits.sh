#!/bin/bash
#source runWiresLimits.sh directory limit
# S Glavin 12/19/17
# setup: all runs that you want to run wires.py have to be in directory (filename format: tracks_{runnumber}.root)

DIR="$1"
LIMIT=$2
FILES="$DIR/tracks_*.root"

for f in $( ls $FILES ); do
	echo $(date +"%T")
        echo python wires_enactLimits.py $f $LIMIT
	python wires_enactLimits.py $f $LIMIT
done


