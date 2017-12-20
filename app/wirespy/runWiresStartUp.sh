#!/bin/bash
#source runWiresStartUp.sh directory
# S Glavin 12/11/17
# setup: all runs that you want to run wires.py have to be in directory (filename format: tracks_{runnumber}.root)

DIR="$1"
FILES="$DIR/tracks_*.root"

for f in $( ls $FILES ); do
        echo $(date +"%T")
        echo python wires_startup.py $f
	python wires_startup.py $f
done


