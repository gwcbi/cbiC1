#!/bin/sh
### monitor a tab-outputed blast job by giving the approximative % done
#to use: blast_monitor.sh <blastresult.blast> <queryfasta.fa>
blast=$1
query=$2
echo "the blast out is: "$blast
echo "the fasta query is: "$query
echo
curquery=$(tail -1 $blast | cut -f 1)
curline=$(fgrep -n $curquery $query |  cut -f 1 -d ':')
nblines=$(wc -l $query | cut -f 1 -d " ")
percent=$(echo "($curline/$nblines) *100" | bc -l | cut -c 1-4)
echo "The blast job is about $percent % done..."
