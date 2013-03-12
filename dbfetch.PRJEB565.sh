#!/bin/bash

# Mon Mar 11 21:43:32 CDT 2013
# download sequences HAAB01000001-HAAB01090174 from EBI database

for ((i=1; i<=90174; i=i+200)); do
	# final command looks like this: 
	# wget -O cs_trans.PRJEB565.dbfetch.fasta "http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=embl&id=HAAB01*&format=fasta&style=raw"
	# assign the printf output to a given variable with the -v argument
	printf -v cmd1 "wget -O HAAB01%06d-HAAB01%06d.fasta \"http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=EMBL&id=" $i $((i+199))
	for ((k=1; k<=200 && k<=$((90174-i+1)); k++ )); do
		if ((k>1))
		then
			printf -v cmd2 ","
			cmd1=$cmd1$cmd2
		fi
		printf -v cmd3 "HAAB01%06d" $((i+k-1))
		cmd1=$cmd1$cmd3
	done
	printf -v cmd4 "&format=fasta&style=raw\""
	cmd1=$cmd1$cmd4

	# run command
	# echo $cmd1
	eval $cmd1
	sleep 5
done
