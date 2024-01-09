#!/bin/bash

#echo "--------START GRAPH BUILD--------"
#echo ""

#echo "Running Megahit on $1"

folder="../data/$1"

fastaFileName=$(ls $folder/genome/*.fasta | xargs -n 1 basename | sed 's/.fasta//g')

fastaFilePath="$folder/genome/$fastaFileName.fasta"

# create data.lib file in  reads folder
#echo "Step 1: Creating data.lib file"

echo "# simulation of $fastaFileName" > $folder/reads/data.lib
echo "se $folder/reads/R1.fastq " >> $folder/reads/data.lib

sleep 1
#echo "Step 2: Building graph"
./megahit_core buildlib $folder/reads/data.lib $folder/graph/outfile_prefix
./megahit_core read2sdbg --host_mem 128000000000 --read_lib_file $folder/graph/outfile_prefix -m 1 -k 23 --num_cpu_threads 24 --o $folder/graph/graph
# 
#echo ""
#echo "--------END GRAPH BUILD---------"

exit 0
