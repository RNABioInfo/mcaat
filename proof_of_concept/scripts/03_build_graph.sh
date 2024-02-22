#!/bin/bash

#echo "--------START GRAPH BUILD--------"
#echo ""

#echo "Running Megahit on $1"

folder="../data/$1"

fastaFileName=$(ls $folder/genome/*.fasta | xargs -n 1 basename | sed 's/.fasta//g')

fastaFilePath="$folder/genome/$fastaFileName.fasta"

# create data.lib file in  reads folder
#echo "Step 1: Creating data.lib file"

echo "# simulation of $fastaFileName" > $folder/graph_fasta/$fastaFileName.fasta.lib
echo "se $fastaFilePath" >> $folder/graph_fasta/$fastaFileName.fasta.lib

echo "$fastaFilePath"
sleep 0.5
#echo "Building graph"
./megahit_core buildlib $folder/graph_fasta/$fastaFileName.fasta.lib $fastaFilePath 
#./megahit_core read2sdbg --host_mem 128000000000 --read_lib_file $folder/graphq_fasta/outfile_prefix -m 1 -k 23 --num_cpu_threads 24 --o $folder/graph/graph
./megahit_core seq2sdbg --host_mem 128000000000  --contig $fastaFilePath  -k 23 --num_cpu_threads 24 --o $folder/graph_fasta/graph
# 
#echo ""
#echo "--------END GRAPH BUILD---------"

exit 0
