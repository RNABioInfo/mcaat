#!/bin/bash
echo "--------START ORGANIZE--------"
echo ""
#check if $1 is provided
if [ $# -eq 0 ]; then
    echo "Usage: $0 folder"
    exit 1
fi

file="../data/Genomes/$1.fasta"
#check if filde exists
if [ ! -f "$file" ]; then
    echo "Error: file $file does not exist"
    exit 1
fi



mkdir ../data/$1
mkdir ../data/$1/genome
mkdir ../data/$1/reads
mkdir ../data/$1/graph
mkdir ../data/$1/cycles
mkdir ../data/$1/reference

mv ../data/Genomes/$1.fasta ../data/$1/genome
echo "Folder structure for genome $1 is sorted!"
echo ""
echo "--------END ORGANIZE---------"
exit 0