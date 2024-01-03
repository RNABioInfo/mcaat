#!/bin/bash

echo "--------START GENOME2READS--------"
echo ""
# assign folder argument to variable - appropriate folder name is data/$1
folder="../data/$1"

# check if folder argument is provided
if [ $# -eq 0 ]; then
    echo "Usage: $0 folder"
    exit 1
fi

# check if folder exists
if [ ! -d "$folder" ]; then
    echo "Error: folder $folder does not exist"
    exit 1
fi

fastaFileName=$(ls $folder/genome/*.fasta | xargs -n 1 basename | sed 's/.fasta//g')

fastaFilePath="$folder/genome/$fastaFileName.fasta"

coverageFilePath="$folder/genome/coverage.txt"

# check if genomes folder contains fasta files
if [ ! -f "$folder/genome/coverage.txt" ]; then
    #create coverage.txt and write name of fasta file in folder
    #select only fasta file name without extension
   
    echo "$fastaFileName  25" > $coverageFilePath
    echo "Created and written coverage.txt"
    
fi
sleep 1
echo "Executing InsilicoSeq on $fastaFileName with coverage 5"

# loop through files in folder
echo "Processing file: $fastaFileName"
iss generate -g $fastaFilePath --coverage_file $coverageFilePath --mode basic --output $folder/reads/sim  --cpus 20
echo ""
echo "---------END GENOME2READS---------"

exit 0