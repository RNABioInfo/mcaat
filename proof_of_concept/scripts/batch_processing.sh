#!/bin/bash
# pick every file from data/Genomes folder and run ./01_organize.sh on it
#for file in ../data/Genomes/*.fasta

#do
    # filename is the name without extension = CP067336.1.fasta 
    # keep all characters before the last dot
    # filename=$(basename -- "$file" | cut -d. -f1).1
    # echo "Filename: $filename"
    #run ./01_organize.sh on file
    #./01_organize.sh $filename
#done

# pick every line from data/Genomes/genomes file 
# and run "python3 ./02_generate_reads.py $line"
#while read line; do
#    echo "For genome: $line"
#    python3 ./02_generate_reads.py $line
#done < ../data/Genomes/genomes 

# pick every line from data/Genomes/genomes file
# and run "03_build_graph.sh $line"
while read line; do
    python3 ./05_compare_to_ref.py $line 
done < ../data/Genomes/genome_
