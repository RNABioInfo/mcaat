#!/bin/bash
# file proof_of_concept/data/Genomes/genomes contains names of genomes
while read line; do
    echo "For genome: $line"
    ./mcaat u_mode $line 100
done < proof_of_concept/data/Genomes/genomes