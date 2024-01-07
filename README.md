# mCAAT

Metagenomic CRISPR Array Analysis Tool

## DEVELOPER Description

Current version is considered for verification purposes. Folder ```proof of concept``` is where the magic happens for now. Current workflow is following:
- go to [crispr_cas_db](https://crisprcas.i2bc.paris-saclay.fr/MainDb/StrainList#)
- pick a genome that contains certain amount of crispr arrays
- download fasta file from there to ```proof_of_concept/data/``` folder 
- make sure to name your fasta file as it is named in database
- from ```proof_of_concept/data ``` navigate to ```proof_of_concept/scripts``` 
- call ```./01_organize.sh <genome_id>``` -> creates a folder with genome name and moves fasta file to ```proof_of_concept/data/<genome_id>/genome/<genome_id>.fasta```
- call ```python3 ./07_generate_reads.py <genome_id>``` -> creates perfect reads from the genome and places them to ```proof_of_concept/data/<genome_id>/reads/R1.fastq```
- call ```./03_build_graph.sh <genome_id>``` -> builds a graph and places under ```.../graph/```
- go to main mcaat folder with compiled to software and call ```./mcaat u_mode <genome_id> 100``` - finds all cycles in a user mode with the length upto 100 -> output folder ```proof_of_concept/data/<genome_id>/cycles```:
    - ```id_paths.txt``` -> saves all cycles as vector of integers
    - ```str_paths.txt``` -> saves all cycles as vector of labels
    - ```multiplicity_distribution.txt``` -> saves multiplicities of every node as vector of integers
- navigate to ```proof_of_concept/scripts``` and execute ```./04_merge_cycles.sh <genome_id>```
- download all the crispr_array files from the [crispr_cas_db](https://crisprcas.i2bc.paris-saclay.fr/MainDb/StrainList#) for current genome_id and save them in ```proof_of_concept/data/<genome_id>/reference/``` folder
- execute ```python3 ./05_compare_to_ref.py``` - creates a file called ```benchmarks.txt``` in cycles folder

## USER Description 
### ```this version is a predicted workflow for users in future!```  

This is a repo for our group's tool - mCAAT. This version of the tool is only intended for proof of concept. However here is a sneak peek of what features are coming soon:  
- libs/megahit - using this library we are reading our succinct de Bruijn graph
- src/data - soon to be imlemented feature, input and output files:
    - in/graphs - input graph, currently input graph is statically read from another folder(you can see/modify the structure of folder naming using main.cpp)
    - in/settings.txt - specifies the settings for megahit like number of cores and k-mer size(this feature is to be implemented soon)
    - out - soon to store all the cycles that are supposed to contain all the spacers; current implementation writes everything to following folder: "proof_of_concept/data/"+genome_name+"/cycles/id_paths.txt"
- helpers - formats all the inputs and outputs for above folders
- postprocessing - soon to be created


## Installation

For now all you need to do is either run __make__ on this directory or simply trust the process and run __./mcaat <genome_name> <maximum_length_of_cycle>__.


### TASKS
#### DONE
- [x] Create a .gitignore
- [x] Create a LICENSE
- [x] Organize CycleFinder
- [x] Find out why writing cycles to file doesn't work
- [X] Create reading and writing mechanism for current workflow
#### TODO
- [ ] Look through all the false positives saved in benchmarks and get line id
- [ ] Compare if multiplicity distribution is any different from the true positives!
- [ ] Test bimodal distribution on multiplicities
- [ ] Isolate cycles into clusters
- [ ] Check if clusters start with same IDs 