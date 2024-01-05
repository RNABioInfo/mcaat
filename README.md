# mCAAT

Metagenomic CRISPR Array Analysis Tool

## Description

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


### TODO List
- [x] Create a project
- [x] Create a README.md
- [x] Create a .gitignore
- [x] Create a LICENSE
- [x] Organize CycleFinder
- [x] Find out why writing cycles to file doesn't work
- [ ] Create reading and writing mechanism for the workflow
