#!/bin/bash
# Path: proof_of_concept/scripts/06_compare_to_reference.sh

# check if $1 is provided
if [ $# -eq 0 ]; then
    echo "Usage: $0 folder"
    exit 1
fi

# check if folder exists
folder=../data/$1
if [ ! -d "$folder" ]; then
    echo "Error: folder $folder does not exist"
    exit 1
fi

# read file_name without a last 2 characters
file_name=${1::-2}

# read in reference genome and store in array
# read file ../data/$1/known_systems/Crispr+$1.fna 
# into array ref_genome
# read only the lines that come below >spacer <int>
# and store them in the array
ref_genome=()
# read all lines into array
all_lines=()

while read line; do
    all_lines+=("$line")
done < ../data/$1/reference/Crispr$file_name\_1.fna

while read line; do
    all_lines+=("$line")
done < ../data/$1/reference/Crispr$file_name\_2.fna

# remove all the elements that come before >spacer <int>
for ((i=0; i < ${#all_lines[@]}; i++)); do
    if [[ ${all_lines[$i]} == ">spacer"* ]]; then
        ref_genome=("${all_lines[@]:$i}")
        break
    fi
done

# create a new array and keep only the lines that do not equal to >spacer <int>
new_ref_genome=()
for ((i=0; i < ${#ref_genome[@]}; i++)); do
    if [[ ${ref_genome[$i]} != ">spacer"* ]]; then
        new_ref_genome+=("${ref_genome[$i]}")
    fi
done

# print length of ref_genome
echo "Length of reference genome: ${#new_ref_genome[@]}"
#print ref_genome line by line
for ((i=0; i < ${#new_ref_genome[@]}; i++)); do
    echo ${new_ref_genome[$i]}
done

# read in the cycles from cycles_merged.txt into an array
cycles=()
while read line; do
    cycles+=("$line")
done < $folder/cycles/cycles_merged.txt

# look if all elements in ref_genome appear in cycles array and using a counter
# count how many of them are found
counter=0
for ((i=0; i < ${#new_ref_genome[@]}; i++)); do
    for ((j=0; j < ${#cycles[@]}; j++)); do
        if [[ ${new_ref_genome[$i]} == ${cycles[$j]} ]]; then
            counter=$((counter+1))
            break
        fi
    done
done

# print the number of elements found in cycles array
echo "Number of elements found in cycles array: $counter"