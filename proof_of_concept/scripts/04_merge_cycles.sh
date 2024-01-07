#!/bin/bash

# running statistics on the graph
folder=../data/$1
file= $folder/cycles/str_paths.txt
echo "--------START MERGING--------"


# read file line by line
while read line; do
    # if line starts with F, skip
    if [[ $line == F* ]]; then
        continue
    fi
    # split line into array
    words=($line)
    # get first word
    first_word=${words[0]}
    
    # get last character of all other words
    last_chars=""
    for ((i=1; i < ${#words[@]}; i++)); do
        last_chars="$last_chars${words[$i]: -1}"
    done
    # write to file
    echo "$first_word$last_chars" >> $folder/cycles/cycles_merged.txt

done < $folder/cycles/str_paths.txt

# output number of occurences and length of cycles
cat $folder/cycles/cycles_merged.txt | awk '{print length, $0}' | sort -n | uniq -c | sort -rn | head -n 10 > $folder/cycles/statistics.txt

echo ""
echo "--------END MERGING---------"
exit 0
