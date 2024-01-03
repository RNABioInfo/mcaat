import sys
import os
import re
import Bio
from collections import Counter
from Bio import SeqIO

def main():
    # check if $1 is provided
    if len(sys.argv) != 2:
        print("Usage: python3 05_compare_to_ref.py reference genome id")
        sys.exit(1)

    # check if folder exists
    folder = "../data/" + sys.argv[1]
    if not os.path.isdir(folder):
        print("Error: folder", folder, "does not exist")
        sys.exit(1)

    # read file_name without a last 2 characters
    file_name = sys.argv[1][:-2]

    ref_genome = []
    all_lines = []
    for filename in os.listdir("../data/" + sys.argv[1] + "/reference/"):
        #use seqIO to read in fasta file
        #change to ../data/ + sys.argv[1] + /reference/ + filename +"fna" to ../data/ + sys.argv[1] + /reference/ + filename +"fasta"
        for record in SeqIO.parse("../data/" + sys.argv[1] + "/reference/" + filename, "fasta"):
            ref_genome.append(str(record.seq))
            all_lines.append(str(record.seq))
    print("Length of reference genome:", len(ref_genome))
    # create a new array and keep only the lines that do not equal to >spacer <int>
    new_ref_genome = []
    ref_gen = set(all_lines)
    new_ref_genome = list(ref_gen)
    # print length of ref_genome
    print("Length of reference genome:", len(new_ref_genome))
    # print ref_genome line by line

    # read in the cycles from cycles_merged.txt into an array
    cycles = []
    with open(folder + "/cycles/cycles_merged.txt", "r") as f:
        for line in f:
            if line.startswith("\n"):
                continue

            cycles.append(line.strip())
    #assign every 23mer a counter of how many times it appears in the cycles array
    counter = dict()
    k_mer = 23
    for line in cycles:
        count = 0
        key = line[0:k_mer]
        if key in counter:
            continue
        for line_one in cycles:
            if line==line_one:
                continue
            for i in range(0, len(line_one)):
                if line[0:k_mer] == line_one[i:i+k_mer]:
                    count += 1
                    
        counter[key] = count
    
    # write the counter dictionary into a file
    with open(folder + "/cycles/counter.txt", "w") as f:
        for key, value in counter.items():
            f.write('%s:%s\n' % (key, value))
    
    # remove all lines from cycles with value 0
    # using above counter dictionary
    cycles_new = []
    for line in cycles:
        if counter[line[0:k_mer]] >= 1:
            cycles_new.append(line)
    # write the new cycles array into a file
    with open(folder + "/cycles/cycles_new.txt", "w") as f:
        for line in cycles_new:
            f.write(line + "\n")
    
    
    

# execute main function
if __name__ == "__main__":
    main()