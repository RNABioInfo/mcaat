
#convert to python above bash script
import sys
import os
import re
import Bio

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
            
            if record.id.startswith("spacer"):
               
                all_lines.append(str(record.seq))
    print("Length of reference genome before duplicates:", len(all_lines))
    # create a new array and keep only the lines that do not equal to >spacer <int>
    new_ref_genome = []
    ref_gen = set(all_lines)
    new_ref_genome = list(ref_gen)
    # print length of ref_genome
    print("Length of reference genome:", len(new_ref_genome))
    # print ref_genome line by line

    # read in the cycles from cycles_merged.txt into an array
    cycles = set()
    with open(folder + "/cycles/cycles_merged.txt", "r") as f:
        for line in f:
            cycles.add(line.strip())
    cycles=list(cycles)
    # look if all elements in ref_genome appear in cycles array and using a counter
    # count how many of them are found
    print("Elements in cycles:", len(cycles))
    counter = 0
    not_found = []
    marked = dict()
    found_counter = 0
    
    for i in range(len(new_ref_genome)):
        marked[new_ref_genome[i]] = []
        temp = []
        for j in range(len(cycles)):
            if new_ref_genome[i] in cycles[j]:
                temp.append(cycles[j])
            
        marked[new_ref_genome[i]] = temp
    
    #create Counter object for found elements and count their occurence
    

    # write marked dictionary to file
    with open(folder + "/cycles/benchmarks.txt", "w") as f:
        for k, v in marked.items():
            
            f.write(k + ":" + str(len(v)) + "")
            f.write("\n")
            if len(v) != 0:
                counter += 1
            
            for i in range(len(v)):
                f.write("\t"+v[i]+"\n")
                
            f.write("\n")

    values = []
    false_positives = set()
    for l in marked.values():
        for m in l:
            values.append(m)

    
    #check if cycles do not appear in values
    for n in range(len(cycles)):
        if cycles[n] not in values:
            false_positives.add(cycles[n])

    #sort false positives by length
    false_positives = sorted(false_positives, key=len)

    #write false_positives to benchmarks.txt
    with open(folder + "/cycles/benchmarks.txt", "a") as f:
        f.write("\n")
        f.write("False positives:\n")
        for o in false_positives:
            f.write("\t"+o+"\n")

    counter = 0
    print("False positives:",len(false_positives))
    
    for i in marked:
        if len(marked[i]) != 0:
            counter += 1
    # print the number of elements found in cycles array
    
    print("Elements found in cycles array:", counter)
    # if th
if __name__ == "__main__":
    main()