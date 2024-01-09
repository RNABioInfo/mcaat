
#convert to python above bash script
import sys
import os

from Bio import SeqIO

def calculate():
    
    if len(sys.argv) != 2:
        print("Usage: python3 05_compare_to_ref.py reference genome id")
        sys.exit(1)


    folder_for_reference = "../data/Genomes/CRISPR_Seq/" + sys.argv[1][:-2] +"/"
    if not os.path.isdir(folder_for_reference):
        print("Error: folder", folder_for_reference, "does not exist")
        sys.exit(1)
        
    file_name = sys.argv[1][:-2]
    print("File name:", file_name)
    spacers = []
    spacers_with_duplicates = []
    for filename in os.listdir(folder_for_reference):
        for record in SeqIO.parse(folder_for_reference + filename, "fasta"):
            if record.id.startswith("spacer"):
               spacers_with_duplicates.append(str(record.seq))
    
    spacers = list(set(spacers_with_duplicates))
    
    print("- Spacers in genome: ", len(spacers))
    
    folder_for_cycles="../data/"+sys.argv[1]
    
    cycles = set()
    with open(folder_for_cycles + "/cycles/cycles_merged.txt", "r") as f:
        for line in f:
            cycles.add(line.strip())
    cycles=list(cycles)
    
    print("- Elements in cycles:", len(cycles))
    counter = 0
    marked = dict()
    
    for i in range(len(spacers)):
        marked[spacers[i]] = []
        temp = []
        for j in range(len(cycles)):
            if spacers[i] in cycles[j]:
                temp.append(cycles[j])
        marked[spacers[i]] = temp
    
    with open(folder_for_cycles + "/cycles/benchmarks.txt", "w") as f:
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

    
    for n in range(len(cycles)):
        if cycles[n] not in values:
            false_positives.add(cycles[n])

    #sort false positives by length
    false_positives = sorted(false_positives, key=len)

    #write false_positives to benchmarks.txt
    with open("../data/"+sys.argv[1] + "/cycles/benchmarks.txt", "a") as f:
        f.write("\n")
        f.write("False positives:\n")
        for o in false_positives:
            f.write("\t"+o+"\n")

    matched_spacer_count = 0
    print("- False positives:",len(false_positives))
    
    for i in marked:
        if len(marked[i]) != 0:
            matched_spacer_count += 1
            
    
    print("- Match:", matched_spacer_count)
    print("- Match percentage:", matched_spacer_count/len(spacers)*100)
    print("- - - - - - - - - - - - - - - - - - - - - - - -")

if __name__ == "__main__":
    calculate()