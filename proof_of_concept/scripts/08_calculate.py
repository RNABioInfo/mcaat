
# read in cycle ID file into a list
# read in multiplicities file into a list
# for each of the multiplicities, find the average in a list, minimum and maximum
# write to a file called averages.txt in the same folder as the multiplicities file

# Path: development/git/mCAAT/proof_of_concept/scripts/09_calculate_multiplicities.py

import sys

def multiplicities_average():
    folder = "../data/"+sys.argv[1]+"/cycles/"
    multiplicities = []
    with open(folder + "multiplicity_distribution.txt", "r") as f:
        
        for line in f:
            multiplicities.append(line.split())
    
    

    for t in range(len(multiplicities)):
        multiplicities_renewed = []
        for i in range(len(multiplicities[t])):
            multiplicities_renewed.append(int(multiplicities[t][i]))
        multiplicities[t] = multiplicities_renewed
    multiplicities_ids = associate_multiplicities()
    with open(folder + "averages.txt", "w") as f:
        for i in range(len(multiplicities)):
            base = "Average: " + str(int(sum(multiplicities[i])/len(multiplicities[i]))) + "\t"+ "Minimum: " + str(min(multiplicities[i])) + "\t"+ "Maximum: " + str(max(multiplicities[i])) + "\t"
            if multiplicities_ids[0].get(multiplicities_ids[1][i]) != None:
                f.write(base+"FP\n")
            else:
                f.write(base+"TP\n")
            
    print("Averages written to", folder + "averages.txt")

def associate_multiplicities():
    # we have two groups of inputs: false_positives and matched
    # we need to associate each false positive with a multiplicity value
    # we can do it by checking at which line false positives occurs in the cycles_merged.txt file
    # and then checking the same line in the multiplicity_distribution.txt file
    # file benchmarks.txt contains the false positives in a following format:
    # False positives:
    #    AAATTTTCCCGGGGGG...CCTAATGTTTCCCGGGG
    #    AAATTTTCTTGGGGGG...CCTAATGTTTCCCGGGG


    folder = "../data/"+sys.argv[1]+"/cycles/"
    false_positives = []
    with open(folder + "false_positives.txt", "r") as f:
        for line in f:
            false_positives.append(line.strip())
    
    fp_cycles_ids_multiplicities = {}
    cycles_merged = []
    count = 0
    with open(folder + "cycles_merged.txt", "r") as f:
        for line in f:
            count += 1
            cycles_merged.append(line.strip())
            if line.strip() in false_positives:
                fp_cycles_ids_multiplicities[line.strip()] = count
                
    
    return fp_cycles_ids_multiplicities, cycles_merged
    

    
if __name__ == "__main__":
    multiplicities_average()