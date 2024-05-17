import os
from PROD_FILTER import *
from string_based_filter import *
#from merge import *
def list_subfolders_in_folder(folder):
    # Input: None
    # Output: list of subfolders in folder
    return [f.path[-10:] for f in os.scandir(folder) if f.is_dir() and not f.name.startswith("_")]

folders = list_subfolders_in_folder("../data/")
folder = sys.argv[2]
''''
for folder in folders:
    if "Genomes" in folder:
        continue 
    merge(folder)
'''
if sys.argv[1] == "1":
    id_paths = read_file("../data/"+folder+"/cycles_genome/","id_paths.txt")
    id_multiplicities = read_file("../data/"+folder+"/cycles_genome/","multiplicity_distribution.txt")
    separate_list_by_first_element(id_paths,id_multiplicities,folder)
    print("DONE FOR "+folder+"!")
else:
    with open("log.txt", "a") as f:
        print("Starting with "+folder)
        f.write("Starting with "+folder+"!\n")
        cycles = read("/vol/d/development/git/mCAAT/proof_of_concept/data/"+folder+"/cycles_genome/results.txt")
        sum_of_quantity = [len(cycles[k_mer]) for k_mer in cycles.keys()]
        f.write("Elements before filter: "+str(sum(sum_of_quantity))+"\n")
        #print first 2 elements in mults
        count = 0
        limits = 55
        filtered_cycles = remove_all_with_difference_less_than_23(cycles,limits)
        f.write("Elements after first round of filtering:"+str(count_elements_in_dictionary(filtered_cycles))+"\n")
        f_c = keep_only_those_with_score_greater_than_0(filtered_cycles)
        write_to_file(f_c, "/vol/d/development/git/mCAAT/proof_of_concept/data/"+folder+"/cycles_genome/results_filtered.txt")
        f.write("Elements after second round of filtering:"+str(count_elements_in_dictionary(f_c))+"\n")
        f.write("Elements removed:"+str(count_elements_in_dictionary(cycles)-count_elements_in_dictionary(f_c))+"\n")


#print all folders to file
'''
with open("folders.txt", "w") as f:
    for folder in folders:
        f.write(folder+"\n")
'''