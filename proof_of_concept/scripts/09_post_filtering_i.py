import sys
from multiprocessing import Pool
import math
from collections import Counter
import statistics

#import matplotlib.pyplot as plt
#import numpy as np
def read_file(folder, file_name):
    output_list = []
    with open(folder + file_name, "r") as f:
        output_list = [[int(n) for n in line.split()] for line in f]
    return output_list

def check_if_subset(list1, list2):
    if len(list1) != len(list2):
        if set(list1).issubset(set(list2)):
            return 1
    return 0

def rearrange_cyles(id_paths,multiplicities,indices):
    new_id_paths = []
    new_multiplicities=[]
    for i in range(len(id_paths)):
        start = indices[i] 
        if start == 0:
            new_id_paths.append(id_paths[i])
            new_multiplicities.append(multiplicities[i])
            continue
        end = -1
        elements_before_index = id_paths[i][0:start]
        path_to_add = id_paths[i][start:end]
        path_to_add.extend(elements_before_index)
        new_id_paths.append(path_to_add)
        
        multiplicities_before_index = multiplicities[i][0:start]
        multiplicities_to_add = multiplicities[i][start:end]
        multiplicities_to_add.extend(multiplicities_before_index)
        new_multiplicities.append(multiplicities_to_add)
        
    return new_id_paths, new_multiplicities

def remove_false_positives(files=""):
    folder = "/vol/d/development/git/mCAAT/proof_of_concept/data/"+sys.argv[1]+"/cycles_genome/"
    histogram_folder = "/vol/d/development/git/mCAAT/proof_of_concept/data/"+sys.argv[1]+"/histograms/"
    if files!="":
        folder = "/vol/d/development/git/mCAAT/proof_of_concept/data/"+files+"/cycles_genome/"
        #histogram_folder = files[1]
    multiplicities = read_file(folder, "multiplicity_distribution.txt")
    id_paths= read_file(folder, "id_paths.txt")
    report = ""
    report += "File name: " + files +"\n"
    report += "Amount of cycles before filtering: " + str(len(id_paths)) + "\n"
    indicies = []
    for i in range(len(multiplicities)):
        indicies.append(multiplicities[i].index(max(multiplicities[i])))
    wrapper = rearrange_cyles(id_paths,multiplicities,indicies)
    new_id_paths = wrapper[0]
    new_multiplicities = wrapper[1]
    
    t=0
    id_paths_no_duplicates = []
    multiplicities_no_duplicates = []
    for i in range(len(new_id_paths)):
        if new_id_paths[i] not in id_paths_no_duplicates:
            id_paths_no_duplicates.append(new_id_paths[i])
            multiplicities_no_duplicates.append(new_multiplicities[i])
    

    count = 0
    association = {}
    absolute_min = {}
    node_cluster = {}
    old_min=min(multiplicities_no_duplicates[0])
    for i in range(len(id_paths_no_duplicates)):
        association[id_paths_no_duplicates[i][0]] = []        
        node_cluster[id_paths_no_duplicates[i][0]] = []

    for i in range(len(id_paths_no_duplicates)):
        association[id_paths_no_duplicates[i][0]].append(i)
        absolute_min[id_paths_no_duplicates[i][0]]=(min(multiplicities_no_duplicates[i]) if min(multiplicities_no_duplicates[i])<old_min else old_min)
        old_min=min(multiplicities_no_duplicates[i])
        node_cluster[id_paths_no_duplicates[i][0]].append(multiplicities_no_duplicates[i])
    id_paths = []
    counting_list = {}
    min_max_proportion = []
    g_cs = {}
    
    for k,v in node_cluster.items():
        extended_node_cluster = []

        for i in v:
            extended_node_cluster.extend(i)    
        g_cs[k]=extended_node_cluster
    
    
    count=0
    for i in range(len(multiplicities_no_duplicates)):
        min_multiplicity = min(multiplicities_no_duplicates[i])
        max_multiplicity = max(multiplicities_no_duplicates[i])
        min_max_proportion.append(max_multiplicity/min_multiplicity)
    

    minimum = min(min_max_proportion)
    #print("Minimum:",minimum)
    multiplicities=[]
    c = 0
    for k,v in association.items():
        
        if len(v) < 300 :
            for i in v:
                multiplicities.append(multiplicities_no_duplicates[i])
                id_paths.append(id_paths_no_duplicates[i])
                
            counting_list[k] = len(v)
    
    with open(folder + "min_max_proportion.txt", "w") as f:
        for i in min_max_proportion:
            f.write(str(i)+"\n")
    id_paths_after_filter = []
    
    for i in id_paths:
        flag = False
        for j in id_paths_no_duplicates:
            if check_if_subset(i, j) == 1:
                flag = True
                break
            else:
                flag = False
        count += 1
        if flag:
            continue
        id_paths_after_filter.append(i)
        if count % 200 == 0:
            print("|",end=" ... ",flush=True)
    
    print("\nPost-filtering:",len(id_paths))
    original_len_ids = len(id_paths_no_duplicates)
    t_count = 0
    lengths = {}
    for k,v in association.items():
        lengths[k] = []
        for i in v:
            lengths[k].append(len(multiplicities_no_duplicates[i]))
    
    lengths_average = {}
    for k,v in lengths.items():
        lengths_average[k] = statistics.median(v)
    new_id_paths = []
    '''
    for k,v in lengths_average.items():
        for i in association[k]:
            if len(id_paths_no_duplicates[i]) <= v-8 or len(id_paths_no_duplicates[i]) >= v+8 or len(id_paths_no_duplicates[i]) > 120 or len(id_paths_no_duplicates[i]) < 46:
                continue    
            new_id_paths.append(id_paths_no_duplicates[i])
    '''
    
    print("Post-post-filtering:",len(id_paths_no_duplicates))
    
    with open(folder + "lengths.txt", "w") as f:
        for k,v in lengths.items():
            f.write(str(k) + ": " + str(v) + "\n")
    
    with open(folder + "counting_list.txt", "w") as f:
        for k,v in counting_list.items():
            f.write(str(k) + ": " + str(v) + "\n")
            
    with open(folder + "id_paths_no_duplicates.txt", "w") as f:
        for cycle in id_paths:
            for element in cycle:
                f.write(str(element) + " ")
            f.write("\n")
    
    #report += line + "\n"
    report += "Post-filtering: " + str(len(id_paths_no_duplicates)) + "\n"
    report += "Post-post-filtering: " + str(len(new_id_paths)) + "\n"
    report += "\n----------------------------------------------------\n"
    with open("report.txt", "a") as f:
        f.write(report)
    print()
#read from ../data/data.txt]
lines = []

with open("/vol/d/development/git/mCAAT/proof_of_concept/data/data.txt", "r") as f:
    lines = f.readlines()
    lines = [line.strip() for line in lines]
print(lines[5])


report = ""
remove_false_positives(sys.argv[1])