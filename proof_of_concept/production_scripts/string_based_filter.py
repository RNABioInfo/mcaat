import sys
def read(file_name_cycles):
    group_results = {}
    key = ""
    value = []
    keys = []
    with open(file_name_cycles) as cf:
        for line in cf:
            if len(line.strip())<=23:
                key = line.strip()
                value = []
                continue
            value.append(line.strip())
            group_results[key] = value
    return group_results

def shift_right_common_word_length(cycles_per_cluster):
    if len(cycles_per_cluster) == 1:
        return 0
    k = 23
    size = len(cycles_per_cluster[0])
    count = 0
    for cycle_length_index in range(k+1, size-k):
        for cycle_index in range(0, len(cycles_per_cluster)):
            if cycles_per_cluster[0][:cycle_length_index] in cycles_per_cluster[cycle_index]:
                continue
            else:
                return count+k
        count+=1
    return count+k
        
def shift_left_common_word_length(cycles_per_cluster):
    if len(cycles_per_cluster) == 1:
        return 0
    k = 23
    size = len(cycles_per_cluster[0])
    count = 0
    for cycle_length_index in range(size-k-1,k+1,-1):
        for cycle_index in range(0, len(cycles_per_cluster)):
            if cycles_per_cluster[0][cycle_length_index:] in cycles_per_cluster[cycle_index]:
                continue
            else:
                return count+k
        count+=1
    return count+k

def count_elements_in_dictionary(cycles):
    sum_of_quantity = [len(cycles[k_mer]) for k_mer in cycles.keys()]
    return sum(sum_of_quantity)

def remove_all_with_difference_less_than_23(cycles,limits=55):
    count=0
    new_dictionary_of_cycles = {}
    new_dictionary_of_multiplicities = {}
    for k_mer in cycles.keys():
        if len(cycles[k_mer]) == 1:
            new_dictionary_of_cycles[k_mer] = cycles[k_mer]
            continue
        value = []
        key = k_mer
        #print("K-mer: "+str(k_mer))
        right_side = shift_left_common_word_length(cycles[k_mer])
        left_side = shift_right_common_word_length(cycles[k_mer])
        if left_side>=47 or right_side>=47:
            count+=1
            continue
        common_elements = right_side+left_side
        for cycle in cycles[k_mer]:
            if len(cycle) - common_elements <= 21 or len(cycle) - common_elements > limits:
                count+=1
            else:
                value.append(cycle)
        if len(value) > 0:
            new_dictionary_of_cycles[key] = value
                
    sum_of_quantity = [len(new_dictionary_of_cycles[k_mer]) for k_mer in new_dictionary_of_cycles.keys()]
    return new_dictionary_of_cycles

def write_to_file(cycles, file_name):
    cycle_set = set()
    for k_mer in cycles.keys():
        for cycle in cycles[k_mer]:
            cycle_set.add(cycle)
    cycles_result = list(cycle_set)
    cycles_result.sort()
    with open(file_name, 'w') as f:
        for cycle in cycles_result:
            f.write(cycle+"\n") 
    
    
def count_singletons(cycles):
    count=0
    for k_mer in cycles.keys():
        if len(cycles[k_mer]) <= 3:
            count+=1
    return count 

import math
def char_comp(cycles):
    if len(cycles) == 1:
        return [0]
    lengths = [len(cycle) for cycle in cycles]
    min_length = min(lengths)
    characters = [cycle[math.ceil(min_length/2)] for cycle in cycles]
    first_character = characters[0]
    vector_of_scores = []
    for character in characters:
        vector_of_scores.append(0 if first_character == character else 1)
    return vector_of_scores

def calculate_mean_of_scores(cycles):
    if len(cycles) == 1:
        return -1
    scores = char_comp(cycles)
    return sum(scores)/len(scores)

def calculate_mean_of_scores_for_all_systems(cycles):
    mean_of_scores = {}
    for k_mer in cycles.keys():
        mean_of_scores[k_mer] = calculate_mean_of_scores(cycles[k_mer])    
    return mean_of_scores

def keep_only_those_with_score_greater_than_0(cycles):
    new_dictionary_of_cycles = {}
    mean_of_scores=calculate_mean_of_scores_for_all_systems(cycles)
    for k_mer in mean_of_scores.keys():
        if mean_of_scores[k_mer] > 0.25 or mean_of_scores[k_mer] == -1:
            new_dictionary_of_cycles[k_mer] = cycles[k_mer]
    return new_dictionary_of_cycles

cycles = read("/vol/d/development/git/mCAAT/proof_of_concept/data/"+sys.argv[1]+"/cycles_genome/results.txt")
sum_of_quantity = [len(cycles[k_mer]) for k_mer in cycles.keys()]
print("Elements before filter: "+str(sum(sum_of_quantity)))
#print first 2 elements in mults
count = 0
limits = 55
filtered_cycles = remove_all_with_difference_less_than_23(cycles,limits)
print("Elements after first round of filtering:",count_elements_in_dictionary(filtered_cycles))
f_c = keep_only_those_with_score_greater_than_0(filtered_cycles)
write_to_file(f_c, "/vol/d/development/git/mCAAT/proof_of_concept/data/"+sys.argv[1]+"/cycles_genome/results_filtered.txt")
print("Elements after second round of filtering:",count_elements_in_dictionary(f_c))
print("Elements removed:",count_elements_in_dictionary(cycles)-count_elements_in_dictionary(f_c))
