import sys
import needle_in_haystack

def read_file(folder, file_name):
    output_list = []
    with open(folder + file_name, "r") as f:
        output_list = [[int(n) for n in line.split()] for line in f]
    return output_list

def write_dict_to_file(dictionary, file_name):
    with open(file_name, "w") as f:
        for key in dictionary.keys():
            f.write(str(key) + " " + "\n")
            for element in dictionary[key]:
                f.write(" ".join(str(e) for e in element) + "\n")
    


def track_first_ids(input_nodes_list) -> list:
    tracker = [0]
    first_element = input_nodes_list[0][0]
    print(len(input_nodes_list))
    for index_of_element in range(len(input_nodes_list)):
        if first_element!=input_nodes_list[index_of_element][0]:
            tracker.append(index_of_element)
        first_element = input_nodes_list[index_of_element][0]
    return tracker

def group_by_first_element(input_nodes_list,input_multiplicity_list,tracker) -> tuple:
    # Input: list of lists of nodes, list of lists of multiplicities, list of indexes of first elements
    # Output: dictionary of lists of nodes grouped by first element, dictionary of lists of multiplicities grouped by first id
    group_cycles = {}
    group_multiplicities = {}
    for i in range(len(tracker)):
        f_i = tracker[i]
        l_i = len(input_nodes_list)
        if i<len(tracker)-1:
            l_i = tracker[i+1]
        group_cycles[input_nodes_list[f_i][0]] = input_nodes_list[f_i:l_i]
        group_multiplicities[input_nodes_list[f_i][0]] = input_multiplicity_list[f_i:l_i]
    return group_cycles,group_multiplicities

def rearrange_cyles(id_paths,multiplicities):
    indicies = []
    for i in range(len(multiplicities)):
        indicies.append(multiplicities[i].index(max(multiplicities[i])))
    new_id_paths = []
    new_multiplicities=[]
    for i in range(len(id_paths)):
        start = indicies[i] 
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


def rearrange(element,mult,id) -> tuple:
    # Input: list of nodes, list of multiplicities, id of first element
    # Output: list of nodes, list of multiplicities
    # get index of id in element
    index_of_id = element.index(id)
    edited_element = element[index_of_id:].extend(element[1:index_of_id])
    edited_multiplicity_element = mult[index_of_id:].extend(mult[1:index_of_id])
    return edited_element, edited_multiplicity_element

def adjust_cycles_to_frequency(group_cycles,group_multiplicities):
    first_ids, all_elements = group_cycles.keys(), group_cycles.values()
    all_mults = group_multiplicities.values()
    for id in first_ids:
        for element, mult in zip(all_elements, all_mults):
            if id in element[1:-1]:
                group_cycles[element[0]].remove(element)
                group_multiplicities[element[0]].remove(mult)
                edited_element, edited_multiplicity_element = rearrange(element, mult, id)
                group_cycles[id].append(edited_element)
                group_multiplicities[id].append(edited_multiplicity_element)
    return group_cycles, group_multiplicities


from collections import Counter
from collections import defaultdict


def count_of_every_element(set_of_elemnts,list_of_lists):
    count_of_element = {}
    for element in set_of_elemnts:
        count_of_element[element] = 0
        for list_of_elements in list_of_lists:
            if element in list_of_elements:
                count_of_element[element]+=1
    return count_of_element

def find_common_elements(arrays):
    # Convert each array to a set, excluding the first and last elements
    
    all_elements_in_first_array = set(arrays[0][1:-1])
    result_indices = []
    for element in all_elements_in_first_array:
       # Check if all arrays contain the element
        if all(element in array for array in arrays[1:]):
            result_indices.append(arrays[0].index(element))
    return len(result_indices)+2

def group_by_start_end(group_cycles):   
    group_common_elements = {}   
    for id in group_cycles.keys():
        if len(group_cycles[id])==1:
            group_common_elements[id] = 2
        group_common_elements[id] = find_common_elements(group_cycles[id]) 
    return group_common_elements

def clean_up(group_common_elements,group_cycles,group_multiplicities):
    for key in list(group_cycles.keys()):
        for cycle in group_cycles[key]:
            length_of_cycle = len(cycle) - group_common_elements[key]
            print(length_of_cycle, len(cycle), group_common_elements[key])
            if length_of_cycle<=22 or length_of_cycle>=75:
                index_at = group_cycles[key].index(cycle)
                group_cycles[key].remove(cycle)
                group_multiplicities[key].remove(group_multiplicities[key][index_at])
    return group_cycles, group_multiplicities

def separate_list_by_first_element(input_nodes_list,input_multiplicity_list,folder) -> None:
    # Input: list of lists of nodes, list of lists of multiplicities
    # Output: void, writes two dictionaries to file
    tracker = track_first_ids(input_nodes_list)
    to_keep = []
    #input_nodes_list, input_multiplicity_list = rearrange_cyles(input_nodes_list,input_multiplicity_list)
    group_cycles, group_multiplicities =  group_by_first_element(input_nodes_list,input_multiplicity_list,tracker)
    group_cycles, group_multiplicities = adjust_cycles_to_frequency(group_cycles,group_multiplicities)
    counter=0
    for key in group_cycles.keys():
        counter+=len(group_cycles[key])
    
    print(counter)
    write_dict_to_file(group_multiplicities,"../data/"+folder+"/cycles_genome/group_multiplicities.txt")
    write_dict_to_file(group_cycles,"../data/"+folder+"/cycles_genome/group_cycles.txt")



#id_paths = read_file("../data/"+sys.argv[1]+"/cycles_genome/","id_paths.txt")
#id_multiplicities = read_file("../data/"+sys.argv[1]+"/cycles_genome/","multiplicity_distribution.txt")
#separate_list_by_first_element(id_paths,id_multiplicities)
#print(common)