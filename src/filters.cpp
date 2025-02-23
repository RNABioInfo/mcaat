#include "filters.h"

char Filters::_FetchLastCharacter(size_t node) {
    vector<uint8_t> seq(sdbg.k());
    sdbg.GetLabel(node, seq.data());
    return "ACGT"[seq[sdbg.k() - 1]];
}
string Filters::_FetchNodeLabel(size_t node) {
    std::string label;            
    uint8_t seq[sdbg.k()];
    uint32_t t = sdbg.GetLabel(node, seq);
    for (int i = sdbg.k() - 1; i >= 0; --i) label.append(1, "ACGT"[seq[i] - 1]);
    reverse(label.begin(), label.end());
    return label;
}
Filters::Filters(SDBG& sdbg, std::unordered_map<uint64_t, std::vector<std::vector<uint64_t>>>& cycles)
    : sdbg(sdbg), cycles(cycles) {}


void Filters::rotateLeft(std::vector<uint64_t>& arr, int k) {
    int n = arr.size();
    if (n == 0) return; // Handle empty array case
    
    k = k % n; // In case k is greater than the size of the array
    if (k == 0) return; // No rotation needed if k is 0 or a multiple of n
    

    // Step 1: Reverse the first k elements
    std::reverse(arr.begin(), arr.begin() + k);
    
    // Step 2: Reverse the remaining n-k elements
    std::reverse(arr.begin() + k, arr.end());
    
    // Step 3: Reverse the entire array
    std::reverse(arr.begin(), arr.end());
}
std::vector<uint64_t> Filters::FindRepeatNodePaths(vector<uint64_t> repeat_nodes,uint64_t start_node) {
    uint64_t start;
    uint64_t end;
    vector<uint64_t> all_the_neighbors;
    
    //LIST
    for(const auto& node : repeat_nodes) {
        uint64_t outgoings[4]; 
        int num_outgoings = this->sdbg.OutgoingEdges(node, outgoings);
        for(int i = 0; i < num_outgoings; i++) {
            all_the_neighbors.push_back(outgoings[i]);
        }
    }
    for(const auto& node : repeat_nodes) 
        if(std::find(all_the_neighbors.begin(), all_the_neighbors.end(), node) == all_the_neighbors.end()) 
            start = node;
        
    

    //previous suggestion is not correct
        int maxSize = 0;
    vector<vector<uint64_t>> cycles_per_group = this->cycles[start_node];
    std::vector<uint64_t> arr;
    auto it = std::find(arr.begin(), arr.end(), start);
    int position_to_rotate = std::distance(arr.begin(), it);
    // Loop through the list of vectors and find the one with the max number of elements
    for (int i = 0; i < cycles_per_group.size(); i++) {
        rotateLeft(this->cycles[start_node][i], position_to_rotate);
        if (this->cycles[start_node][i].size() > maxSize) {
            maxSize = this->cycles[start_node][i].size();
            arr = this->cycles[start_node][i];  // Update the vector with the maximum size
        }
    }
    

    this->cycles[start_node] = cycles_per_group;
    rotateLeft(arr, position_to_rotate);

    arr.resize(repeat_nodes.size());
    // END OF ROTATION
    return arr;
}

pair<vector<uint64_t>, vector<vector<uint64_t>>> Filters::_FindCRISPRArrayNodes(uint64_t start_node) {
    std::unordered_map<uint64_t, int> element_count;
    
    if (cycles.find(start_node) == cycles.end()) {
        std::cerr << "Error: start_node not found in cycles." << std::endl;
        return {{}, {}};
    }
    auto data = cycles[start_node];
    if (data.size() < 2) {
        return {{}, {}};
    }
    //delete every last element in vector of cycles[start_node] so that it is modified
    for (auto &vec : data) {
        vec.pop_back();
    }
    cycles[start_node] = data;
    
    int threshold = static_cast<int>(data.size()); // 85% threshold

    for (const auto& vec : data) {
        std::unordered_set<uint64_t> uniqueElements(vec.begin(), vec.end());
        for (const auto& element : uniqueElements) {
            element_count[element]++;
        }
    }

    if (data.empty() || data[0].empty()) {
        std::cerr << "Error: data or data[0] is empty." << std::endl;
        return {{}, {}};
    }
    std::vector<uint64_t> repeat_nodes;
    for (const auto& [element, count] : element_count) {
        if (count >= threshold) {
            repeat_nodes.push_back(element);
        }
    }

    if (repeat_nodes.size() >= 27) {
        return {{}, {}};
    }
    
    std::vector<std::vector<uint64_t>> spacer_nodes;
    

    repeat_nodes = FindRepeatNodePaths(repeat_nodes,start_node);
    
    for (auto& vec : this->cycles[start_node]) {
        // IMPORTANT: DO NOT ADD UPPER BOUNDARY

        if (vec.size() - repeat_nodes.size() >= 23) {
            std::vector<uint64_t> spacers(vec.begin() + repeat_nodes.size(), vec.end());
            
            spacer_nodes.push_back(spacers);
        // IMPORTANT: DO NOT ADD UPPER BOUNDARY
        }
    }
    if(repeat_nodes.size() == 0 || spacer_nodes.size() < 3) {
        return {{}, {}};
    }
    return {repeat_nodes, spacer_nodes};

}
bool Filters::isDistanceValid(const std::vector<int>& positions, int minDistance, int maxDistance) {
    for (size_t i = 1; i < positions.size(); ++i) {
        int distance = positions[i] - positions[i - 1];
        if (distance < minDistance || distance > maxDistance) {
            return false;
        }
    }
    return true;
}

// Function to find the most frequent sequence with distance constraints
std::string Filters::findMostFrequentSequence(const std::string& input, int minLength, int maxLength, int minDistance, int maxDistance) {
    std::unordered_map<std::string, std::vector<int>> substringPositions;
    int inputLength = input.size();

    // Generate substrings of lengths between minLength and maxLength
    for (int length = minLength; length <= maxLength; ++length) {
        for (int i = 0; i <= inputLength - length; ++i) {
            std::string substring = input.substr(i, length);
            substringPositions[substring].push_back(i);
        }
    }

    // Find the most frequent substring that satisfies the distance condition
    std::string mostFrequentSequence;
    int maxFrequency = 0;

    for (const auto& entry : substringPositions) {
        const std::string& substring = entry.first;
        const std::vector<int>& positions = entry.second;

        if (positions.size() > 1 && isDistanceValid(positions, minDistance, maxDistance)) {
            int frequency = positions.size();
            if (frequency > maxFrequency) {
                mostFrequentSequence = substring;
                maxFrequency = frequency;
            }
        }
    }

    return mostFrequentSequence;
}

unordered_map<string, string> Filters::ListArrays(int& number_of_spacers) {
    unordered_map<string, string> CRISPRArrays;
    int counter = 0;
    for (const auto& [start_node, _] : cycles) {
        auto CRISPRArrayNodes = _FindCRISPRArrayNodes(start_node);
        auto spacers_nodes = CRISPRArrayNodes.second;
        vector<uint64_t> repeat_nodes = CRISPRArrayNodes.first;
        if (!CRISPRArrayNodes.first.empty() && !CRISPRArrayNodes.second.empty()) {
            string repeat = _FetchNodeLabel(repeat_nodes[0]);
            
            for (size_t i = 1; i < repeat_nodes.size(); i++) {
                std::string node_label = _FetchNodeLabel(repeat_nodes[i]);
                    // Method 1: Using back() method
                char lastChar = node_label.back();  // Get the last character
                 std::string lastCharStr(1, lastChar);  // Convert char to string
              
                repeat += lastCharStr;
            }
            vector<vector<uint64_t>> cycles_nodes =this->cycles[start_node];
            //print the first spacer, first cycle and repeat nodes
            
            vector<string> spacers_temp;
            vector<string> spacers;
            string all_cycles_togehter;
            for (const auto& cycle : cycles_nodes) {
                std::string cycle_str =_FetchNodeLabel(cycle[0]);
                for (size_t i = 1; i < cycle.size(); i++) {
                    uint64_t node = cycle[i];
                    std::string node_label = _FetchNodeLabel(node);
                    // Method 1: Using back() method
                    char lastChar = node_label.back();  // Get the last character
                    std::string lastCharStr(1, lastChar);  // Convert char to string
                    cycle_str += lastCharStr;
                    
                }
               
                all_cycles_togehter += cycle_str.substr(0, cycle_str.size()-21);
            }
            
            size_t start = 0;
            size_t end;

            // Iterate through the string and find substrings
            while ((end = all_cycles_togehter.find(repeat, start)) != std::string::npos) {
                std::string part = all_cycles_togehter.substr(start, end - start);
                if (!part.empty()) {
                    spacers_temp.push_back(part);
                }
                start = end + repeat.size();
            }

            // Add the remaining part after the last delimiter
            if (start < all_cycles_togehter.size()) {
                spacers_temp.push_back(all_cycles_togehter.substr(start));
            }
            for(const auto& spacer : spacers_temp) {
                //cout<<spacer.size()<<endl;
                if(spacer.size() < 23 || spacer.size() > 50)
                    continue;
                
                spacers.push_back(spacer.substr(0, spacer.size()));
                number_of_spacers++;
                
            }
            if(spacers.size() < 3){
                number_of_spacers -= spacers.size();
                continue;
            }

            CRISPRArrays[repeat] = all_cycles_togehter;
        }
    }
    
    return CRISPRArrays;
}

int Filters::WriteToFile(const string& filename) {
    
    ofstream file(filename);
    if (!file.is_open()) {
        throw runtime_error("Failed to open file: " + filename);
    }
    int number_of_spacers = 0;
    
    auto CRISPRArrays = ListArrays(number_of_spacers);
    
    for (const auto& [repeat, spacers] : CRISPRArrays) {
        
        file << spacers << "\n";

        
    }
    file.close();
    
    return number_of_spacers;

}
