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
        
    

    // ROTATE THE CYCLE TO START AT THE PROPER START
    vector<uint64_t> arr = this->cycles[start_node][0];
    auto it = std::find(arr.begin(), arr.end(), start);
    int position_to_rotate = std::distance(arr.begin(), it);
    rotateLeft(arr, position_to_rotate);
    cout<<"Position to rotate: "<<position_to_rotate<<endl;
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

    const auto& data = cycles[start_node];
    if (data.size() < 2) {
        return {{}, {}};
    }
    cout<<"StartNode: "<< start_node << endl;
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
        if (count >= threshold*0.95) {
            repeat_nodes.push_back(element);
        }
    }

    if (repeat_nodes.size() >= 27) {
        return {{}, {}};
    }
    
    std::vector<std::vector<uint64_t>> spacer_nodes;
    for (const auto& vec : data) {
        if (vec.size() - repeat_nodes.size() >= 23) {
            std::vector<uint64_t> spacers;
            for (const auto& element : vec) {
                if (std::find(repeat_nodes.begin(), repeat_nodes.end(), element) == repeat_nodes.end()) {
                    spacers.push_back(element);
                }
                else {
                    break;
                }
            }
            spacer_nodes.push_back(spacers);
        }
    }

    repeat_nodes = FindRepeatNodePaths(repeat_nodes,start_node);
    if(repeat_nodes.size() == 0 || spacer_nodes.size() < 3) {
        return {{}, {}};
    }
    return {repeat_nodes, spacer_nodes};

}

unordered_map<string, vector<string>> Filters::ListArrays() {
    unordered_map<string, vector<string>> CRISPRArrays;
    for (const auto& [start_node, _] : cycles) {
        auto CRISPRArrayNodes = _FindCRISPRArrayNodes(start_node);
        vector<uint64_t> repeat_nodes = CRISPRArrayNodes.first;
        if (!CRISPRArrayNodes.first.empty() && !CRISPRArrayNodes.second.empty()) {
            string repeat = _FetchNodeLabel(repeat_nodes[0]);
            for (size_t i = 1; i < repeat_nodes.size(); ++i) {
                std::string node_label = _FetchNodeLabel(repeat_nodes[i]);
                    // Method 1: Using back() method
                char lastChar = node_label.back();  // Get the last character
                 std::string lastCharStr(1, lastChar);  // Convert char to string
                
                repeat += lastCharStr;
            }
            vector<string> spacers;
            for (const auto& spacer : CRISPRArrayNodes.second) {
                /*
                string spacer_str;
                for (const auto& node : spacer) {
                    spacer_str += _FetchNodeLabel(node);
                }
                spacers.push_back(spacer_str);
                */
            }

            CRISPRArrays[repeat] = spacers;
        }
    }
    return CRISPRArrays;
}

void Filters::WriteToFile(const string& filename) {
    
    ofstream file(filename);
    if (!file.is_open()) {
        throw runtime_error("Failed to open file: " + filename);
    }

    auto CRISPRArrays = ListArrays();
    for (const auto& [repeat, spacers] : CRISPRArrays) {
        file << repeat << "\n";
        for (const auto& spacer : spacers) {
            file << spacer << "\n";
        }
    }
}
