#include "crispr_array/filters.h"

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

    std::vector<uint64_t> Filters::FindRepeatNodePaths(vector<uint64_t> repeat_nodes, uint64_t start_node) {
        uint64_t start = 0;
        vector<uint64_t> all_the_neighbors;

        for (const auto& node : repeat_nodes) {
            uint64_t outgoings[4];
            int num_outgoings = this->sdbg.OutgoingEdges(node, outgoings);
            for (int i = 0; i < num_outgoings; i++) {
                all_the_neighbors.push_back(outgoings[i]);
            }
        }
        for (const auto& node : repeat_nodes)
            if (std::find(all_the_neighbors.begin(), all_the_neighbors.end(), node) == all_the_neighbors.end())
                start = node;

        int maxSize = 0;
        vector<vector<uint64_t>> cycles_per_group = this->cycles[start_node];
        std::vector<uint64_t> arr;
        auto it = std::find(arr.begin(), arr.end(), start);
        int position_to_rotate = std::distance(arr.begin(), it);

        for (int i = 0; i < cycles_per_group.size(); i++) {
            if (cycles_per_group[i].size() > maxSize) {
                maxSize = cycles_per_group[i].size();
                arr = cycles_per_group[i];
            }
        }

        arr.resize(repeat_nodes.size());
        return arr;
    }

    pair<vector<uint64_t>, vector<vector<uint64_t>>> Filters::_FindCRISPRArrayNodes(uint64_t start_node) {
        std::unordered_map<uint64_t, int> element_count;
        string start_node_label = std::to_string(start_node);
        if (cycles.find(start_node) == cycles.end()) {
            std::cerr << "Logging: " << start_node_label << " has been removed from consideration" << std::endl;
            return {{}, {}};
        }
        auto data = cycles[start_node];
        if (data.size() < 2) {
            return {{}, {}};
        }
        for (auto& vec : data) {
            vec.pop_back();
        }
        cycles[start_node] = data;

        int threshold = static_cast<int>(data.size());

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

        repeat_nodes = FindRepeatNodePaths(repeat_nodes, start_node);

        for (auto& vec : this->cycles[start_node]) {
            if (vec.size() - repeat_nodes.size() >= 23) {
                std::vector<uint64_t> spacers(vec.begin() + repeat_nodes.size(), vec.end());
                spacer_nodes.push_back(spacers);
            }
        }
        if (repeat_nodes.size() == 0 || spacer_nodes.size() < 3) {
            return {{}, {}};
        }
        return {repeat_nodes, spacer_nodes};
    }

    unordered_map<string, vector<string>> Filters::ListArrays(int& number_of_spacers) {
        unordered_map<string, vector<string>> CRISPRArrays;
        for (const auto& [start_node, _] : cycles) {
            auto CRISPRArrayNodes = _FindCRISPRArrayNodes(start_node);
            auto spacers_nodes = CRISPRArrayNodes.second;
            vector<uint64_t> repeat_nodes = CRISPRArrayNodes.first;
            if (!CRISPRArrayNodes.first.empty() && !CRISPRArrayNodes.second.empty()) {
                string repeat = _FetchNodeLabel(repeat_nodes[0]);
                for (size_t i = 1; i < repeat_nodes.size(); i++) {
                    std::string node_label = _FetchNodeLabel(repeat_nodes[i]);
                    char lastChar = node_label.back();
                    std::string lastCharStr(1, lastChar);
                    repeat += lastCharStr;
                }
                vector<vector<uint64_t>> cycles_nodes = this->cycles[start_node];
                vector<string> spacers_temp;
                vector<string> spacers;
                string all_cycles_togehter;
                for (const auto& cycle : cycles_nodes) {
                    std::string cycle_str = _FetchNodeLabel(cycle[0]);
                    for (size_t i = 1; i < cycle.size(); i++) {
                        uint64_t node = cycle[i];
                        std::string node_label = _FetchNodeLabel(node);
                        char lastChar = node_label.back();
                        std::string lastCharStr(1, lastChar);
                        cycle_str += lastCharStr;
                    }
                    all_cycles_togehter += cycle_str.substr(0, cycle_str.size() - 21);
                }
                size_t start = 0;
                size_t end;
                while ((end = all_cycles_togehter.find(repeat, start)) != std::string::npos) {
                    std::string part = all_cycles_togehter.substr(start, end - start);
                    if (!part.empty()) {
                        spacers_temp.push_back(part);
                    }
                    start = end + repeat.size();
                }
                if (start < all_cycles_togehter.size()) {
                    spacers_temp.push_back(all_cycles_togehter.substr(start));
                }
                for (const auto& spacer : spacers_temp) {
                    if (spacer.size() < 23 || spacer.size() > 50)
                        continue;
                    spacers.push_back(spacer.substr(0, spacer.size()));
                    number_of_spacers++;
                }
                if (spacers.size() < 2) {
                    number_of_spacers -= spacers.size();
                    continue;
                }
                CRISPRArrays[repeat] = spacers;
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
            file << "Repeat: " << repeat << endl;
            file << "Number of Spacers: " << spacers.size() << endl;
            file << "Spacers:" << endl;
            for (const auto& spacer : spacers) {
                file << spacer << endl;
            }
            file << "----------------------------------" << endl;
        }
        file.close();
        return number_of_spacers;
    }