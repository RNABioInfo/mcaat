#include "tmp_utils.h"

void trim_string(string& s) {
    s.erase(0, s.find_first_not_of(" \t\n\r"));
    s.erase(s.find_last_not_of(" \t\n\r") + 1);
}

pair<string, optional<string>> get_fastq_files_from_settings(
    const Settings& settings
) {
    const bool two_fastq_files_provided = settings.input_files.find(" ") != string::npos;
    if (two_fastq_files_provided) {
        const size_t space_pos = settings.input_files.find(" ");
        string input_file1 = settings.input_files.substr(0, space_pos);
        string input_file2 = settings.input_files.substr(space_pos + 1);

        trim_string(input_file1);
        trim_string(input_file2);

        return std::make_pair(input_file1, optional<string>(input_file2));
    } else {
        return std::make_pair(settings.input_files, std::nullopt);
    }
}

vector<vector<uint64_t>> cycles_map_to_cycles(
    const unordered_map<uint64_t, vector<vector<uint64_t>>>& cycles_map
) {
    vector<vector<uint64_t>> cycles;

    for (const auto& [_, inner_cycles] : cycles_map) {
        for (const auto& cycle : inner_cycles) {
            cycles.push_back(cycle);
        }
    }

    return cycles;
}

int get_cycle_count(
    const vector<vector<uint64_t>>& cycles
) {
    int count = 0;
    for (const auto& cycle : cycles) {
        count += cycle.size();
    }
    return count;
}

void print_sdbg_graph_to_dot_file_convert(const string& lib_file_path) {
    SDBG sdbg;
    sdbg.LoadFromFile(lib_file_path.c_str());


    std::cout << "--------------DOT-FILE------------------" << std::endl;
    std::cout << "digraph TinyGraph {\n";

    for (size_t i = 0; i < sdbg.size(); ++i) {
        std::cout << "    " << i << ";\n";
    }

    for (size_t node = 0; node < sdbg.size(); ++node) {
        int outdegree = sdbg.EdgeOutdegree(node);
        if (outdegree > 0) {
            uint64_t* outgoings = new uint64_t[outdegree];
            if (sdbg.OutgoingEdges(node, outgoings) != -1) {
                for (int i = 0; i < outdegree; ++i) {
                    const uint64_t neighbor = outgoings[i];
                    if (!sdbg.IsValidEdge(neighbor)) continue;

                    std::cout << "    " << node << " -> " << neighbor << ";\n";
                }
            }
            delete[] outgoings;
        }
    }

    std::cout << "}\n";
    std::cout << "--------------DOT-FILE-END--------------" << std::endl;
}


string fetch_node_label(SDBG& sdbg, const size_t& node) {
    string label;
    uint8_t seq[sdbg.k()];
    uint32_t t = sdbg.GetLabel(node, seq);
    for (int i = sdbg.k() - 1; i >= 0; --i) label.append(1, "ACGT"[seq[i] - 1]);
    reverse(label.begin(), label.end());
    return label;
}

pair<vector<uint64_t>, string> Find_CRISPR_repeat_nodes_and_sequence(
    SDBG& sdbg,
    const vector<vector<uint64_t>>& ordered_cycles
) {
    const int threshold = static_cast<int>(0.9 * static_cast<float>(ordered_cycles.size()) + 0.5);
    
    unordered_map<uint64_t, int> element_count;
    for (const auto& cycle : ordered_cycles) {
        for (const auto& element : cycle) {
            element_count[element]++;
        }
    }

    // Get all repeat_nodes
    vector<uint64_t> repeat_nodes;
    unordered_set<uint64_t> repeat_nodes_set;
    for (const auto& cycle : ordered_cycles) {
        bool passed_spacers = false;
        bool going_through_repeat_region = false;
        int loop_counter = 0;

        // The spacers might wrap around, therefore go through cycle twice:
        // Starting with S -> SSSSSRRRRRSSSSS SSSSSRRRRRSSSSS
        // Starting with R -> RRRRRSSSSSRRRRR RRRRRSSSSSRRRRR
        while (loop_counter < 2) {
            for (const auto& node : cycle) {
                bool is_spacer_node = element_count[node] < threshold;
                if (is_spacer_node) {
                    passed_spacers = true;
                    if (going_through_repeat_region) {
                        break;
                    }
                } else if (passed_spacers) {
                    going_through_repeat_region = true;
                    bool already_in_set = repeat_nodes_set.insert(node).second;
                    if (!already_in_set) {
                        repeat_nodes.push_back(node);
                    }
                }
            }

            loop_counter++;
        }
        for (const auto& [element, count] : element_count) {
            if (static_cast<float>(count) >= threshold) {
                repeat_nodes.push_back(element);
            }
        }
    }

    // Choose repeat node sequence that occur the most
    int repeat_nodes_count = 0;
    int best_idx = 0;
    for (int i = 0; i < ordered_cycles.size(); ++i) {
        const auto cycle = ordered_cycles.at(i);

        int current_repeat_nodes_count = 0;
        for (const auto& node : cycle) {
            bool is_repeat = std::find(repeat_nodes.begin(), repeat_nodes.end(), node) != repeat_nodes.end();
            if (is_repeat) {
                current_repeat_nodes_count += element_count[node];
            }
        }

        if (current_repeat_nodes_count > repeat_nodes_count) {
            repeat_nodes_count = current_repeat_nodes_count;
            best_idx = i;
        }
    }

    // Fetch best sequence of repeat nodes
    string repeat = "";
    bool passed_spacers = false;
    bool going_through_repeat_region = false;
    int loop_counter = 0;

    // The spacers might wrap around, therefore go through cycle twice:
    // Starting with S -> SSSSSRRRRRSSSSS SSSSSRRRRRSSSSS
    // Starting with R -> RRRRRSSSSSRRRRR RRRRRSSSSSRRRRR
    while (loop_counter < 2) {
        for (const auto& node : ordered_cycles.at(best_idx)) {
            bool is_repeat_node = std::find(repeat_nodes.begin(), repeat_nodes.end(), node) != repeat_nodes.end();
            if (!is_repeat_node) {
                passed_spacers = true;
                if (going_through_repeat_region) {
                    // Idk why, but one letter is missing
                    char new_char = fetch_node_label(sdbg, node).back();
                    string new_char_str(1, new_char);
                    repeat += new_char_str;
                    break;
                }
            } else if (passed_spacers) {
                going_through_repeat_region = true;
                if (repeat == "") {
                    repeat = fetch_node_label(sdbg, node);
                } else {
                    char new_char = fetch_node_label(sdbg, node).back();
                    string new_char_str(1, new_char);
                    repeat += new_char_str;
                }
            }
        }

        loop_counter++;
    }

    return {repeat_nodes, repeat};
}

tuple<string, vector<string>, string> get_systems(
    SDBG& sdbg,
    vector<vector<uint64_t>>& ordered_cycles
) {
    int smallest_cycle_size = ordered_cycles.at(0).size();
    for (const auto& cycle : ordered_cycles) {
        if (cycle.size() < smallest_cycle_size) {
            smallest_cycle_size = cycle.size();
        }
    }

    /* Repeat node extend to the right */
    int extension_to_right = 0;
    for (int i = 0; i < smallest_cycle_size - 1; ++i) {
        unordered_set<char> next_bp;
        for (const auto& cycle : ordered_cycles) {
            const uint64_t current_node = cycle.at(i);
            const string label = fetch_node_label(sdbg, current_node);
            next_bp.insert(label.at(0));
        }

        const bool bp_branched_off = next_bp.size() > 1;
        if (bp_branched_off) {
            unordered_set<char> next_next_bp;
            for (const auto& cycle: ordered_cycles) {
                const uint64_t current_node = cycle.at(i + 1);
                const string label = fetch_node_label(sdbg, current_node);
                next_next_bp.insert(label.at(0));
            }
            
            const bool is_point_mutation = next_next_bp.size() == 1;
            if (!is_point_mutation) {
                extension_to_right = i;
                break;
            }
        }
    }

    /* Repeat node extend to the left */
    int extension_to_left = 0;
    for (int i = 0; i < smallest_cycle_size - 1; ++i) {
        unordered_set<char> prev_bp;
        for (const auto& cycle: ordered_cycles) {
            const uint64_t current_node = cycle.at(cycle.size() - i - 1);
            const string label = fetch_node_label(sdbg, current_node);
            prev_bp.insert(label.back());
        }

        const bool bp_branched_off = prev_bp.size() > 1;
        if (bp_branched_off) {
            unordered_set<char> prev_prev_bp;
            for (const auto& cycle: ordered_cycles) {
                const uint64_t inner_current_node = cycle.at(cycle.size() - i - 2);
                const string inner_label = fetch_node_label(sdbg, inner_current_node);
                prev_prev_bp.insert(inner_label.back());
            }
            
            const bool is_point_mutation = prev_prev_bp.size() == 1;
            if (!is_point_mutation) {
                extension_to_left = i;
                break;
            }
        }
    }

    const int repeat_length = extension_to_left + extension_to_right - sdbg.k();

    /* Turn cycles into repeats and spacers */
    vector<string> spacers;
    vector<string> repeats;
    for (const auto& cycle : ordered_cycles) {
        string spacer;
        string repeat;

        const int offset_repeat = cycle.size() - extension_to_left;

        for (int i = 0; i < cycle.size(); ++i) {
            const uint64_t node = cycle.at((offset_repeat + i) % cycle.size());
            const string label = fetch_node_label(sdbg, node);
            
            if (i < repeat_length) {
                repeat += label.back();
            } else {
                spacer += label.back();
            }
        }

        spacers.push_back(spacer);
        repeats.push_back(repeat);
    }

    /* Get consensus/most-common repeat */
    unordered_map<string, int> repeat_count;
    for (const string& repeat : repeats) {
        repeat_count[repeat]++;
    }

    string consensus_repeat = "";
    int max_count = 0;
    for (const auto& [repeat, count] : repeat_count) {
        if (count > max_count) {
            max_count = count;
            consensus_repeat = repeat;
        }
    }

    /* Reconstruct sequences */
    string full_sequence = "";
    for (int i = 0; i < spacers.size(); ++i) {
        const string repeat = repeats.at(i);
        const string spacer = spacers.at(i);

        if (repeat == consensus_repeat) {
            full_sequence += repeat;
            full_sequence += spacer;
        }
    }

    /* crispr sequence end with a repat (usually mutated), tmp taking consensus repeat */
    full_sequence += consensus_repeat;

    return std::make_tuple(consensus_repeat, spacers, full_sequence);
}
