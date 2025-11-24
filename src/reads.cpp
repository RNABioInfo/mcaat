#include "reads.h"

vector<string> extract_sequences_from_fastq_file(const string& fastq_file_path) {
    vector<string> sequences;

    try {
        klibpp::SeqStreamIn iss(fastq_file_path.c_str());
        auto file_records = iss.read();
        for (const auto& record : file_records) {
            sequences.push_back(record.seq);
        }
    } catch (const std::exception& e) {
        std::cerr << "Error reading file \"" << fastq_file_path;
        std::cerr << "\" sequences because: " << e.what() << std::endl;
    }

    return sequences;
}

void reverse_pair_ends_sequence(string& sequence) {
    std::reverse(sequence.begin(), sequence.end());

    for (char& base : sequence) {
        switch (base) {
            case 'A': base = 'T'; break;
            case 'T': base = 'A'; break;
            case 'C': base = 'G'; break;
            case 'G': base = 'C'; break;
        }
    }
}

uint64_t k_mer_to_node_id(
    const SDBG& sdbg,
    const string k_mer
) {
    if (k_mer.size() != sdbg.k()) {
        return 0;
    }

    uint8_t seq[sdbg.k()];
    for (int i = 0; i < sdbg.k(); ++i) {
        const char c = k_mer.at(i);
        if (c == 'A') {
            seq[i] = 1;
        } else if (c == 'C') {
            seq[i] = 2;
        } else if (c == 'G') {
            seq[i] = 3;
        } else {
            seq[i] = 4;
        }
    }
    return sdbg.IndexBinarySearch(seq);
}

vector<uint64_t> get_read_from_sequence(
    const SDBG& sdbg,
    const unordered_set<uint64_t>& nodes_of_cycles,
    const string& sequence
) {
    const uint32_t K = sdbg.k();

    if (sequence.size() <= 2 * K) {
        return {};
    }

    const string start_k_mer = sequence.substr(0, K);
    const string end_k_mer = sequence.substr(sequence.size() - K, K);

    const uint64_t start_node_id = k_mer_to_node_id(sdbg, start_k_mer);
    const uint64_t end_node_id = k_mer_to_node_id(sdbg, end_k_mer);

    if (nodes_of_cycles.find(start_node_id) == nodes_of_cycles.end()
    && nodes_of_cycles.find(end_node_id) == nodes_of_cycles.end()) {
        return {};
    }

    vector<uint64_t> read = {start_node_id};
    for (int i = 1; i < sequence.size() - K; ++i) {
        const string k_mer = sequence.substr(i, K);
        const uint64_t node_id = k_mer_to_node_id(sdbg, k_mer);

        read.push_back(node_id);
    }
    read.push_back(end_node_id);

    return read;
}

vector<vector<uint64_t>> get_reads(
    const SDBG& sdbg,
    const string& fastq_file_1,
    const optional<string>& fastq_file_2,
    const vector<vector<uint64_t>> cycles
) {

    unordered_set<uint64_t> nodes_of_cycles;
    for (const auto& cycle : cycles) {
        for (const auto& node : cycle) {
            nodes_of_cycles.insert(node);
        }
    }

    vector<vector<uint64_t>> reads;
    {
        vector<string> sequences = extract_sequences_from_fastq_file(fastq_file_1);
        for (const auto& sequence : sequences) {
            vector<uint64_t> read = get_read_from_sequence(sdbg, nodes_of_cycles, sequence);
            if (!read.empty()) {
                reads.push_back(read);
            }
        }
    }

    if (fastq_file_2.has_value()) {
        vector<string> sequences = extract_sequences_from_fastq_file(fastq_file_2.value());
        for (const auto& sequence : sequences) {
            string rev_sequence = sequence;
            reverse_pair_ends_sequence(rev_sequence);

            vector<uint64_t> read = get_read_from_sequence(sdbg, nodes_of_cycles, rev_sequence);
            if (!read.empty()) {
                reads.push_back(read);
            }
        }
    }

    return reads;
}
