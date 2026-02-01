#include "main_run_and_debug.h"

vector<vector<uint64_t>> run_and_debug_finding_of_relevant_reads(
    const vector<vector<uint64_t>> &cycles,
    const Settings& settings,
    const SDBG &sdbg
) {
    chrono::_V2::system_clock::time_point start_time = chrono::high_resolution_clock::now();
    
    auto fastq_files = get_fastq_files_from_settings(settings);
    auto reads = get_reads(
        sdbg,
        fastq_files.first,
        fastq_files.second,
        cycles
    );
    cout << "    ▸ Found " << reads.size() << " reads" << endl;
    if (reads.size() == 0) {
        cout << "══════════════════════════════════════════════" << endl;
        return reads;
    }

    chrono::_V2::system_clock::time_point end_time = chrono::high_resolution_clock::now();
    cout << "\n⏳ Time elapsed: ";
    cout << std::fixed << std::setprecision(2);
    cout << chrono::duration<double>(end_time - start_time).count();
    cout << " seconds" << endl;

    return reads;
}

vector<tuple<string, string, vector<string>, float, float>> run_and_debug_spacer_ordering(
    const Settings& settings,
    SDBG& sdbg,
    const vector<vector<uint64_t>>& cycles
) {
    chrono::_V2::system_clock::time_point start_time = chrono::high_resolution_clock::now();
    
    cout << "  ▸ Splitting into subproblems" << endl;
    const size_t reads_size = get_reads_size_of_fastq_files(settings);
    auto subgraphs = get_crispr_regions_extended_by_k(sdbg, reads_size, cycles);

    cout << "  ▸ Extracting relevant reads from the whole file" << endl;
    const auto relevant_nodes = combine_nodes_of_subgraphs(subgraphs);
    
    auto fastq_files = get_fastq_files_from_settings(settings);
    vector<vector<uint64_t>> reads;

    {
        typedef decltype(open_fastq_file_stream("")) StreamType;
        std::unique_ptr<StreamType> file_stream_ptr = std::make_unique<StreamType>(fastq_files.first.c_str());
        bool both_fastq_files_streamed = false;
        long number_of_reads = 1; // 1 to pass first if condition

        do {
            if (number_of_reads == 0 && fastq_files.second.has_value()) {
                std::cout << "switching fastq_file" << std::endl;
                both_fastq_files_streamed = true;
                file_stream_ptr = std::make_unique<StreamType>(fastq_files.second.value().c_str());
            }
            number_of_reads = 0;

            try {
                auto file_records = file_stream_ptr->read(settings.chunk_size);
                for (const auto& record : file_records) {
                    string sequence = record.seq;
                    if (both_fastq_files_streamed) {
                        reverse_pair_ends_sequence(sequence);
                    }

                    vector<uint64_t> read = get_read_from_sequence(sdbg, relevant_nodes, sequence);
                    if (!read.empty()) {
                        number_of_reads++;
                        reads.push_back(read);
                    }
                }
            } catch (const std::exception& e) {
                std::cerr << "Error reading ";
                if (both_fastq_files_streamed) {
                    std::cerr << " second";
                } else {
                    std::cerr << " first";
                }
                std::cerr << " fastq file sequences because: " << e.what() << std::endl;
            }
        } while (number_of_reads > 0 || !both_fastq_files_streamed);
    }

    cout << "  🔄 Filtering subproblems:" << endl;
    vector<Graph> remaining_subgraphs;
    vector<vector<vector<uint64_t>>> remaining_reads;
    vector<vector<vector<size_t>>> remaining_cycles;
    for (size_t idx = 0; idx < subgraphs.size(); ++idx) {
        const auto& subgraph = subgraphs[idx];
        const auto relevant_reads = get_relevant_reads(subgraph.nodes, reads);
        auto relevant_cycles = get_relevant_cycles(subgraph, cycles);

        get_minimum_cycles_for_full_coverage(relevant_cycles);
        
        // The assembly of megahit always assembles the graph in its reverse complement.
        // We discard the reverse complement by assuming that it won't have any relevant reads
        if (relevant_reads.size() == 0 || relevant_cycles.size() < 3) {
            continue;
        }

        remaining_subgraphs.push_back(subgraph);
        remaining_reads.push_back(relevant_reads);
        remaining_cycles.push_back(relevant_cycles);
    }
    cout << "  ✅ Filtered out " << subgraphs.size() - remaining_subgraphs.size();
    cout << "/" << subgraphs.size() << " subproblems" << endl;

    
    cout << "  🔄 Solving " << remaining_subgraphs.size();
    cout << " subproblems..." << endl;
    vector<tuple<string, string, vector<string>, float, float>> found_systems;
    for (size_t idx = 0; idx < remaining_subgraphs.size(); ++idx) {
        const auto& subgraph = remaining_subgraphs[idx];
        const auto& relevant_reads = remaining_reads[idx];
        const auto& relevant_cycles = remaining_cycles[idx];
        
        cout << "    Subproblem " << idx + 1 << "/";
        cout << remaining_subgraphs.size() << ":" << endl;
        
        cout << "      🛈 Graph with " << subgraph.nodes.size();
        cout << " nodes and " << subgraph.edge_count() << " edges" << endl;

        cout << "      🛈 Reads with " << relevant_reads.size() << "/";
        cout << reads.size() << " used" << endl;

        cout << "      🛈 Cycles with " << relevant_cycles.size() << "/";
        cout << get_cycle_count(cycles) << " used" << endl;

        float confidence_cycle_resolution = 1.0;
        float confidence_topological_sort = 1.0;
        auto cycle_order = order_cycles(
            subgraph,
            relevant_reads,
            relevant_cycles,
            confidence_cycle_resolution,
            confidence_topological_sort
        );

        cout << "      ▸ The order is ";
        for (auto node : cycle_order) {
            cout << node << " ";
        }
        cout << endl;

        cout << "      ▸ Cycles were resolved with a confidence of ";
        cout << std::fixed << std::setprecision(2);
        cout << (confidence_cycle_resolution * 100) << "%" << endl;
        cout << "      ▸ Topological sort has a confidence of ";
        cout << (confidence_topological_sort * 100) << "%" << endl;

        cout << "      ▸ Turning the cycle order into a node order" << endl;
        auto ordered_cycles = get_ordered_cycles(cycle_order, relevant_cycles);
        
        if (ordered_cycles.size() < 2) {
            cout << "      ▸ Node order is to short and is not processed further" << endl;
            continue;
        }

        cout << "      ▸ Starting the filter process:" << endl;
        auto [repeat, spacers, full_sequence] = get_systems(sdbg, ordered_cycles);
        
        cout << "        ▸ Number of spacers: " << spacers.size() << endl;

        found_systems.push_back(std::make_tuple(
            full_sequence,
            repeat,
            spacers,
            confidence_cycle_resolution,
            confidence_topological_sort
        ));
    }
    cout << "  ✅ Completed each subproblem" << endl;

    chrono::_V2::system_clock::time_point end_time = chrono::high_resolution_clock::now();
    cout << "\n⏳ Time elapsed: ";
    cout << std::fixed << std::setprecision(2);
    cout << chrono::duration<double>(end_time - start_time).count();
    cout << " seconds" << endl;

    return found_systems;
}

void run_and_debug_benchmark_results(
    const Settings& settings,
    const vector<tuple<string, string, vector<string>, float, float>>& found_systems
) {
    chrono::_V2::system_clock::time_point start_time = chrono::high_resolution_clock::now();

    vector<string> benchmark_sequences;
    ifstream benchmark_in(settings.benchmark_file);
    if (!benchmark_in) {
        cerr << "Error: Could not open benchmark file: " << settings.benchmark_file << endl;
    } else {
        string line;
        while (getline(benchmark_in, line)) {
            if (!line.empty()) {
                benchmark_sequences.push_back(line);
            }
        }
        benchmark_in.close();
        cout << "Loaded " << benchmark_sequences.size() << " benchmark sequences." << endl;
    }

    cout << "  ▸ " << found_systems.size();
    cout << " crispr sequences are found and benchmarked using ";
    cout << benchmark_sequences.size() << " sequences" << endl;

    size_t no_match_count = 0;
    float average_sequence_similarity = 0.0;
    for (const auto& [sequence, repeat, spacers, confidence_cycle_resolution, confidence_topological_sort] : found_systems) {
        // This leads to overestimation by choosing greedy
        const auto expected_sequence = get_most_similar_sequence(
            sequence, benchmark_sequences
        );
        if (expected_sequence == "") {
            cout << "    ▸ No expected match for sequence: ";
            cout << sequence << endl;
            no_match_count++;
            continue;
        }
        
        const auto sequence_similarity = get_string_similarity(
            sequence, expected_sequence
        );
        const auto amount_of_duplicate_spacers = get_number_of_duplicate_spacers(
            spacers, expected_sequence
        );

        cout << "    ▸ ≥" << std::fixed << std::setprecision(2);
        cout << (sequence_similarity * 100) << "% sequence similarity, ";
        cout << "with ";
        cout << spacers.size() << " spacers, ";
        cout << amount_of_duplicate_spacers << " duplicate spacers, confidence of cycle resolution: ";
        cout << (confidence_cycle_resolution * 100) << "%, confidence of topological sort: ";
        cout << (confidence_topological_sort * 100) << "%, and the repeat: ";
        cout << repeat << ", and sequence: ";
        cout << sequence << endl;
        average_sequence_similarity += sequence_similarity;
    }
    average_sequence_similarity /= static_cast<float>(found_systems.size() - no_match_count);

    cout << "  ▸ The average sequence similarity is ";
    cout << std::fixed << std::setprecision(2);
    cout << (average_sequence_similarity * 100);
    cout << "% with " << no_match_count << "/" << found_systems.size();
    cout << " ignored" << endl;

    chrono::_V2::system_clock::time_point end_time = chrono::high_resolution_clock::now();
    cout << "\n⏳ Time elapsed: ";
    cout << std::fixed << std::setprecision(2);
    cout << chrono::duration<double>(end_time - start_time).count();
    cout << " seconds" << endl;
}

void run_and_debug_results(
    const vector<tuple<string, string, vector<string>, float, float>>& found_systems
) {
    cout << "Each result has their own confidence score that can give ";
    cout << "some guidance of how accurate the prediction is." << endl;
    cout << "Take these predictions with a grain of salt:" << endl;
    cout << "  🔴: Many uncertainties, ";
    cout << "e.g. no clear repeat sequence, high spacer contradictions" << endl;
    cout << "  🟠: Some uncertainties, ";
    cout << "e.g. some spacer positions are unclear and were intuitively guessed" << endl;
    cout << "  🟡: Minor uncertainties" << endl;
    cout << "  🟢: Highly confident with the result" << endl;
    cout << endl;
    cout << "----------------------------------------------" << endl;

    int red = 0;
    int orange = 0;
    int yellow = 0;
    int green = 0;
    for (const auto& [sequence, repeat, spacers, confidence_cycle_resolution, confidence_topological_sort] : found_systems) {
        if (repeat.size() <= 23 || confidence_cycle_resolution < 0.5 || confidence_topological_sort < 0.5) {
            cout << "  🔴 ";
            ++red;
        } else if (confidence_cycle_resolution < 0.75 || confidence_topological_sort < 0.75) {
            cout << "  🟠 ";
            ++orange;
        } else if (confidence_cycle_resolution < 0.85 || confidence_topological_sort < 0.85) {
            cout << "  🟡 ";
            ++yellow;
        } else {
            cout << "  🟢 ";
            ++green;
        }

        cout << "repeat: " << repeat << ", sequence: " << sequence << endl;
    }

    const int total = red + orange + yellow + green;
    cout << endl;
    cout << "  ▸ " << found_systems.size() << " CRISPR Arrays were found with ";
    cout << "🔴 (" << red << "/" << total << "), ";
    cout << "🟠 (" << orange << "/" << total << "), ";
    cout << "🟡 (" << yellow << "/" << total << "), ";
    cout << "🟢 (" << green << "/" << total << ")" << endl;
}
