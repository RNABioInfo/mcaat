#ifndef INCLUDE_READS_H_
#define INCLUDE_READS_H_

#include <iostream>
#include <cstdint>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <filesystem>
#include <optional>

#include "sdbg/sdbg.h"
#include "sequence/io/sequence_lib.h"
#include "idba/sequence.h"
#include "kseq++/seqio.hpp"

using namespace std;

/**
 * @brief Reads all the sequences in the fastq file
 * 
 * Uses kseqpp to parse the fastq files and filters out only the sequences
 * 
 * @param fastq_file_path
 * @return sequences (vector<string>)
 */
vector<string> extract_sequences_from_fastq_file(const string& fastq_file_path);

/**
 * @internal
 * @brief Reverses the string and flips the nucleotides
 * 
 * The String is reversed and every 'A' - 'T' and 'C' - 'G' are flipped.
 * Any other character is not modified
 * 
 * @param sequence
 */
void reverse_pair_ends_sequence(string& sequence);

/**
 * @internal
 * @brief Get the read from sequence
 * 
 * for each k_mer contained in sequence the node id is found and put into the result in order
 * 
 * @param sdbg Used to convert k_mer to node id
 * @param nodes_of_cycles filtering of reads
 * @param sequence Used to get node id chain
 * 
 * @return {}, if sequence <= 2*K || start and end node_id are not in nodes_of_cycles
 * @return node id chain (vector<uint64_t>), else
 */
vector<uint64_t> get_read_from_sequence(
    const SDBG& sdbg,
    const unordered_set<uint64_t>& nodes_of_cycles,
    const string& sequence
);

/**
 * @brief Get the reads that are relevant in a compressed lossless format
 * 
 * Takes any read where the start or the end is part of the cycles
 * 
 * @pre fastq files provided have to be valid
 * 
 * @param sdbg The graph used to map sequence k-mer to node id
 * @param fastq_file_1 Normal fastq read file
 * @param fastq_file_2 Optional second fastq file with reversed complement sequences
 * @param cycles 
 * @param biggest_jump_size 
 * 
 * @return Reads that start or end in the cycles (vector<uint64_t>)
 */
vector<vector<uint64_t>> get_reads(
    const SDBG& sdbg,
    const string& fastq_file_1,
    const optional<string>& fastq_file_2,
    const vector<vector<uint64_t>> cycles
);

#endif
