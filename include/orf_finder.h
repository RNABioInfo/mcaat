/**
 * @file orf_finder.h
 * @brief ORF (Open Reading Frame) finder for sequence De Bruijn graphs
 * 
 * This module finds ORFs in SDBG starting from repeat nodes.
 * It uses BFS to locate start codons at a specific distance range,
 * then scans forward to find the corresponding stop codons.
 * 
 * Start codons: ATG, GTG, TTG
 * Stop codons: TAA, TAG, TGA
 */

#ifndef INCLUDE_ORF_FINDER_H_
#define INCLUDE_ORF_FINDER_H_

#include <cstdint>
#include <vector>
#include <optional>
#include <string>
#include "sdbg/sdbg.h"

/**
 * @brief Information about a discovered ORF
 */
struct ORFInfo {
    uint64_t start_node;  // Node containing the start codon
    uint64_t end_node;    // Node containing the stop codon
    int distance;         // Distance in bp between start and end (ORF length)
    
    ORFInfo() : start_node(SDBG::kNullID), end_node(SDBG::kNullID), distance(0) {}
    ORFInfo(uint64_t s, uint64_t e, int d) : start_node(s), end_node(e), distance(d) {}
};

/**
 * @brief ORF finder class for discovering Open Reading Frames in SDBG
 * 
 * This class provides functionality to find ORFs starting from repeat nodes
 * in a sequence De Bruijn graph. It uses BFS-based search strategies to
 * locate start codons at specific distances and scan for stop codons.
 */
class ORFFinder {
private:
    const SDBG& sdbg;
    static constexpr int MAX_ORF_LENGTH = 5000;  // Maximum ORF length in bp

    /**
     * @brief Convert node sequence to DNA string
     * @param node_id The node to convert
     * @return DNA string representation of the node's k-mer
     */
    std::string NodeToSequence(uint64_t node_id) const;

    /**
     * @brief Check if a sequence contains a start codon (ATG, GTG, TTG)
     * @param sequence The DNA sequence to check
     * @return true if a start codon is found
     */
    bool ContainsStartCodon(const std::string& sequence) const;

    /**
     * @brief Check if a sequence contains a stop codon (TAA, TAG, TGA)
     * @param sequence The DNA sequence to check
     * @return true if a stop codon is found
     */
    bool ContainsStopCodon(const std::string& sequence) const;

    /**
     * @brief Scan forward from start_node to find stop codon
     * @param start_node Node where ORF starts
     * @param max_orf_length Maximum ORF length to search
     * @return ORFInfo if stop codon found, empty optional otherwise
     */
    std::optional<ORFInfo> ScanForStopCodon(
        uint64_t start_node,
        int max_orf_length
    ) const;

public:
    /**
     * @brief Constructor
     * @param sdbg Reference to the sequence De Bruijn graph
     */
    explicit ORFFinder(const SDBG& sdbg);

    /**
     * @brief Find the first ORF starting at distance D from a repeat node
     * 
     * This function performs a BFS from the repeat node to find nodes at distance
     * [min_distance, max_distance] that contain start codons (ATG, GTG, TTG).
     * For each start codon found, it scans forward to find the first stop codon
     * (TAA, TAG, TGA) and returns information about the first complete ORF.
     * 
     * The BFS stores only the last layer of neighbors to minimize memory usage.
     * 
     * @param repeat_node The starting repeat node
     * @param min_distance Minimum distance in bp from repeat to start codon
     * @param max_distance Maximum distance in bp from repeat to start codon
     * @return ORFInfo containing start_node, end_node, and distance if found,
     *         or empty optional if no ORF found
     */
    std::optional<ORFInfo> FindFirstORF(
        uint64_t repeat_node,
        int min_distance,
        int max_distance
    ) const;
};

#endif  // INCLUDE_ORF_FINDER_H_
