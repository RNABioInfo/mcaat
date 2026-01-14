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
    uint64_t start_node;        // Node containing the start codon
    uint64_t end_node;          // Node containing the stop codon
    int distance_from_repeat;   // Distance in bp from repeat to start codon
    int orf_length;             // Distance in bp from start to stop codon (ORF length)
    
    ORFInfo() : start_node(SDBG::kNullID), end_node(SDBG::kNullID), 
                distance_from_repeat(0), orf_length(0) {}
    ORFInfo(uint64_t s, uint64_t e, int d_repeat, int d_orf) 
        : start_node(s), end_node(e), distance_from_repeat(d_repeat), orf_length(d_orf) {}
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
     * @brief Check if a sequence contains an in-frame stop codon
     * @param sequence The DNA sequence to check  
     * @param frame_offset The reading frame offset (0, 1, or 2)
     * @return true if an in-frame stop codon is found
     */
    bool ContainsInFrameStopCodon(const std::string& sequence, int frame_offset) const;

    /**
     * @brief Scan forward from start_node to find stop codon
     * @param start_node Node where ORF starts
     * @param distance_from_repeat Distance from original repeat to this start node
     * @param max_orf_length Maximum ORF length to search
     * @return ORFInfo if stop codon found, empty optional otherwise
     */
    std::optional<ORFInfo> ScanForStopCodon(
        uint64_t start_node,
        int distance_from_repeat,
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
     * If min_orf_length is specified, only ORFs with length >= min_orf_length
     * are returned. Shorter ORFs are skipped and search continues.
     * 
     * The BFS stores only the last layer of neighbors to minimize memory usage.
     * 
     * @param repeat_node The starting repeat node
     * @param min_distance Minimum distance in bp from repeat to start codon
     * @param max_distance Maximum distance in bp from repeat to start codon
     * @param min_orf_length Minimum ORF length in bp (default 0, no filter)
     * @return ORFInfo containing start_node, end_node, and distance if found,
     *         or empty optional if no ORF found
     */
    std::optional<ORFInfo> FindFirstORF(
        uint64_t repeat_node,
        int min_distance,
        int max_distance,
        int min_orf_length = 0
    ) const;
    
    /**
     * @brief Find ALL ORFs starting at distance D from a repeat node
     * 
     * This function performs a BFS from the repeat node to find ALL nodes at distance
     * [min_distance, max_distance] that contain start codons (ATG, GTG, TTG).
     * For each start codon found, it scans forward to find the stop codon
     * and returns ALL complete ORFs found.
     * 
     * If min_orf_length is specified, only ORFs with length >= min_orf_length
     * are included. Shorter ORFs are filtered out.
     * 
     * @param repeat_node The starting repeat node
     * @param min_distance Minimum distance in bp from repeat to start codon
     * @param max_distance Maximum distance in bp from repeat to start codon
     * @param min_orf_length Minimum ORF length in bp (default 0, no filter)
     * @param max_traverse Maximum nodes to traverse (default 5000)
     * @return Vector of all ORFInfo found (sorted by distance from repeat)
     */
    std::vector<ORFInfo> FindAllORFs(
        uint64_t repeat_node,
        int min_distance,
        int max_distance,
        int min_orf_length = 0,
        int max_traverse = 5000
    ) const;
};

#endif  // INCLUDE_ORF_FINDER_H_
