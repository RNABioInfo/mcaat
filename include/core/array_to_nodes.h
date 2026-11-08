#ifndef ARRAY_TO_NODES_H
#define ARRAY_TO_NODES_H

/**
 * @file array_to_nodes.h
 * @brief Bridge from MCAAT CRISPR array output files to SDBG node paths.
 *
 * Prefers validated CRISPR_graph_paths.tsv observations produced by PostProcessor.
 * For legacy runs, parses CRISPR_Arrays_N.txt and maps
 * each repeat consensus and spacer sequence back to SDBG node IDs via
 * SDBG::IndexBinarySearch, producing CRISPRPostProcessor::FilteredArray
 * structs ready for CasWorkflow::DetectAllCassettesFromFiltered.
 *
 * Output file format (written by PostProcessor):
 *   # comment lines
 *   >Array_N  spacers=M
 *   REPEATCONSENSUS
 *           REPEAT_VARIANT<TAB>SPACER
 *           REPEAT_VARIANT<TAB>SPACER
 *   (blank line)
 *   >Array_N+1 ...
 */

#include "cas/crispr_postprocessor.h"
#include "sdbg/sdbg.h"
#include "crispr_array/sdbg_repeat_graph.h"
#include <string>
#include <vector>
#include <fstream>
#include <filesystem>
#include <algorithm>
#include <cstdint>

// ── Encoding ─────────────────────────────────────────────────────────────────

/**
 * @brief Encode a DNA character to SDBG alphabet (A=1, C=2, G=3, T=4).
 * @return 0 for any non-ACGT character (signals an invalid base).
 */
inline uint8_t EncodeBase(char c) {
    switch (c) {
        case 'A': case 'a': return 1;
        case 'C': case 'c': return 2;
        case 'G': case 'g': return 3;
        case 'T': case 't': return 4;
        default:             return 0;
    }
}

// ── Core mapping ─────────────────────────────────────────────────────────────

/**
 * @brief Map a DNA sequence to an ordered list of SDBG node IDs.
 *
 * Slides a k-mer window of width sdbg.k() across @p seq, encodes each
 * k-mer, and calls SDBG::IndexBinarySearch to find the corresponding
 * edge/node.  K-mers that contain non-ACGT bases or are absent from the
 * graph are silently skipped (path may have gaps).
 *
 * @param seq  DNA string (any case, ACGT only).
 * @param sdbg The succinct de Bruijn graph to query.
 * @return     Vector of valid node IDs, one per found k-mer.
 */
inline std::vector<uint64_t> SequenceToNodePath(const std::string& seq,
                                                 const SDBG& sdbg) {
    const uint32_t k = sdbg.k();
    std::vector<uint64_t> path;
    if (seq.size() < k) return path;

    path.reserve(seq.size() - k + 1);
    std::vector<uint8_t> encoded(k);

    for (size_t i = 0; i + k <= seq.size(); ++i) {
        bool valid = true;
        for (uint32_t j = 0; j < k; ++j) {
            uint8_t c = EncodeBase(seq[i + j]);
            if (c == 0) { valid = false; break; }
            encoded[j] = c;
        }
        if (!valid) continue;

        uint64_t node = sdbg.IndexBinarySearch(encoded.data());
        if (node != SDBG::kNullID && sdbg.IsValidEdge(node)) {
            path.push_back(node);
        }
    }
    return path;
}

// ── File parser ───────────────────────────────────────────────────────────────

/**
 * @brief Parse one CRISPR_Arrays_N.txt file and populate FilteredArray structs.
 *
 * For each array in the file:
 *   - repeat_sequence  = consensus line (non-indented line after '>Array_N')
 *   - repeat_path      = SequenceToNodePath(repeat_sequence)
 *   - spacer_node_paths= one entry per spacer (tab-separated field on indented lines)
 *
 * Arrays whose repeat_path is empty (e.g. consensus not in the graph) are
 * included with an empty repeat_path so CasWorkflow can decide what to do.
 *
 * @param arrays_file Path to a single CRISPR_Arrays_N.txt file.
 * @param sdbg        The succinct de Bruijn graph.
 * @return            Vector of FilteredArray, one per array in the file.
 */
inline std::vector<CRISPRPostProcessor::FilteredArray>
BuildFilteredArrays(const std::string& arrays_file, const SDBG& sdbg) {
    std::vector<CRISPRPostProcessor::FilteredArray> result;

    std::ifstream f(arrays_file);
    if (!f.is_open()) return result;

    CRISPRPostProcessor::FilteredArray current;
    bool in_array    = false;
    bool repeat_set  = false;

    // Flush accumulated array into result (even if repeat_path is empty)
    auto flush = [&]() {
        if (in_array) {
            result.push_back(std::move(current));
        }
        current      = CRISPRPostProcessor::FilteredArray{};
        in_array     = false;
        repeat_set   = false;
    };

    std::string line;
    while (std::getline(f, line)) {
        // Strip trailing CR/LF
        while (!line.empty() && (line.back() == '\r' || line.back() == '\n'))
            line.pop_back();

        if (line.empty())  continue;
        if (line[0] == '#') continue;   // comment / file header

        if (line[0] == '>') {
            // Start of a new array block
            flush();
            in_array   = true;
            repeat_set = false;
            continue;
        }

        if (!in_array) continue;

        bool is_indented = (line[0] == ' ' || line[0] == '\t');

        if (!repeat_set && !is_indented) {
            // First non-indented, non-header line = repeat consensus
            // Trim any trailing whitespace
            size_t end = line.find_last_not_of(" \t");
            current.repeat_sequence = (end != std::string::npos)
                                      ? line.substr(0, end + 1)
                                      : line;
            current.repeat_path = SequenceToNodePath(current.repeat_sequence, sdbg);
            repeat_set = true;
        } else if (repeat_set && is_indented) {
            // Indented line: "        REPEAT_VARIANT<TAB>SPACER"
            // The spacer is always the second tab-delimited field.
            size_t tab = line.find('\t');
            if (tab == std::string::npos) continue;  // malformed line

            std::string spacer = line.substr(tab + 1);
            // Trim any trailing whitespace from spacer
            size_t s_end = spacer.find_last_not_of(" \t\r\n");
            if (s_end != std::string::npos) spacer = spacer.substr(0, s_end + 1);
            if (spacer.empty()) continue;

            current.spacer_node_paths.push_back(
                SequenceToNodePath(spacer, sdbg));
        }
    }
    flush();

    return result;
}

// Read exact observed edge paths when graph-based post-processing supplied them.
// Numeric IDs alone are not trusted: validate every label, edge and boundary.
inline std::vector<CRISPRPostProcessor::FilteredArray>
BuildFilteredArraysFromGraphPaths(const std::string& file,const SDBG& sdbg) {
    using namespace mcaat::arrays;
    std::ifstream input(file);
    if (!input) throw std::runtime_error("Cannot open graph path table: "+file);
    const std::string expected="array_id\tobservation_id\trepeat\tspacer\trepeat_vertices\tspacer_vertices"
        "\trepeat_to_spacer_vertices\tspacer_to_repeat_vertices\tcycle_vertices\tcycle_edges\tsupports_repeat_path";
    std::string line;
    if (!std::getline(input,line) || line!=expected) throw std::runtime_error("Unsupported graph path table header");
    SdbgRepeatGraph graph(sdbg);
    Options options;options.k=sdbg.k();
    RepeatFilter<SdbgRepeatGraph> spelling(graph,options);
    std::vector<CRISPRPostProcessor::FilteredArray> result;
    std::map<std::string,size_t> indices;
    std::map<std::string,std::set<std::string>> seen_spacers;
    while (std::getline(input,line)) {
        std::vector<std::string> fields;
        std::istringstream stream(line);std::string field;
        while(std::getline(stream,field,'\t')) fields.push_back(field);
        if(fields.size()!=11) throw std::runtime_error("Malformed graph path table row");
        const auto& id=fields[0];const auto& repeat=fields[2];const auto& spacer=fields[3];
        auto vertices=parse_path(fields[8]),edges=parse_path(fields[9]);
        if(repeat.size()<sdbg.k() || vertices.size()!=repeat.size()+spacer.size() ||
           vertices.size()!=std::set<Node>(vertices.begin(),vertices.end()).size() ||
           graph.vertices(edges)!=vertices || spelling.spell(vertices).substr(0,vertices.size())!=repeat+spacer)
            throw std::runtime_error("Graph path table does not match graph or sequences: "+id);
        size_t repeat_nodes=repeat.size()-sdbg.k()+1;
        auto parts=partition(vertices,repeat_nodes,sdbg.k());
        if(parse_path(fields[4])!=parts.repeat || parse_path(fields[5])!=parts.spacer ||
           parse_path(fields[6])!=parts.repeat_to_spacer || parse_path(fields[7])!=parts.spacer_to_repeat ||
           (fields[10]!="0" && fields[10]!="1"))
            throw std::runtime_error("Invalid repeat/spacer boundary assignment: "+id);
        if(!indices.count(id)) {indices[id]=result.size();result.emplace_back();}
        auto& array=result[indices.at(id)];
        if(fields[10]=="1") {
            if(!array.repeat_sequence.empty() && array.repeat_sequence!=repeat)
                throw std::runtime_error("Inconsistent supported repeat: "+id);
            if(array.repeat_path.empty()) {
                array.repeat_sequence=repeat;
                array.repeat_path.assign(edges.begin(),edges.begin()+repeat_nodes);
            }
        }
        if(seen_spacers[id].insert(spacer).second) {
            std::vector<Node> path;
            for(size_t i=repeat.size();i+sdbg.k()<=edges.size();++i) path.push_back(edges[i]);
            array.spacer_node_paths.push_back(std::move(path));
        }
    }
    for(const auto& a:result) if(a.repeat_path.empty() || a.spacer_node_paths.size()<2)
        throw std::runtime_error("Graph path table contains an unsupported or singleton array");
    return result;
}

// ── Directory scan ────────────────────────────────────────────────────────────

/**
 * @brief Collect FilteredArrays from all CRISPR_Arrays_*.txt files in a directory.
 *
 * Files are processed in lexicographic order so array numbering is
 * deterministic across runs.
 *
 * @param output_dir Directory that contains the CRISPR_Arrays_N.txt files
 *                   (typically Settings::output_folder).
 * @param sdbg       The succinct de Bruijn graph.
 * @return           Combined vector of all FilteredArray structs found.
 */
inline std::vector<CRISPRPostProcessor::FilteredArray>
BuildFilteredArraysFromDir(const std::string& output_dir, const SDBG& sdbg) {
    namespace fs = std::filesystem;

    std::vector<CRISPRPostProcessor::FilteredArray> all;
    if (!fs::exists(output_dir) || !fs::is_directory(output_dir)) return all;

    const fs::path paths=fs::path(output_dir)/"CRISPR_graph_paths.tsv";
    if (fs::exists(paths)) return BuildFilteredArraysFromGraphPaths(paths.string(),sdbg);

    // Collect matching files
    std::vector<fs::path> files;
    for (const auto& entry : fs::directory_iterator(output_dir)) {
        if (!entry.is_regular_file()) continue;
        const std::string name = entry.path().filename().string();
        if (name.find("CRISPR_Arrays_") != std::string::npos &&
            entry.path().extension() == ".txt") {
            files.push_back(entry.path());
        }
    }
    std::sort(files.begin(), files.end());

    for (const auto& p : files) {
        auto batch = BuildFilteredArrays(p.string(), sdbg);
        for (auto& fa : batch)
            all.push_back(std::move(fa));
    }

    return all;
}

#endif  // ARRAY_TO_NODES_H
