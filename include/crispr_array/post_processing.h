#ifndef POST_PROCESSING_H
#define POST_PROCESSING_H

#include "core/settings.h"
#include "crispr_array/sdbg_repeat_graph.h"
#include "crispr_array/array_output.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <filesystem>
#include <iomanip>
#include <chrono>

// Recover repeat/spacer boundaries from ordered, observed graph paths. SPOA
// consensus and arbitrary shared-k-mer spacer unions do not determine identity.
class PostProcessor {
    Settings& settings;
    std::map<std::string,std::vector<std::string>> crispr_arrays;
    std::vector<mcaat::arrays::System> graph_systems;
    using Path=mcaat::arrays::Path;
    using Node=mcaat::arrays::Node;
    using System=mcaat::arrays::System;
    using Graph=mcaat::arrays::SdbgRepeatGraph;
    std::map<Node,std::vector<Path>> read_cycles(Graph& graph,const std::string& file) {
        std::map<Node,std::vector<Path>> groups;
        auto started=std::chrono::steady_clock::now(),last_report=started;
        size_t cycle_count=0,kmer_count=0;
        if (file.empty()) {
            std::cerr<<"[arrays] Loading in-memory cycle paths from "<<settings.cycles.size()<<" start groups\n";
            for (const auto& [start,cycles]:settings.cycles) {
                for (const auto& cycle:cycles) {
                    groups[start].push_back(graph.vertices(cycle));++cycle_count;kmer_count+=cycle.size();
                }
                auto now=std::chrono::steady_clock::now();
                if (now-last_report>=std::chrono::seconds(5)) {
                    std::cerr<<"[arrays] Loaded "<<cycle_count<<" cycles from "<<groups.size()<<'/'
                             <<settings.cycles.size()<<" start groups\n";last_report=now;
                }
            }
            std::cerr<<"[arrays] Loaded "<<cycle_count<<" cycles in "
                     <<std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count()<<" s\n";
            return groups;
        }
        std::ifstream in(file);
        if (!in) throw std::runtime_error("Cannot read cycles: "+file);
        std::error_code error;
        auto total_bytes=std::filesystem::file_size(file,error);
        std::cerr<<"[arrays] Loading saved cycles: "<<file;
        if (!error) std::cerr<<" ("<<total_bytes/1048576.0<<" MiB)";
        std::cerr<<'\n';
        std::string line;
        size_t line_number=0,bytes=0;
        while (std::getline(in,line)) {
            ++line_number;bytes+=line.size()+1;
            if (line.empty() || line[0]=='#') continue;
            std::istringstream words(line);std::string word,previous;Path cycle;
            Node first_edge=mcaat::arrays::missing_node;
            while (words>>word) {
                ++kmer_count;
                if (cycle.empty()) {
                    Node v=graph.find(word);
                    if (v==mcaat::arrays::missing_node)
                        throw std::runtime_error("Unknown k-mer in cycles, line "+std::to_string(line_number));
                    cycle.push_back(v);
                } else {
                    if (word.size()!=graph.k() || previous.compare(1,graph.k()-1,word,0,graph.k()-1)!=0)
                        throw std::runtime_error("Non-overlapping k-mers in cycles, line "+std::to_string(line_number));
                    auto [edge,vertex]=graph.step(cycle.back(),word.back());
                    if (vertex==mcaat::arrays::missing_node)
                        throw std::runtime_error("Missing graph transition in cycles, line "+std::to_string(line_number));
                    if (cycle.size()==1) first_edge=edge;
                    cycle.push_back(vertex);
                }
                previous.swap(word);
            }
            if (cycle.size()>1 && cycle.front()==cycle.back()) cycle.pop_back();
            if (cycle.size()<2) throw std::runtime_error("Invalid cycle, line "+std::to_string(line_number));
            if (graph.edge(cycle.back(),cycle.front())==mcaat::arrays::missing_node)
                throw std::runtime_error("Non-closing cycle, line "+std::to_string(line_number));
            groups[first_edge].push_back(std::move(cycle));++cycle_count;
            if ((line_number&1023)==0) {
                auto now=std::chrono::steady_clock::now();
                if (now-last_report>=std::chrono::seconds(5)) {
                    std::cerr<<"[arrays] Loaded "<<cycle_count<<" cycles, "<<kmer_count<<" k-mers, "
                             <<bytes/1048576.0<<" MiB";
                    if (!error) std::cerr<<" / "<<total_bytes/1048576.0<<" MiB";
                    std::cerr<<'\n';last_report=now;
                }
            }
        }
        if (in.bad()) throw std::runtime_error("I/O error while reading cycles: "+file);
        std::cerr<<"[arrays] Loaded "<<cycle_count<<" cycles, "<<kmer_count<<" k-mers in "
                 <<std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count()
                 <<" s; "<<graph.index_queries<<" graph index searches, "<<graph.lookup_hits<<" cached starts\n";
        return groups;
    }
public:
    explicit PostProcessor(Settings& s):settings(s) {}
    // Legacy sequence-only view; use graph systems or TSV paths for identity.
    std::map<std::string,std::vector<std::string>> getSystems() const {return crispr_arrays;}
    const std::vector<System>& getGraphSystems() const {return graph_systems;}

    // An explicit cycles file enables re-filtering without re-enumeration. In
    // the pipeline, use settings.cycles even when empty (never stale disk data).
    void run_analysis(const std::string& cycles_file="") {
        using namespace mcaat::arrays;
        if (!settings.sdbg) throw std::runtime_error("Graph must be loaded before array post-processing");
        const auto& dna=settings.dna_sequence_settings;
        if (dna.repeat_min_length<1 || dna.repeat_max_length<1 ||
            dna.spacer_min_length<1 || dna.spacer_max_length<1)
            throw std::invalid_argument("repeat/spacer length limits must be positive");
        Options options{size_t(settings.sdbg->k()),dna.repeat_support,
            size_t(dna.repeat_min_length),size_t(dna.repeat_max_length),
            size_t(dna.spacer_min_length),size_t(dna.spacer_max_length)};
        options.quality=dna.quality;
        Graph graph(*settings.sdbg);
        RepeatFilter<Graph> filter(graph,options);
        auto groups=read_cycles(graph,cycles_file);
        auto filter_start=std::chrono::steady_clock::now(),last_report=filter_start;
        std::string stage;
        graph_systems=filter.run(groups,[&](const char* phase,size_t done,size_t total) {
            auto now=std::chrono::steady_clock::now();
            if (stage!=phase || now-last_report>=std::chrono::seconds(5) || (total && done==total)) {
                std::cerr<<"[arrays] "<<phase<<": "<<done;
                if (total) std::cerr<<'/'<<total<<" groups";
                else std::cerr<<" comparisons";
                std::cerr<<'\n';stage=phase;last_report=now;
            }
        });
        std::cerr<<"[arrays] Filtering completed in "
                 <<std::chrono::duration<double>(std::chrono::steady_clock::now()-filter_start).count()<<" s\n";
        groups.clear(); // Release input paths before serializing retained systems.
        std::cerr<<"[arrays] Writing "<<graph_systems.size()<<" retained systems\n";
        crispr_arrays.clear();
        for (const auto& system:graph_systems) {
            std::set<std::string> variants;
            for (const auto& row:system.observations) variants.insert(row.spacer);
            auto& view=crispr_arrays[filter.spell(system.repeat_path)];
            view.insert(view.end(),variants.begin(),variants.end());
        }
        write_array_outputs(graph,filter,options,graph_systems,settings.output_folder);
    }
};
#endif
