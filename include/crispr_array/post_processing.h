#ifndef POST_PROCESSING_H
#define POST_PROCESSING_H

#include "core/settings.h"
#include "crispr_array/sdbg_repeat_graph.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <filesystem>
#include <iomanip>

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
    static std::ofstream output(const std::filesystem::path& path) {
        std::ofstream out(path);
        out.exceptions(std::ios::badbit | std::ios::failbit);
        return out;
    }
    std::map<Node,std::vector<Path>> read_cycles(Graph& graph,const std::string& file) {
        std::map<Node,std::vector<Path>> groups;
        if (file.empty()) {
            for (const auto& [start,cycles]:settings.cycles)
                for (const auto& cycle:cycles) groups[start].push_back(graph.vertices(cycle));
            return groups;
        }
        std::ifstream in(file);
        if (!in) throw std::runtime_error("Cannot read cycles: "+file);
        std::string line;
        size_t line_number=0;
        while (std::getline(in,line)) {
            ++line_number;
            if (line.empty() || line[0]=='#') continue;
            std::istringstream words(line);std::string word;Path cycle;
            while (words>>word) {
                Node v=graph.find(word);
                if (v==mcaat::arrays::missing_node)
                    throw std::runtime_error("Unknown k-mer in cycles, line "+std::to_string(line_number));
                cycle.push_back(v);
            }
            if (cycle.size()>1 && cycle.front()==cycle.back()) cycle.pop_back();
            if (cycle.size()<2) throw std::runtime_error("Invalid cycle, line "+std::to_string(line_number));
            auto edges=graph.edges(cycle); // Requires every transition, including closure.
            groups[edges.front()].push_back(std::move(cycle));
        }
        return groups;
    }
    static bool owned_file(const std::string& name) {
        if (name.rfind("CRISPR_Arrays_",0)==0 && name.size()>18 &&
            name.substr(name.size()-4)==".txt") {
            auto number=name.substr(14,name.size()-18);
            return !number.empty() && number.find_first_not_of("0123456789")==std::string::npos;
        }
        return name=="CRISPR_graph_paths.tsv" || name=="CRISPR_graph_systems.tsv" ||
               name=="CRISPR_graph_candidates.tsv" || name=="CRISPR_graph_summary.tsv";
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
        namespace fs=std::filesystem;
        if (!settings.sdbg) throw std::runtime_error("Graph must be loaded before array post-processing");
        const auto& dna=settings.dna_sequence_settings;
        if (dna.repeat_min_length<1 || dna.repeat_max_length<1 ||
            dna.spacer_min_length<1 || dna.spacer_max_length<1)
            throw std::invalid_argument("repeat/spacer length limits must be positive");
        Options options{size_t(settings.sdbg->k()),dna.repeat_support,
            size_t(dna.repeat_min_length),size_t(dna.repeat_max_length),
            size_t(dna.spacer_min_length),size_t(dna.spacer_max_length)};
        Graph graph(*settings.sdbg);
        RepeatFilter<Graph> filter(graph,options);
        auto groups=read_cycles(graph,cycles_file);
        graph_systems=filter.run(groups);
        const auto& systems=graph_systems;
        crispr_arrays.clear();
        fs::create_directories(settings.output_folder);
        const fs::path destination(settings.output_folder);
        const fs::path staging=destination/(".mcaat-postprocess-"+std::to_string(
            std::chrono::steady_clock::now().time_since_epoch().count()));
        fs::create_directory(staging);
        try {
            auto paths=output(staging/"CRISPR_graph_paths.tsv");
            paths<<"array_id\tobservation_id\trepeat\tspacer\trepeat_vertices\tspacer_vertices"
                   "\trepeat_to_spacer_vertices\tspacer_to_repeat_vertices\tcycle_vertices"
                   "\tcycle_edges\tsupports_repeat_path\n";
            auto summary=output(staging/"CRISPR_graph_systems.tsv");
            summary<<"array_id\trepeat\tsupport_cycles\tcandidate_cycles\tassigned_cycles"
                     "\trepeat_support_threshold\trepeat_vertices\tsource_start_edges\talternative_repeat_paths\n";
            auto candidates=output(staging/"CRISPR_graph_candidates.tsv");
            candidates<<"array_id\tcycle_vertices\tboundary_assignment\n";
            auto arrays=output(staging/"CRISPR_Arrays_1.txt");
            auto header=[&]() {
                arrays<<"# MCAAT graph-based CRISPR array candidates\n"
                      <<"# Repeat support: "<<std::setprecision(6)<<options.repeat_support
                      <<" (distinct cycles per candidate group; not sequence identity)\n"
                      <<"# Spacer rows have no inferred genomic order. See CRISPR_graph_*.tsv for paths.\n\n";
            };
            header();
            size_t array_number=0,file_number=1,lines=4,spacer_count=0,unresolved=0;
            for (const auto& system:systems) {
                const auto repeat=filter.spell(system.repeat_path);
                const auto id="Array_"+std::to_string(++array_number);
                std::set<std::string> unique;
                for (const auto& row:system.observations) unique.insert(row.spacer);
                if (lines+unique.size()+3>10000 && lines>4) {
                    arrays.close();arrays=output(staging/("CRISPR_Arrays_"+std::to_string(++file_number)+".txt"));
                    header();lines=4;
                }
                arrays<<'>'<<id<<"  spacers="<<unique.size()<<'\n'<<repeat<<'\n';
                unique.clear();size_t observation=0;
                std::set<Path> assigned;
                for (const auto& row:system.observations) {
                    assigned.insert(cycle_key(row.cycle));
                    if (unique.insert(row.spacer).second) {
                        arrays<<"        "<<(row.repeat==repeat?std::string(repeat.size(),'-'):row.repeat)
                              <<'\t'<<row.spacer<<'\n';
                        crispr_arrays[repeat].push_back(row.spacer);
                        ++spacer_count;
                    }
                    const auto parts=partition(row.cycle,row.repeat_nodes,options.k);
                    paths<<id<<'\t'<<++observation<<'\t'<<row.repeat<<'\t'<<row.spacer<<'\t'
                         <<path_text(parts.repeat)<<'\t'<<path_text(parts.spacer)<<'\t'
                         <<path_text(parts.repeat_to_spacer)<<'\t'<<path_text(parts.spacer_to_repeat)<<'\t'
                         <<path_text(row.cycle)<<'\t'<<path_text(graph.edges(row.cycle))<<'\t'
                         <<(row.supports_path?1:0)<<'\n';
                }
                arrays<<'\n';lines+=unique.size()+3;
                summary<<id<<'\t'<<repeat<<'\t'<<system.support_count<<'\t'
                       <<system.candidate_cycles.size()<<'\t'<<system.observations.size()<<'\t'
                       <<std::setprecision(17)<<options.repeat_support<<'\t'<<path_text(system.repeat_path)<<'\t'
                       <<path_text(system.sources)<<'\t'<<system.alternate_repeat_paths<<'\n';
                for (const auto& c:system.candidate_cycles) {
                    bool located=assigned.count(cycle_key(c));unresolved+=!located;
                    candidates<<id<<'\t'<<path_text(c)<<'\t'<<(located?"assigned":"unresolved")<<'\n';
                }
            }
            auto stats=output(staging/"CRISPR_graph_summary.tsv");
            const auto& s=filter.stats;
            stats<<"metric\tvalue\nrepeat_support\t"<<std::setprecision(17)<<options.repeat_support
                 <<"\nk\t"<<options.k<<"\ninput_cycles\t"<<s.input_cycles
                 <<"\ninvalid_cycles\t"<<s.invalid_cycles<<"\nduplicate_cycles_within_groups\t"<<s.duplicate_cycles
                 <<"\ninput_groups\t"<<s.input_groups<<"\ngroups_without_repeat\t"<<s.groups_without_repeat
                 <<"\nduplicate_systems_removed\t"<<s.duplicate_systems<<"\noverlapping_systems_merged\t"<<s.merged_systems
                 <<"\narrays\t"<<systems.size()<<"\nspacers\t"<<spacer_count<<"\nunresolved_candidate_cycles\t"<<unresolved<<'\n';
            paths.close();summary.close();candidates.close();arrays.close();stats.close();
            // Stage successful output before replacing only files owned by this processor.
            for (const auto& item:fs::directory_iterator(destination))
                if (item.is_regular_file() && owned_file(item.path().filename().string())) fs::remove(item.path());
            for (const auto& item:fs::directory_iterator(staging)) fs::rename(item.path(),destination/item.path().filename());
            fs::remove(staging);
            std::ostringstream support_text;support_text<<std::setprecision(6)<<options.repeat_support;
            std::cout<<"  ▸ Graph filtering: "<<systems.size()<<" arrays, "<<spacer_count<<" spacers; "
                     <<s.duplicate_systems<<" duplicates removed, "<<s.merged_systems<<" overlaps merged\n"
                     <<"  ▸ Repeat support: "<<support_text.str()<<"; invalid cycles: "<<s.invalid_cycles
                     <<"; unresolved boundaries: "<<unresolved<<std::endl;
        } catch (...) {std::error_code ec;fs::remove_all(staging,ec);throw;}
    }
};
#endif
