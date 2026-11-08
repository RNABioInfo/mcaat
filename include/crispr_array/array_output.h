#ifndef MCAAT_ARRAY_OUTPUT_H
#define MCAAT_ARRAY_OUTPUT_H
#include "crispr_array/graph_repeat_filter.h"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <chrono>

namespace mcaat::arrays {
inline std::ofstream array_output(const std::filesystem::path& path) {
    std::ofstream out(path);out.exceptions(std::ios::badbit|std::ios::failbit);return out;
}
inline bool owned_array_file(const std::string& name) {
    if (name.rfind("CRISPR_Arrays_",0)==0 && name.size()>18 && name.substr(name.size()-4)==".txt") {
        auto number=name.substr(14,name.size()-18);
        return !number.empty() && number.find_first_not_of("0123456789")==std::string::npos;
    }
    static const std::set<std::string> names={"CRISPR_graph_paths.tsv","CRISPR_graph_systems.tsv",
        "CRISPR_graph_candidates.tsv","CRISPR_graph_summary.tsv","CRISPR_spacer_variants.tsv",
        "CRISPR_graph_decisions.tsv","CRISPR_graph_alternatives.tsv","CRISPR_graph_review_cycles.tsv","CRISPR_graph_shared_cycles.tsv"};
    return names.count(name);
}

template<class Graph>
void write_array_outputs(Graph& graph,RepeatFilter<Graph>& filter,const Options& options,
                         const std::vector<System>& systems,const std::string& folder) {
    namespace fs=std::filesystem;
    fs::create_directories(folder);const fs::path destination(folder);
    const fs::path staging=destination/(".mcaat-postprocess-"+std::to_string(
        std::chrono::steady_clock::now().time_since_epoch().count()));
    fs::create_directory(staging);
    try {
        auto flags_text=[](const Quality& q) {
            std::string text;for (const auto& flag:q.flags) {if (!text.empty()) text+=',';text+=flag;}
            return text.empty()?std::string("-"):text;
        };
        auto paths=array_output(staging/"CRISPR_graph_paths.tsv");
        const std::string path_header="observation_id\trepeat\tspacer\trepeat_vertices\tspacer_vertices"
            "\trepeat_to_spacer_vertices\tspacer_to_repeat_vertices\tcycle_vertices\tcycle_edges\tsupports_repeat_path\n";
        paths<<"array_id\t"<<path_header;
        auto emit_path=[&](std::ostream& out,const Observation& row,size_t observation) {
            const auto parts=partition(row.cycle,row.repeat_nodes,options.k);
            out<<observation<<'\t'<<row.repeat<<'\t'<<row.spacer<<'\t'
               <<path_text(parts.repeat)<<'\t'<<path_text(parts.spacer)<<'\t'
               <<path_text(parts.repeat_to_spacer)<<'\t'<<path_text(parts.spacer_to_repeat)<<'\t'
               <<path_text(row.cycle)<<'\t'<<path_text(graph.edges(row.cycle))<<'\t'<<(row.supports_path?1:0)<<'\n';
        };
        auto summary=array_output(staging/"CRISPR_graph_systems.tsv");
        summary<<"array_id\trepeat\tsupport_cycles\tcandidate_cycles\tassigned_cycles\trepeat_support_threshold"
                 "\trepeat_vertices\tsource_start_edges\talternative_repeat_paths\tcandidate_id\tspacer_families\tsequence_variants\tquality_flags\n";
        auto candidates=array_output(staging/"CRISPR_graph_candidates.tsv");
        candidates<<"array_id\tcycle_vertices\tboundary_assignment\n";
        auto variants=array_output(staging/"CRISPR_spacer_variants.tsv");
        variants<<"array_id\tfamily_id\trepresentative\tvariant\tis_representative\n";
        auto arrays=array_output(staging/"CRISPR_Arrays_1.txt");
        auto header=[&]() {
            arrays<<"# MCAAT graph-based CRISPR array candidates\n"
                  <<"# Repeat support: "<<std::setprecision(6)<<options.repeat_support
                  <<" (distinct cycles per candidate group; not sequence identity)\n"
                  <<"# Rows: all distinct observed spacers, without inferred genomic order. Similarity families and paths are in TSV files.\n\n";
        };header();
        size_t array_number=0,file_number=1,lines=4,families=0,variant_count=0,observations=0,unresolved=0;
        std::map<size_t,std::pair<std::string,const System*>> accepted;
        std::map<std::string,std::map<std::string,size_t>> cycle_arrays;
        std::set<std::string> distinct_spacers;
        for (const auto& system:systems) {
            const auto repeat=filter.spell(system.repeat_path),id="Array_"+std::to_string(++array_number);
            accepted.emplace(system.candidate_id,std::make_pair(id,&system));
            if (lines+system.quality.variants+3>10000 && lines>4) {
                arrays.close();arrays=array_output(staging/("CRISPR_Arrays_"+std::to_string(++file_number)+".txt"));header();lines=4;
            }
            arrays<<'>'<<id<<"  spacers="<<system.quality.variants<<"  spacer_families="<<system.quality.families.size()
                  <<"  flags="<<flags_text(system.quality)<<'\n'<<repeat<<'\n';
            std::map<std::string,const Observation*> observed;
            for (const auto& row:system.observations) observed.emplace(row.spacer,&row);
            for (const auto& [spacer,row]:observed)
                arrays<<"        "<<(row->repeat==repeat?std::string(repeat.size(),'-'):row->repeat)<<'\t'<<spacer<<'\n';
            size_t family_id=0;
            for (const auto& family:system.quality.families) {
                if (!observed.count(family.representative)) throw std::runtime_error("spacer representative is not observed");
                ++family_id;
                for (const auto& variant:family.variants)
                    variants<<id<<'\t'<<family_id<<'\t'<<family.representative<<'\t'<<variant<<'\t'<<(variant==family.representative?1:0)<<'\n';
            }
            arrays<<'\n';lines+=system.quality.variants+3;
            families+=system.quality.families.size();variant_count+=system.quality.variants;observations+=system.observations.size();
            std::set<Path> assigned;size_t observation=0;
            for (const auto& row:system.observations) {
                assigned.insert(cycle_key(row.cycle));paths<<id<<'\t';emit_path(paths,row,++observation);
                cycle_arrays[circular_sequence_key(row.repeat+row.spacer)].emplace(id,system.candidate_id);
                distinct_spacers.insert(std::min(row.spacer,reverse_complement(row.spacer)));
            }
            summary<<id<<'\t'<<repeat<<'\t'<<system.support_count<<'\t'<<system.candidate_count<<'\t'<<system.observations.size()<<'\t'
                   <<std::setprecision(17)<<options.repeat_support<<'\t'<<path_text(system.repeat_path)<<'\t'<<path_text(system.sources)<<'\t'
                   <<system.alternate_repeat_paths<<'\t'<<system.candidate_id<<'\t'<<system.quality.families.size()<<'\t'<<system.quality.variants
                   <<'\t'<<flags_text(system.quality)<<'\n';
            unresolved+=system.candidate_count-system.observations.size();
            for (const auto& c:system.candidate_cycles)
                candidates<<id<<'\t'<<path_text(c)<<'\t'<<(assigned.count(cycle_key(c))?"assigned":"unresolved")<<'\n';
        }
        auto shared=array_output(staging/"CRISPR_graph_shared_cycles.tsv");
        shared<<"canonical_circular_sequence\tarray_ids\tcandidate_ids\n";
        size_t shared_count=0;
        for (const auto& [sequence,owners]:cycle_arrays) if (owners.size()>1) {
            ++shared_count;std::string names;Path ids;
            for (const auto& [name,id]:owners) {if (!names.empty()) names+=',';names+=name;ids.push_back(id);}
            shared<<sequence<<'\t'<<names<<'\t'<<path_text(ids)<<'\n';
        }
        auto decisions=array_output(staging/"CRISPR_graph_decisions.tsv");
        decisions<<"candidate_id\tstatus\tarray_id\trepeat\toriginal_cycles\toriginal_assigned\tretained_assigned\tremoved_shared_cycles"
                   "\tcompeting_candidates\tspacer_families\tsequence_variants\tmean_spacer_similarity\tmean_repeat_spacer_similarity\tsource_start_edges\tremoved_strand_duplicates\tquality_flags\tcycles_outside_length_range\tobservations_outside_repeat_spacer_lengths\n";
        auto alternatives=array_output(staging/"CRISPR_graph_alternatives.tsv");
        alternatives<<"candidate_id\treason\t"<<path_header;
        auto review=array_output(staging/"CRISPR_graph_review_cycles.tsv");
        review<<"candidate_id\tstatus\tcycle_vertices\tboundary_assignment\n";
        size_t review_observations=0,known_review_cycles=0;
        for (const auto& d:filter.decisions) {
            const auto& s=d.original;auto found=accepted.find(s.candidate_id);
            Path competing(d.competing_candidates.begin(),d.competing_candidates.end());
            decisions<<s.candidate_id<<'\t'<<d.status<<'\t'<<(found==accepted.end()?"-":found->second.first)<<'\t'
                <<filter.spell(s.repeat_path)<<'\t'<<s.candidate_count<<'\t'<<s.observations.size()<<'\t'
                <<(found==accepted.end()?0:found->second.second->observations.size())<<'\t'<<d.removed_cycles<<'\t'<<path_text(competing)<<'\t'
                <<s.quality.families.size()<<'\t'<<s.quality.variants<<'\t';
            if (s.quality.similarities_computed)
                decisions<<std::setprecision(8)<<s.quality.spacer_similarity<<'\t'<<s.quality.repeat_spacer_similarity;
            else decisions<<"NA\tNA";
            decisions<<'\t'<<path_text(s.sources)<<'\t'<<d.duplicate_strand_cycles<<'\t'
                     <<flags_text(found==accepted.end()?s.quality:found->second.second->quality)<<'\t'<<d.out_of_range_cycles
                     <<'\t'<<d.invalid_length_observations<<'\n';
            std::set<Path> assigned;size_t observation=0;
            for (const auto& row:s.observations) {
                const auto key=cycle_key(row.cycle);assigned.insert(key);++observation;
                if (!d.retained_paths.count(key)) {
                    alternatives<<s.candidate_id<<'\t'<<
                        ((row.cycle.size()<options.cycle_min || row.cycle.size()>options.cycle_max)?"cycle_length_range":
                        ((row.repeat.size()<options.repeat_min || row.repeat.size()>options.repeat_max ||
                          row.spacer.size()<options.spacer_min || row.spacer.size()>options.spacer_max)?"repeat_spacer_length_range":
                        (found==accepted.end()?d.status:"redundant_cycle_alternative")))<<'\t';
                    emit_path(alternatives,row,observation);++review_observations;
                }
            }
            for (const auto& cycle:s.candidate_cycles) if (!d.retained_paths.count(cycle_key(cycle))) {
                review<<s.candidate_id<<'\t'<<d.status<<'\t'<<path_text(cycle)<<'\t'
                      <<(assigned.count(cycle_key(cycle))?"assigned":"unresolved")<<'\n';++known_review_cycles;
            }
        }
        auto stats=array_output(staging/"CRISPR_graph_summary.tsv");const auto& s=filter.stats;
        stats<<"metric\tvalue\nrepeat_support\t"<<std::setprecision(17)<<options.repeat_support<<"\nk\t"<<options.k
             <<"\narray_filter_mode\t"<<(options.quality.strict?"strict":"conservative")
             <<"\ncycle_min_length\t"<<options.cycle_min<<"\ncycle_max_length\t"<<options.cycle_max
             <<"\nrepeat_min_length\t"<<options.repeat_min<<"\nrepeat_max_length\t"<<options.repeat_max
             <<"\nspacer_min_length\t"<<options.spacer_min<<"\nspacer_max_length\t"<<options.spacer_max
             <<"\nspacer_variant_similarity\t"<<options.quality.variant_similarity<<"\nmax_spacer_similarity\t"<<options.quality.max_spacer_similarity
             <<"\nmax_repeat_spacer_similarity\t"<<options.quality.max_repeat_spacer_similarity<<"\nmax_cycles_exclusive\t"<<options.quality.max_cycles
             <<"\ninput_cycles\t"<<s.input_cycles<<"\ninvalid_cycles\t"<<s.invalid_cycles<<"\nduplicate_cycles_within_groups\t"<<s.duplicate_cycles
             <<"\ninput_groups\t"<<s.input_groups<<"\ngroups_without_repeat\t"<<s.groups_without_repeat
             <<"\nduplicate_systems_removed\t"<<s.duplicate_systems<<"\noverlapping_systems_merged\t"<<s.merged_systems
             <<"\ncurated_candidates\t"<<filter.decisions.size()<<"\nquality_rejected\t"<<s.quality_rejected<<"\noverlap_rejected\t"<<s.overlap_rejected
             <<"\nshared_cycles_removed\t"<<s.shared_cycles_removed<<"\nduplicate_strand_cycles_removed\t"<<s.duplicate_strand_cycles
             <<"\ncycles_outside_length_range\t"<<s.out_of_range_cycles+s.curation_out_of_range_cycles
             <<"\nobservations_outside_repeat_spacer_lengths\t"<<s.invalid_length_observations
             <<"\nredundant_candidates_removed\t"<<s.redundant_candidates
             <<"\nredundant_observations_removed\t"<<s.redundant_observations
             <<"\nsequence_subset_checks\t"<<s.sequence_subset_checks
             <<"\nshared_circular_sequences\t"<<shared_count<<"\nunique_circular_sequences\t"<<cycle_arrays.size()
             <<"\nunique_spacer_sequences_across_arrays\t"<<distinct_spacers.size()
             <<"\narrays\t"<<systems.size()<<"\nspacers\t"<<variant_count
             <<"\nspacer_families\t"<<families<<"\nsequence_variants\t"<<variant_count<<"\nassigned_observations\t"<<observations
             <<"\nreview_observations\t"<<review_observations<<"\nknown_review_cycles\t"<<known_review_cycles<<"\nunresolved_candidate_cycles\t"<<unresolved<<'\n';
        paths.close();summary.close();candidates.close();variants.close();arrays.close();decisions.close();alternatives.close();review.close();shared.close();stats.close();
        // Finish all serialization before replacing processor-owned files.
        for (const auto& item:fs::directory_iterator(destination))
            if (item.is_regular_file() && owned_array_file(item.path().filename().string())) fs::remove(item.path());
        for (const auto& item:fs::directory_iterator(staging)) fs::rename(item.path(),destination/item.path().filename());
        fs::remove(staging);
        std::cout<<"  ▸ Graph filtering: "<<systems.size()<<" candidates, "<<families<<" spacer families, "<<variant_count<<" sequence variants\n"
                 <<"  ▸ Review: "<<s.quality_rejected<<" quality rejections, "<<s.overlap_rejected<<" overlap rejections, "
                 <<s.redundant_candidates<<" redundant candidates removed; "<<review_observations<<" alternative observations\n";
    } catch (...) {std::error_code ec;fs::remove_all(staging,ec);throw;}
}
} // namespace mcaat::arrays
#endif
