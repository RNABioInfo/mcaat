// Included inside extension_module.cpp's namespace after its shared IO helpers.
struct SuppliedTarget {
    std::string id,sequence,error;
    std::vector<TargetContext> contexts;
    std::set<std::string> source_cycles;
    std::vector<std::tuple<std::string,std::string,size_t,size_t,char>> sources;
    TargetExtension extension;
};
static int run_retained_targets(ExtensionGraph& graph,const std::vector<Query>& queries,
                               const fs::path& input,const fs::path& arrays,const fs::path& graph_path,
                               const fs::path& destination,size_t threads,size_t max_bases,const Options& options) {
    const std::string marker="mcaat supplied protospacers v1\n";
    if(fs::exists(destination)) {
        if(!fs::is_directory(destination))throw std::invalid_argument("--output must be a directory");
        if(!fs::is_empty(destination)) {
            std::ifstream previous(destination/"retained_protospacer_mode.txt");
            std::string value((std::istreambuf_iterator<char>(previous)),{});
            if(value!=marker)throw std::invalid_argument("use a new output directory for --protospacers; do not overwrite a discovery run");
        }
    }
    if(fs::weakly_canonical(input)==fs::weakly_canonical(destination/"Protospacers.fasta"))
        throw std::invalid_argument("output must not overwrite the supplied protospacer FASTA");
    std::vector<SuppliedTarget> targets;std::set<std::string> ids;size_t input_bases=0;
    klibpp::SeqStreamIn in(input.string().c_str());klibpp::KSeq record;
    while(in>>record) {
        auto sequence=upper(record.seq);
        if(record.name.empty() || record.name.find_first_of("\t\r\n ")!=std::string::npos || !is_dna(sequence))
            throw std::invalid_argument("protospacer FASTA needs unique IDs and nonempty ACGT sequences");
        if(!ids.insert(record.name).second)throw std::invalid_argument("duplicate protospacer ID: "+record.name);
        input_bases+=sequence.size();if(input_bases>max_bases)throw std::runtime_error("protospacer input exceeds --max-retained-bases");
        SuppliedTarget target;target.id=record.name;target.sequence=std::move(sequence);targets.push_back(std::move(target));
    }
    if(in.err() || in.tqs() || targets.empty())throw std::invalid_argument("invalid or empty protospacer FASTA");
    // Match to full corrected source spacers solely for repeat-context exclusion.
    // Neither source lookup nor subsequent extension changes the supplied targets.
    for(auto& target:targets) {
        std::set<std::tuple<std::string,std::string,size_t,size_t>> unique;
        std::set<std::tuple<std::string,std::string,size_t,size_t,char>> sources;
        for(const auto& q:queries)for(bool flip:{false,true}) {
            auto full=flip?reverse_complement(q.sequence):q.sequence;
            for(size_t pos=full.find(target.sequence);pos!=std::string::npos;pos=full.find(target.sequence,pos+1)) {
                for(const auto& cycle:q.source_cycles)
                    target.source_cycles.insert(flip?reverse_complement(cycle):cycle);
                for(const auto& s:q.sources) {
                    bool source_flip=s.reverse!=flip;
                    auto repeat=source_flip?reverse_complement(s.repeat):s.repeat;
                    if(unique.emplace(full,repeat,pos,pos+target.sequence.size()).second)
                        target.contexts.push_back({full,repeat,pos,pos+target.sequence.size()});
                    size_t begin=source_flip?full.size()-pos-target.sequence.size():pos;
                    sources.emplace(s.id,s.array,begin+s.trim_left,begin+target.sequence.size()+s.trim_left,source_flip?'-':'+');
                    if(unique.size()>1024)break;
                }
                if(unique.size()>1024)break;
            }
            if(unique.size()>1024)break;
        }
        target.sources.assign(sources.begin(),sources.end());
    }
    omp_set_num_threads(threads);std::atomic<bool> cancelled{false};
    std::atomic<size_t> retained{input_bases},completed{0};size_t announced=0;auto last=std::chrono::steady_clock::now();
    #pragma omp parallel for schedule(dynamic)
    for(size_t i=0;i<targets.size();++i) {
        auto& target=targets[i];
        try {
            if(cancelled)continue;
            target.extension=extend_target_with_fallback(graph,target.sequence,target.contexts,options,&target.source_cycles);
            if(target.extension.has_extension && retained.fetch_add(target.extension.selected.sequence.size())+target.extension.selected.sequence.size()>max_bases)
                throw std::runtime_error("retained targets/extensions exceed --max-retained-bases");
        }catch(const std::exception& e){target.error=e.what();cancelled=true;}
        ++completed;
        #pragma omp critical(proto_retained_progress)
        {
            auto now=std::chrono::steady_clock::now();
            size_t done=completed.load();
            if(done!=announced && (done==targets.size() || now-last>std::chrono::seconds(5))) {
                std::cout<<"Processed "<<done<<'/'<<targets.size()<<" supplied protospacers\n"<<std::flush;last=now;announced=done;
            }
        }
    }
    for(const auto& t:targets)if(!t.error.empty())throw std::runtime_error(t.id+": "+t.error);
    fs::create_directories(destination);
    auto staging=destination/(".protospacer-staging-"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    if(!fs::create_directory(staging))throw std::runtime_error("cannot create extension staging directory");
    try {
        auto fasta=output(staging/"Protospacers.fasta"),contigs=output(staging/"ContiguousPaths.fasta");
        auto candidates=output(staging/"protospacer_candidates.tsv"),extensions=output(staging/"protospacer_extensions.tsv");
        auto regions=output(staging/"protospacer_regions.tsv"),associations=output(staging/"protospacer_associations.tsv");
        auto sources=output(staging/"protospacer_sources.tsv"),branches=output(staging/"protospacer_branch_checks.tsv");
        auto adapters=output(staging/"protospacer_adapter_stops.tsv");
        auto placements=output(staging/"protospacer_target_placements.tsv");
        auto placement_checks=output(staging/"protospacer_placement_checks.tsv");
        placement_checks<<"seed_id\tcore\tphase\treason\textended\tlength\tscore\tquery_offset\tquery_begin\tquery_end\texpanded\tlimited\tselected\tsequence\n";
        placements<<"seed_id\tplacement_mode\tobserved_target\tquery_begin\tquery_end\tidentity\tquery_coverage\tmismatches\tplacement_candidates\tplacement_work\tplacement_limited\n";
        candidates<<"seed_id\tsequence\tlength\tretention_status\tsource_observations\textension_status\n";
        extensions<<"seed_id\tretention_status\textension_status\treason\tregion_id\tlength\tleft_bp\tright_bp\tleft_stop\tright_stop\tcompared_paths\tscore\trunner_score\tambiguous\tsearch_limited\texpanded_states\tobserved_stops\n";
        regions<<"region_id\tlength\tsequence\tvertices\tedges\tviral_status\n";
        associations<<"seed_id\tregion_id\tseed_begin\tseed_end\tstrand\textension_status\tobserved_target\tquery_begin\tquery_end\tidentity\tquery_coverage\n";
        sources<<"seed_id\tsource_id\tarray_id\tquery_begin\tquery_end\tstrand\n";
        branches<<"seed_id\tside\tvertex\tedge\ttarget\tprobe_steps\tobservations\tscore\tstatus\tchosen\tdecision\n";
        adapters<<"seed_id\tside\tmotif\tclean_flank_bp\toverlaps_seed\tclean_endpoint_vertex\tdetected_vertex\n";
        std::map<std::string,std::string> written;size_t extended=0,ambiguous=0,limited=0,branch_count=0,adapter_count=0;
        for(const auto& t:targets) {
            const auto& e=t.extension;const auto& c=e.selected;std::string region="-";
            fasta<<'>'<<t.id<<'\n'<<t.sequence<<'\n';
            candidates<<t.id<<'\t'<<t.sequence<<'\t'<<t.sequence.size()<<"\tretained_input\t"<<t.sources.size()<<'\t'<<e.status<<'\n';
            for(const auto& [id,array,begin,end,strand]:t.sources)sources<<t.id<<'\t'<<id<<'\t'<<array<<'\t'<<begin<<'\t'<<end<<'\t'<<strand<<'\n';
            branch_count+=e.diagnostics.branch_checks;adapter_count+=e.diagnostics.adapter_stops;
            for(const auto& b:e.diagnostics.branches)branches<<t.id<<'\t'<<(b.backward?"left":"right")<<'\t'<<b.vertex<<'\t'<<b.edge<<'\t'<<b.target<<'\t'<<b.steps<<'\t'<<b.observations<<'\t'<<b.score<<'\t'<<b.status<<'\t'<<b.chosen<<'\t'<<b.decision<<'\n';
            for(const auto& a:e.diagnostics.adapters)adapters<<t.id<<'\t'<<(a.backward?"left":"right")<<'\t'<<a.motif<<'\t'<<a.clean_flank_bp<<'\t'<<a.overlaps_seed<<'\t'<<a.clean_endpoint<<'\t'<<a.detected_node<<'\n';
            if(e.has_extension) {
                size_t target_begin=c.seed_begin+e.placement.query_offset,target_end=target_begin+t.sequence.size();
                auto observed=c.sequence.substr(target_begin,t.sequence.size());
                if(observed!=e.placement.sequence.substr(e.placement.query_offset,t.sequence.size()))throw std::runtime_error("reported target placement disagrees with graph path");
                ++extended;auto key=canonical(c.sequence);region=identifier(key,"region_");
                auto [it,added]=written.emplace(region,c.sequence);
                if(canonical(it->second)!=key)throw std::runtime_error("region hash collision");
                if(added) {
                    Path edges;for(size_t i=1;i<c.path.size();++i)edges.push_back(graph.edge(c.path[i-1],c.path[i]));
                    contigs<<'>'<<region<<'\n'<<c.sequence<<'\n';
                    regions<<region<<'\t'<<c.sequence.size()<<'\t'<<c.sequence<<'\t'<<path_text(c.path)<<'\t'<<path_text(edges)<<"\tnot_assessed\n";
                }
                bool flip=it->second!=c.sequence;
                associations<<t.id<<'\t'<<region<<'\t'<<(flip?c.sequence.size()-target_end:target_begin)<<'\t'<<(flip?c.sequence.size()-target_begin:target_end)<<'\t'<<(flip?'-':'+')<<'\t'<<e.status<<'\t'<<observed<<'\t'<<e.placement.query_begin<<'\t'<<e.placement.query_end<<'\t'<<e.placement.identity<<'\t'<<e.placement.coverage<<'\n';
                placements<<t.id<<'\t'<<(e.variant_placement?"graph_variant":"exact")<<'\t'<<observed<<'\t'<<e.placement.query_begin<<'\t'<<e.placement.query_end<<'\t'<<e.placement.identity<<'\t'<<e.placement.coverage<<'\t'<<e.placement.mismatches;
            }else {
                placements<<t.id<<"\tunplaced\t-\t0\t0\t0\t0\t0";
            }
            placements<<'\t'<<e.placement_candidates<<'\t'<<e.placement_work<<'\t'<<e.placement_limited<<'\n';
            for(const auto& check:e.placement_checks)placement_checks<<t.id<<'\t'<<check.core<<'\t'<<check.phase<<'\t'<<check.reason<<'\t'<<check.extended<<'\t'<<check.sequence.size()<<'\t'<<check.score<<'\t'<<check.query_offset<<'\t'<<check.query_begin<<'\t'<<check.query_end<<'\t'<<check.expanded<<'\t'<<check.limited<<'\t'<<(e.has_extension && check.sequence==c.sequence)<<'\t'<<check.sequence<<'\n';
            ambiguous+=e.ambiguous;limited+=e.search_limited;
            extensions<<t.id<<"\tretained_input\t"<<e.status<<'\t'<<e.reason<<'\t'<<region<<'\t'
                <<(e.has_extension?c.sequence.size():0)<<'\t'<<(e.has_extension?c.left_length+e.placement.query_offset:0)<<'\t'<<(e.has_extension?c.right_length+c.seed_end-c.seed_begin-e.placement.query_offset-t.sequence.size():0)<<'\t'
                <<(e.has_extension?c.left.stop:"-")<<'\t'<<(e.has_extension?c.right.stop:"-")<<'\t'<<e.compared<<'\t';
            if(e.has_extension)extensions<<c.score;else extensions<<"NA";
            extensions<<'\t';if(e.has_runner)extensions<<e.runner_score;else extensions<<"NA";
            extensions<<'\t'<<e.ambiguous<<'\t'<<e.search_limited<<'\t'<<e.expanded<<'\t';
            bool first=true;for(const auto& stop:e.stops)if(!stop.empty()){if(!first)extensions<<',';extensions<<stop;first=false;}
            extensions<<'\n';
        }
        auto params=output(staging/"protospacer_parameters.json");
        params<<"{\n  \"schema\": 12,\n  \"extension_algorithm\": \"equal_depth_beam_local_window_v4\",\n  \"source_cycle_scope\": \"all_observations_from_source_arrays\",\n  \"stability_observations\": \"all_external_edges\",\n  \"coverage_local_window\": \"all_external_edges\",\n  \"length_tie_policy\": \"fixed_margin_and_within_flank_stability\",\n  \"weak_edge_fraction_weight\": 1,\n  \"mode\": \"extend_supplied_protospacers\",\n  \"protospacers\": "<<json_string(fs::absolute(input).string())
            <<",\n  \"arrays\": "<<json_string(fs::absolute(arrays).string())<<",\n  \"graph\": "<<json_string(fs::absolute(graph_path).string())
            <<",\n  \"graph_k\": "<<graph.k()
            <<",\n  \"target_variant_fallback\": "<<(options.target_variant_fallback?"true":"false")
            <<",\n  \"max_target_placement_states\": "<<options.max_target_placement_states
            <<",\n  \"target_anchor_policy\": \"one_substitution_then_masked_kmer_if_unextended\",\n  \"masked_rescue_stability\": \"zero_truncated_poisson_noise_adjusted\",\n  \"variant_min_identity\": 0.90,\n  \"variant_min_query_coverage\": 0.85"
            <<",\n  \"input_records\": "<<targets.size()<<",\n  \"retained_records\": "<<targets.size()<<",\n  \"extended_records\": "<<extended
            <<",\n  \"not_extended_records\": "<<targets.size()-extended<<",\n  \"ambiguous_records\": "<<ambiguous<<",\n  \"regions\": "<<written.size()
            <<",\n  \"search_limited_records\": "<<limited<<",\n  \"threads\": "<<threads<<",\n  \"max_length\": "<<options.max_bp<<",\n  \"max_flank_length\": "<<options.max_flank_bp
            <<",\n  \"beam_width\": "<<options.beam_width<<",\n  \"max_node_visits\": "<<options.max_node_visits<<",\n  \"max_states_per_side\": "<<options.max_states<<",\n  \"max_query_states\": "<<options.max_query_states
            <<",\n  \"branch_guard\": "<<(options.branch_guard?"true":"false")<<",\n  \"resolve_graph_ties\": "<<(options.resolve_graph_ties?"true":"false")
            <<",\n  \"junction_check\": "<<(options.junction_check?"true":"false")<<",\n  \"joint_flank_check\": "<<(options.joint_flank_check?"true":"false")
            <<",\n  \"coverage_model\": "<<json_string(options.coverage_model)<<",\n  \"coverage_window\": "<<options.coverage_window
            <<",\n  \"max_coverage_penalty\": "<<options.max_coverage_penalty<<",\n  \"array_repeat_identity\": "<<options.array_repeat_identity
            <<",\n  \"min_external_flank_bp\": "<<(options.min_external_flank?options.min_external_flank:graph.k()+1)
            <<",\n  \"branch_lookahead\": "<<options.branch_lookahead<<",\n  \"branch_score_margin\": "<<options.branch_margin<<",\n  \"prediction_score_margin\": "<<options.prediction_margin
            <<",\n  \"correct_spacer_boundaries\": "<<(options.boundaries.enabled?"true":"false")
            <<",\n  \"boundary_support\": "<<options.boundaries.support<<",\n  \"boundary_min_spacers\": "<<options.boundaries.min_spacers<<",\n  \"boundary_min_bases\": "<<options.boundaries.min_bases
            <<",\n  \"max_retained_bases\": "<<max_bases<<",\n  \"adapter_patterns\": "<<adapter_json(options)
            <<",\n  \"branch_checks\": "<<branch_count<<",\n  \"adapter_encounters\": "<<adapter_count<<"\n}\n";
        auto mode=output(staging/"retained_protospacer_mode.txt");mode<<marker;
        fasta.close();contigs.close();candidates.close();extensions.close();regions.close();associations.close();sources.close();branches.close();adapters.close();placements.close();placement_checks.close();params.close();mode.close();
        for(const auto& file:fs::directory_iterator(staging))fs::rename(file.path(),destination/file.path().filename());
        fs::remove(staging);
        std::cout<<"Retained "<<targets.size()<<'/'<<targets.size()<<" supplied protospacers unchanged; "<<extended<<" with an extension, "<<targets.size()-extended<<" without an extension\n";
    }catch(...){std::error_code error;fs::remove_all(staging,error);throw;}
    return 0;
}
