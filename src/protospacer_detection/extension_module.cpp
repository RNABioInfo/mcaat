#include "protospacer_detection/sdbg_extension_graph.h"
#include "protospacer_detection/extension_evidence.h"
#include <kseq++/seqio.hpp>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <chrono>
#include <atomic>
#include <omp.h>

namespace mcaat::proto {
namespace fs=std::filesystem;
using Fields=std::map<std::string,std::string>;
static std::vector<std::string> split(const std::string& line) {
    std::vector<std::string> out;size_t begin=0;
    while(true){auto end=line.find('\t',begin);out.push_back(line.substr(begin,end-begin));if(end==std::string::npos)break;begin=end+1;}return out;
}
static void table(const fs::path& path,const std::function<void(const Fields&)>& callback) {
    std::ifstream in(path);if(!in)throw std::runtime_error("cannot open "+path.string());
    std::string line;if(!std::getline(in,line))throw std::runtime_error("empty table: "+path.string());
    if(!line.empty() && line.back()=='\r')line.pop_back();auto header=split(line);
    if(std::set<std::string>(header.begin(),header.end()).size()!=header.size())throw std::runtime_error("duplicate TSV columns");
    while(std::getline(in,line)) {
        if(!line.empty() && line.back()=='\r')line.pop_back();if(line.empty())continue;
        auto values=split(line);if(values.size()!=header.size())throw std::runtime_error("wrong column count in "+path.string());
        Fields row;for(size_t i=0;i<header.size();++i)row[header[i]]=values[i];callback(row);
    }
    if(in.bad())throw std::runtime_error("read failed: "+path.string());
}
static std::string identifier(const std::string& value,const char* prefix) {
    uint64_t hash=14695981039346656037ULL;for(unsigned char c:value){hash^=c;hash*=1099511628211ULL;}
    std::ostringstream s;s<<prefix<<std::hex<<std::setw(16)<<std::setfill('0')<<hash;return s.str();
}
static std::string path_text(const Path& p) {std::ostringstream out;for(size_t i=0;i<p.size();++i){if(i)out<<',';out<<p[i];}return out.str().empty()?"-":out.str();}
struct Source {std::string id,array,observation,repeat,spacer;bool reverse=false;};
struct Entry {Seed seed;std::string id;std::vector<Candidate> paths;};
struct Decision {
    std::string seed,outcome,reason,left_stop,right_stop,replacement;
    size_t begin=0,end=0,expanded=0;
};
struct Query {
    std::string sequence,status="unmapped",error;
    std::vector<Source> sources;
    std::set<Node> array_edges;
    std::vector<Entry> entries;
    std::vector<Decision> decisions;
    std::set<std::string> viable_seeds;
    size_t compared_paths=0,compared_regions=0;
    bool has_runner=false,comparable_runner=false,score_tie=false;
    double runner_score=0;
    size_t eligible=0;bool truncated=false;
    SearchDiagnostics adapter_diagnostics;
    std::vector<std::string> adapter_seeds;
};
static void record(Query& q,const Entry& e,const Candidate& c,const std::string& outcome,
                   const std::string& reason,const std::string& replacement="-") {
    q.decisions.push_back({e.id,outcome,reason,c.left.stop,c.right.stop,replacement,e.seed.begin,e.seed.end,c.expanded});
}
static void trim_left(Candidate& c) {
    size_t remove=c.seed_begin;c.path.erase(c.path.begin(),c.path.begin()+remove);
    c.sequence.erase(0,remove);c.seed_end-=remove;c.seed_begin=0;c.left_length=0;
    c.left.coverage={};c.external_left=0;
}
static void trim_right(Candidate& c,size_t k) {
    c.sequence.resize(c.seed_end);c.path.resize(c.seed_end-k+1);c.right_length=0;
    c.right.coverage={};c.external_right=0;
}
static void compact_intervals(Query& q) {
    // Equivalent complete paths at the same query alignment can arise from
    // nested seed intervals. Retain maximal intervals for that exact context;
    // a different extension sequence or alignment remains a separate candidate.
    std::sort(q.entries.begin(),q.entries.end(),[](const auto& a,const auto& b) {
        auto al=a.seed.end-a.seed.begin,bl=b.seed.end-b.seed.begin;
        return al!=bl?al>bl:a.id<b.id;
    });
    RetainedIntervalContexts contexts;
    for(auto& e:q.entries) {
        std::vector<Candidate> keep;
        for(auto& c:e.paths) {
            auto replacement=contexts.replacement(c.sequence,c.seed_begin,e.seed.begin,e.seed.end,e.id);
            if(!replacement.empty())record(q,e,c,"redundant","nested_interval_same_path",replacement);
            else keep.push_back(std::move(c));
        }
        e.paths=std::move(keep);
    }
    q.entries.erase(std::remove_if(q.entries.begin(),q.entries.end(),[](const auto& e){return e.paths.empty();}),q.entries.end());
}
static void select_prediction(Query& q,const Options& options) {
    if(q.entries.empty())return;
    size_t best_entry=0,best_path=0;std::set<std::string> regions;
    for(size_t i=0;i<q.entries.size();++i) {
        const auto& e=q.entries[i];q.viable_seeds.insert(e.id);
        for(size_t j=0;j<e.paths.size();++j) {
            const auto& c=e.paths[j];++q.compared_paths;regions.insert(canonical(c.sequence));
            const auto& best=q.entries[best_entry];
            if(better_prediction(c,e.seed.end-e.seed.begin,best.paths[best_path],best.seed.end-best.seed.begin,options)) {best_entry=i;best_path=j;}
        }
    }
    q.compared_regions=regions.size();
    const auto& winner=q.entries[best_entry].paths[best_path];auto winner_key=canonical(winner.sequence);
    const Candidate* runner=nullptr;size_t runner_matched=0;
    for(const auto& e:q.entries)for(const auto& c:e.paths) {
        if(canonical(c.sequence)==winner_key)continue;
        auto matched=e.seed.end-e.seed.begin;
        if(!runner || better_prediction(c,matched,*runner,runner_matched,options)){runner=&c;runner_matched=matched;}
    }
    if(runner) {
        q.has_runner=true;q.runner_score=runner->score;
        q.comparable_runner=options.pair_weight==0 || bool(winner.pair_support)==bool(runner->pair_support);
        q.score_tie=q.comparable_runner && std::abs(winner.score-runner->score)<1e-12;
    }
    Entry selected=std::move(q.entries[best_entry]);Candidate chosen=std::move(selected.paths[best_path]);
    selected.paths.clear();selected.paths.push_back(std::move(chosen));
    q.entries.clear();q.entries.push_back(std::move(selected));
}
static std::vector<Query> load_queries(const fs::path& folder) {
    std::map<std::string,Query> queries;std::set<std::string> ids;
    auto add=[&](Source source) {
        source.repeat=upper(source.repeat);source.spacer=upper(source.spacer);
        if(!is_dna(source.repeat) || !is_dna(source.spacer))throw std::runtime_error("array observation contains non-ACGT sequence: "+source.id);
        if(!ids.insert(source.id).second)throw std::runtime_error("duplicate source observation ID: "+source.id);
        auto key=canonical(source.spacer);source.reverse=key!=source.spacer;
        auto& q=queries[key];q.sequence=key;q.sources.push_back(std::move(source));
    };
    if(fs::exists(folder/"CRISPR_graph_paths.tsv")) {
        table(folder/"CRISPR_graph_paths.tsv",[&](const Fields& row) {
            Source s;s.array=row.at("array_id");s.observation=row.at("observation_id");s.repeat=row.at("repeat");s.spacer=row.at("spacer");
            s.id=s.array+":"+s.observation;add(std::move(s));
        });
    } else {
        std::vector<fs::path> files;
        for(auto& f:fs::directory_iterator(folder)) {
            auto name=f.path().filename().string();
            if(f.is_regular_file() && name.rfind("CRISPR_Arrays_",0)==0 && f.path().extension()==".txt")files.push_back(f.path());
        }
        std::sort(files.begin(),files.end());
        for(const auto& file:files) {
            std::ifstream in(file);if(!in)throw std::runtime_error("cannot read array report");
            std::string line,array,repeat;size_t observation=0;
            while(std::getline(in,line)) {
                if(!line.empty() && line.back()=='\r')line.pop_back();if(line.empty() || line[0]=='#')continue;
                if(line[0]=='>'){std::istringstream h(line.substr(1));h>>array;repeat.clear();observation=0;continue;}
                if(array.empty())continue;
                if(line[0]!=' ' && line[0]!='\t'){std::istringstream r(line);r>>repeat;continue;}
                std::string variant,spacer,extra;std::istringstream fields(line);
                if(!(fields>>variant>>spacer) || (fields>>extra))throw std::runtime_error("malformed array observation in "+file.string());
                // MCAAT's text format writes an all-dash repeat when the entry
                // is identical to the consensus printed above it.
                if(variant.find_first_not_of('-')==std::string::npos)variant=repeat;
                Source s{file.filename().string()+":"+array+":"+std::to_string(++observation),array,std::to_string(observation),variant,spacer};add(std::move(s));
            }
            if(in.bad())throw std::runtime_error("array report read failed");
        }
    }
    if(queries.empty())throw std::runtime_error("no spacer observations found in --arrays");
    std::vector<Query> result;for(auto& [seq,q]:queries)result.push_back(std::move(q));return result;
}
static std::vector<PamRule> load_pams(const fs::path& path) {
    std::vector<PamRule> rules;if(path.empty())return rules;
    table(path,[&](const Fields& row) {
        PamRule r{row.at("array_id"),row.at("subtype"),row.at("strand"),row.at("side"),upper(row.at("motif"))};
        if(r.array.empty() || r.subtype.empty() || r.motif.empty() ||
           (r.strand!="+" && r.strand!="-" && r.strand!="?") || (r.side!="5prime" && r.side!="3prime" && r.side!="none"))throw std::runtime_error("invalid PAM rule");
        if(r.side=="none"){if(r.motif!="-")throw std::runtime_error("PAM side none requires motif -");}
        else for(char c:r.motif)iupac(c,'A');rules.push_back(std::move(r));
    });return rules;
}
static std::pair<std::string,std::string> pam_status(const Query& q,const Entry& e,const Candidate& c,const Source& source,const std::vector<PamRule>& rules) {
    std::string status="unknown_no_model",subtype="-";
    for(const auto& r:rules)if(r.array==source.array) {
        auto value=check_pam(c.sequence,c.seed_begin,c.seed_end,source.reverse,r,e.seed.begin==0 && e.seed.end==q.sequence.size());
        if(status=="match")continue;
        if(value=="match" || status=="unknown_no_model" || value=="mismatch"){status=value;subtype=r.subtype;}
    }
    return {status,subtype};
}
static std::vector<std::string> junctions(ExtensionGraph& graph,const Candidate& c) {
    std::set<std::string> result;
    for(size_t i=1;i+1<c.path.size();++i)
        if(graph.neighbors(c.path[i],false).size()>1 || graph.neighbors(c.path[i],true).size()>1)
            result.insert(canonical(c.sequence.substr(i-1,graph.k()+2)));
    return {result.begin(),result.end()};
}
static std::string mate_name(std::string s) {
    if(s.size()>2 && s[s.size()-2]=='/' && (s.back()=='1' || s.back()=='2'))s.resize(s.size()-2);return s;
}
static void read_evidence(std::vector<Query>& queries,ExtensionGraph& graph,const fs::path& first,const fs::path& second,
                          size_t min_insert,size_t max_insert,size_t max_index,size_t max_read_work,const Options& options) {
    std::map<std::string,size_t> ids;std::vector<ReadRegion> regions;
    for(auto& q:queries)for(auto& e:q.entries)for(auto& c:e.paths) {
        auto key=canonical(c.sequence);auto [it,inserted]=ids.emplace(key,regions.size());
        if(inserted)regions.push_back({c.sequence,junctions(graph,c)});
    }
    ReadEvidence evidence(graph.k(),max_index,max_read_work);evidence.index(std::move(regions));
    klibpp::SeqStreamIn r1(first.string().c_str());
    std::unique_ptr<klibpp::SeqStreamIn> r2;if(!second.empty())r2=std::make_unique<klibpp::SeqStreamIn>(second.string().c_str());
    klibpp::KSeq a,b;
    auto last_progress=std::chrono::steady_clock::now();
    while(r1>>a) {
        auto s1=upper(a.seq);evidence.read(s1);
        if(r2) {
            if(!(*r2>>b))throw std::runtime_error("R2 ends before R1 or contains malformed FASTQ");
            if(mate_name(a.name)!=mate_name(b.name))throw std::runtime_error("read-pair names disagree: "+a.name+" / "+b.name);
            auto s2=upper(b.seq);evidence.read(s2);evidence.pair(s1,s2,min_insert,max_insert);
        }
        if(evidence.reads%100000==0) {
            auto now=std::chrono::steady_clock::now();
            if(now-last_progress>std::chrono::seconds(5)) {
                std::cout<<"Scanned "<<evidence.reads<<" reads for candidate support\n"<<std::flush;last_progress=now;
            }
        }
    }
    if(r1.err() || r1.tqs())throw std::runtime_error("invalid/truncated R1 FASTA/FASTQ");
    if(r2 && ((*r2>>b) || r2->err() || r2->tqs()))throw std::runtime_error("R2 has extra records or malformed FASTQ");
    for(auto& q:queries)for(auto& e:q.entries)for(auto& c:e.paths) {
        auto id=ids.at(canonical(c.sequence));const auto& r=evidence.regions[id];
        c.junctions=r.junctions.size();
        for(const auto& word:r.junctions){size_t count=evidence.junction_count(word);c.read_support+=count;c.supported_junctions+=count>0;}
        c.pair_support=r.pairs;c.unique_pairs=r.unique_pairs;
        double support=c.junctions?double(c.supported_junctions)/c.junctions:0;
        // Overlapping/nested regions can support the same real fragment. Do not
        // discard that evidence just because it is not unique to one region.
        c.score+=options.read_weight*support+options.pair_weight*std::log1p(c.pair_support);
    }
    std::cout<<"Read evidence: "<<evidence.reads<<" reads, "<<evidence.pairs<<" pairs; exact support over retained alternatives\n";
}
static std::ofstream output(const fs::path& path) {
    std::ofstream out(path);out.exceptions(std::ios::badbit|std::ios::failbit);out<<std::setprecision(17);return out;
}
static std::string json_string(const std::string& value) {
    std::string out="\"";for(unsigned char c:value){if(c=='"' || c=='\\'){out+='\\';out+=c;}else if(c<32){std::ostringstream hex;hex<<"\\u"<<std::hex<<std::setw(4)<<std::setfill('0')<<int(c);out+=hex.str();}else out+=c;}return out+'"';
}
static size_t integer(const std::string& flag,const std::string& value,bool zero=false) {
    if(value.empty() || value.find_first_not_of("0123456789")!=std::string::npos)throw std::invalid_argument(flag+" needs an integer");
    size_t n=std::stoull(value);if((!zero && !n) || n>1000000000)throw std::invalid_argument(flag+" out of range");return n;
}
static double number(const std::string& flag,const std::string& value) {
    size_t used=0;double n=std::stod(value,&used);if(used!=value.size() || !std::isfinite(n) || n<0)throw std::invalid_argument(flag+" needs a finite nonnegative number");return n;
}

int run_extension_module(int argc,char** argv) {
    Options options;fs::path arrays,graph_path,destination,reads1,reads2,pam_file;
    size_t threads=1,min_insert=0,max_insert=1000,max_index=5000000,max_read_work=1000000,max_retained_bases=50000000;
    for(int i=1;i<argc;++i) {
        std::string flag=argv[i];
        if(flag=="--help" || flag=="-h") {
            std::cout<<"mcaat phage-curate --arrays DIR --graph PREFIX [--output DIR]\n"
                "Exact protospacer seeds and bounded bidirectional extension; partial results retained.\n"
                "Primary outputs require an external flank. Unresolved seeds are reported separately.\n"
                "Select ONE prediction per distinct spacer across all seed intervals, ranked by multiplicity stability.\n"
                "--threads 1 --beam-width 32 --max-length 2500 --max-states 50000\n"
                "--max-paths 16 (internal shortlist per seed; final output is one per spacer)\n"
                "--max-seeds 128 --max-query-states 200000 --min-match-bp 0 (0 = graph k)\n"
                "--max-retained-bases 50000000 (abort explicitly before excessive result storage)\n"
                "--coverage-model auto|poisson|negative-binomial --coverage-window 7\n"
                "--array-repeat-identity 0.85 (edit similarity for source-array repeat copies)\n"
                "Default adapter guard: confirmed Illumina cores and exact Nextera, both strands.\n"
                "--adapter-sequences SEQ[,SEQ...] (replace defaults with exact motifs; none disables)\n"
                "--max-coverage-penalty 20\n"
                "Optional evidence: --reads1 FILE [--reads2 FILE] --min-insert 0 --max-insert 1000\n"
                "--max-read-index 5000000 --max-read-work 1000000 --read-weight 1 --pair-weight 2\n"
                "--pam-rules FILE (TSV: array_id subtype strand side motif) --pam-weight 1\n"
                "PAM strand: +,-,? relative to submitted spacer; side: 5prime, 3prime, or none (motif -).\n"
                "Read evidence reranks graph-generated alternatives; it cannot recover pruned paths.\n";
            return 0;
        }
        if(++i>=argc)throw std::invalid_argument("missing value for "+flag);std::string value=argv[i];
        if(flag=="--arrays" || flag=="-a")arrays=value;
        else if(flag=="--graph" || flag=="-g")graph_path=value;
        else if(flag=="--output" || flag=="--output-folder")destination=value;
        else if(flag=="--threads")threads=integer(flag,value);
        else if(flag=="--beam-width")options.beam_width=integer(flag,value);
        else if(flag=="--max-length")options.max_bp=integer(flag,value);
        else if(flag=="--max-states")options.max_states=integer(flag,value);
        else if(flag=="--max-query-states")options.max_query_states=integer(flag,value);
        else if(flag=="--max-retained-bases")max_retained_bases=integer(flag,value);
        else if(flag=="--max-paths")options.max_paths=integer(flag,value);
        else if(flag=="--max-seeds")options.max_seeds=integer(flag,value);
        else if(flag=="--min-match-bp")options.min_match=integer(flag,value,true);
        else if(flag=="--coverage-window")options.coverage_window=integer(flag,value);
        else if(flag=="--coverage-model")options.coverage_model=value;
        else if(flag=="--array-repeat-identity")options.array_repeat_identity=number(flag,value);
        else if(flag=="--adapter-sequences") {
            options.adapters.clear();
            if(value!="none") {
                size_t begin=0;
                while(true) {
                    auto end=value.find(',',begin);options.adapters.push_back({upper(value.substr(begin,end-begin))});
                    if(end==std::string::npos)break;begin=end+1;
                }
            }
        }
        else if(flag=="--max-coverage-penalty")options.max_coverage_penalty=number(flag,value);
        else if(flag=="--reads1")reads1=value;
        else if(flag=="--reads2")reads2=value;
        else if(flag=="--min-insert")min_insert=integer(flag,value,true);
        else if(flag=="--max-insert")max_insert=integer(flag,value);
        else if(flag=="--max-read-index")max_index=integer(flag,value);
        else if(flag=="--max-read-work")max_read_work=integer(flag,value);
        else if(flag=="--read-weight")options.read_weight=number(flag,value);
        else if(flag=="--pair-weight")options.pair_weight=number(flag,value);
        else if(flag=="--pam-rules")pam_file=value;
        else if(flag=="--pam-weight")options.pam_weight=number(flag,value);
        else throw std::invalid_argument("unknown phage-curate option: "+flag);
    }
    options.validate();
    if(arrays.empty() || graph_path.empty())throw std::invalid_argument("phage-curate requires --arrays and --graph");
    if(!fs::is_directory(arrays))throw std::invalid_argument("--arrays must be a directory");
    if(!reads2.empty() && reads1.empty())throw std::invalid_argument("--reads2 requires --reads1");
    if(min_insert>max_insert || threads>1024)throw std::invalid_argument("invalid insert range or thread count");
    for(const auto& f:{reads1,reads2,pam_file})if(!f.empty() && !fs::is_regular_file(f))throw std::invalid_argument("missing input: "+f.string());
    if(destination.empty())destination=arrays;
    auto queries=load_queries(arrays);auto rules=load_pams(pam_file);
    for(const auto& q:queries)if(q.sequence.size()>options.max_bp)throw std::invalid_argument("--max-length must accommodate every complete input spacer");
    if(fs::is_directory(graph_path))graph_path/="graph";
    if(!fs::is_regular_file(graph_path.string()+".sdbg_info"))throw std::invalid_argument("missing graph metadata: "+graph_path.string());
    SDBG sdbg;std::cout<<"Loading graph...\n";sdbg.LoadFromFile(graph_path.string().c_str());ExtensionGraph graph(sdbg);
    if(options.max_bp<graph.k())throw std::invalid_argument("--max-length is shorter than graph k");
    std::cout<<"Processing "<<queries.size()<<" distinct spacer sequences with "<<threads<<" threads\n";
    omp_set_num_threads(threads);
    std::atomic<size_t> completed{0},retained_bases{0};std::atomic<bool> cancelled{false};auto last_progress=std::chrono::steady_clock::now();
    #pragma omp parallel for schedule(dynamic)
    for(size_t i=0;i<queries.size();++i) {
        auto& q=queries[i];
        try { [&]() {
            if(cancelled)return;
            if(q.sequence.size()<graph.k()){q.status="shorter_than_k";return;}
            if(AdapterGuard(options.adapters).contains(q.sequence)){q.status="adapter_in_spacer";return;}
            Path full=graph.trace(q.sequence);if(full.empty()){q.status="not_fully_mapped";return;}
            std::set<std::pair<std::string,std::string>> contexts;
            std::set<std::string> repeat_set;
            for(const auto& source:q.sources)if(contexts.emplace(source.repeat,source.spacer).second) {
                auto edges=graph.array_context(source.repeat,source.spacer);
                q.array_edges.insert(edges.first.begin(),edges.first.end());
                repeat_set.insert(source.reverse?reverse_complement(source.repeat):source.repeat);
            }
            std::vector<std::string> repeats(repeat_set.begin(),repeat_set.end());
            auto isolated=isolate(graph,full,{},options);q.eligible=isolated.eligible;q.truncated=isolated.truncated;
            q.status=isolated.seeds.empty()?"no_external_flank":"candidate";
            size_t used=0;
            for(auto& seed:isolated.seeds) {
                if(cancelled)return;
                Entry entry;entry.id=identifier(q.sequence+":"+std::to_string(seed.begin)+":"+std::to_string(seed.end),"seed_");
                entry.seed=std::move(seed);
                auto left_context=flank_contexts(q.sequence,entry.seed.begin,entry.seed.end,repeats,true);
                auto right_context=flank_contexts(q.sequence,entry.seed.begin,entry.seed.end,repeats,false);
                if(used+2>options.max_query_states) {
                    Candidate c;c.path=entry.seed.path;c.sequence=graph.spell(c.path);c.seed_end=c.sequence.size();
                    c.left.stop=c.right.stop="query_work_limit";c.search_limited=true;entry.paths.push_back(std::move(c));
                } else {
                    auto local=options;local.max_states=std::min(options.max_states,(options.max_query_states-used)/2);
                    entry.paths=reconstruct(graph,entry.seed,{},local,&q.array_edges,&left_context,&right_context,&q.adapter_diagnostics);
                    q.adapter_seeds.resize(q.adapter_diagnostics.adapters.size(),entry.id);
                    if(!entry.paths.empty())used+=entry.paths.front().expanded;
                }
                // Recheck the EMITTED flanks: the joined length budget may have
                // trimmed away the bases that originally established divergence.
                for(auto& c:entry.paths) {
                    auto left=c.sequence.substr(0,c.seed_begin);std::reverse(left.begin(),left.end());
                    if(!left.empty()) {
                        auto state=context_state(left,left_context,options.array_repeat_identity);
                        if(state!=ContextState::external) {
                            c.left.context=state==ContextState::array?"array":"incomplete";
                            if(state==ContextState::array)c.left.stop="source_array_context";
                            trim_left(c);
                        }
                    }
                    auto right=c.sequence.substr(c.seed_end);
                    if(!right.empty()) {
                        auto state=context_state(right,right_context,options.array_repeat_identity);
                        if(state!=ContextState::external) {
                            c.right.context=state==ContextState::array?"array":"incomplete";
                            if(state==ContextState::array)c.right.stop="source_array_context";
                            trim_right(c,graph.k());
                        }
                    }
                }
                for(auto& c:entry.paths)for(size_t pos=1;pos<c.path.size();++pos) {
                    if(q.array_edges.count(graph.edge(c.path[pos-1],c.path[pos])))continue;
                    if(pos<=c.seed_begin)++c.external_left;
                    if(pos>c.seed_end-graph.k())++c.external_right;
                }
                for(auto& c:entry.paths) {
                    // An array-only prefix/suffix is not a reconstructed target
                    // flank. Keep the seed and independently supported side.
                    if(!c.external_left && c.left_length) {
                        trim_left(c);c.left.context="array";c.left.stop="array_context_only";
                    }
                    if(!c.external_right && c.right_length) {
                        trim_right(c,graph.k());c.right.context="array";c.right.stop="array_context_only";
                    }
                    c.score=stability(c.left.coverage,c.right.coverage).score();
                }
                // A work-limited or failed-join seed is diagnostic information,
                // never a detection in the primary FASTA/association outputs.
                entry.paths.erase(std::remove_if(entry.paths.begin(),entry.paths.end(),[&](const auto& c) {
                    if(c.external_left || c.external_right)return false;
                    bool unresolved=c.search_limited || c.left.context=="incomplete" || c.right.context=="incomplete";
                    std::string reason=(c.left.context=="array" || c.right.context=="array")?"source_array_context":"no_external_flank";
                    if(unresolved)reason=c.left.stop=="query_work_limit"?"query_work_limit":
                        c.left.stop=="incompatible_flanks"?"incompatible_flanks":
                        (c.left.stop=="length_limit" || c.right.stop=="length_limit")?"length_limit":"incomplete_search_or_context";
                    if(c.left.stop=="adapter_in_seed" || c.right.stop=="adapter_in_seed" ||
                       c.left.stop=="sequencing_adapter" || c.right.stop=="sequencing_adapter") {
                        reason="sequencing_adapter";
                        unresolved=c.search_limited;
                    }
                    record(q,entry,c,unresolved?"unresolved":"rejected",reason);return true;
                }),entry.paths.end());
                std::set<Path> unique_paths;
                entry.paths.erase(std::remove_if(entry.paths.begin(),entry.paths.end(),[&](const auto& c){return !unique_paths.insert(c.path).second;}),entry.paths.end());
                size_t bases=0;for(const auto& c:entry.paths)bases+=c.sequence.size();
                if(retained_bases.fetch_add(bases)+bases>max_retained_bases)throw std::runtime_error("results exceed --max-retained-bases; increase the budget or lower --beam-width/--max-seeds");
                if(!entry.paths.empty()){q.viable_seeds.insert(entry.id);q.entries.push_back(std::move(entry));}
            }
            if(!q.entries.empty())q.status="candidate";
            else if(q.truncated || std::any_of(q.decisions.begin(),q.decisions.end(),[](const auto& d){return d.outcome=="unresolved";}))q.status="unresolved";
            else q.status="no_external_flank";
        }(); }catch(const std::exception& e){q.error=e.what();cancelled=true;}
        ++completed;
        #pragma omp critical(proto_progress)
        {
            auto now=std::chrono::steady_clock::now();
            size_t done=completed.load();
            if(done==queries.size() || now-last_progress>std::chrono::seconds(5)) {
                std::cout<<"Extended "<<done<<'/'<<queries.size()<<" queries\n"<<std::flush;last_progress=now;
            }
        }
    }
    for(const auto& q:queries)if(!q.error.empty())throw std::runtime_error("spacer "+identifier(q.sequence,"query_")+": "+q.error);
    if(!reads1.empty())read_evidence(queries,graph,reads1,reads2,min_insert,max_insert,max_index,max_read_work,options);
    for(auto& q:queries)for(auto& e:q.entries) {
        for(auto& c:e.paths) {
            bool match=false;for(const auto& s:q.sources)match|=pam_status(q,e,c,s,rules).first=="match";
            if(match){c.pam_bonus=options.pam_weight;c.score+=c.pam_bonus;}
            c.score-=1.-double(e.seed.end-e.seed.begin)/q.sequence.size();
        }
        std::sort(e.paths.begin(),e.paths.end(),[&](const Candidate& a,const Candidate& b){
            return better_prediction(a,e.seed.end-e.seed.begin,b,e.seed.end-e.seed.begin,options);
        });
        if(e.paths.size()>options.max_paths){e.paths.resize(options.max_paths);for(auto& c:e.paths)c.search_limited=true;}
    }
    // Only suppress a redundant interval once its replacement has survived
    // evidence ranking and output limits.
    for(auto& q:queries)compact_intervals(q);
    for(auto& q:queries)select_prediction(q,options);
    fs::create_directories(destination);
    auto staging=destination/(".protospacer-stage-"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    fs::create_directory(staging);
    try {
        auto query_out=output(staging/"protospacer_queries.tsv");
        query_out<<"source_id\tarray_id\tobservation_id\tquery_id\tspacer\tmapping_status\teligible_intervals\tseed_limit_reached\tretained_intervals\tunresolved_intervals\n";
        auto seed_out=output(staging/"protospacer_candidates.tsv");
        seed_out<<"seed_id\tquery_id\tquery_begin\tquery_end\tmatched_bases\tquery_coverage\tidentity\tsequence\tsource_count\tseed_vertices\n";
        auto seed_fasta=output(staging/"Protospacers.fasta");
        auto matched_fasta=output(staging/"MatchedSpacers.fasta");
        auto unresolved_fasta=output(staging/"UnresolvedSeeds.fasta");
        auto unresolved_out=output(staging/"protospacer_unresolved.tsv");
        unresolved_out<<"seed_id\tquery_id\tquery_begin\tquery_end\tsequence\treasons\n";
        auto decisions=output(staging/"protospacer_decisions.tsv");
        decisions<<"query_id\tseed_id\tquery_begin\tquery_end\toutcome\treason\tleft_stop\tright_stop\texpanded_states\treplacement_seed_id\n";
        auto selection=output(staging/"protospacer_selection.tsv");
        auto adapter_stops=output(staging/"protospacer_adapter_stops.tsv");
        adapter_stops<<"query_id\tseed_id\tside\tmotif\tclean_flank_bp\toverlaps_seed\tclean_endpoint_vertex\tdetected_vertex\n";
        selection<<"query_id\tseed_id\tregion_id\tcompared_paths\tcompared_regions\ttotal_score\tstability_score\twithin_flank_variance\tadjacent_change\tleft_right_difference\tsample_penalty\tleft_samples\tright_samples\tleft_multiplicity\tright_multiplicity\tquery_coverage\tmultiplicity_evidence\trunner_up_score\tscore_gap\tscore_tie\tsearch_limited\n";
        auto regions_out=output(staging/"protospacer_regions.tsv");
        regions_out<<"region_id\tlength\tsequence\tvertices\tedges\tviral_status\n";
        auto fasta=output(staging/"ContiguousPaths.fasta");
        auto associations=output(staging/"protospacer_associations.tsv");
        associations<<"source_id\tarray_id\tseed_id\tregion_id\tquery_begin\tquery_end\tseed_begin\tseed_end\tstrand\tstatus\tleft_stop\tright_stop\tleft_coverage\tright_coverage\tleft_dispersion\tright_dispersion\tcoverage_samples\tcoverage_score\ttotal_score\talternative_rank\tsearch_limited\texpanded_states\tread_evidence\tsupported_junctions\tjunctions\tread_support\tpair_support\tunique_pair_support\tpam_status\tsubtype\texternal_left_edges\texternal_right_edges\tambiguity\tleft_context\tright_context\n";
        std::map<std::string,std::pair<std::string,std::string>> regions;size_t seeds=0,links=0,sources=0;
        size_t matched_queries=0,full_matched_queries=0,unresolved_queries=0,unresolved_seeds=0,rejected_context_seeds=0,redundant_paths=0;
        size_t adapter_branches=0,adapter_rows=0,adapter_queries=0;
        for(const auto& q:queries) {
            auto query_id=identifier(q.sequence,"query_");
            adapter_branches+=q.adapter_diagnostics.adapter_stops;
            adapter_queries+=q.status=="adapter_in_spacer";
            for(size_t j=0;j<q.adapter_diagnostics.adapters.size();++j) {
                const auto& a=q.adapter_diagnostics.adapters[j];++adapter_rows;
                adapter_stops<<query_id<<'\t'<<q.adapter_seeds[j]<<'\t'<<(a.backward?"left":"right")<<'\t'<<a.motif<<'\t'
                    <<a.clean_flank_bp<<'\t'<<a.overlaps_seed<<'\t'<<a.clean_endpoint<<'\t'<<a.detected_node<<'\n';
            }
            if(!q.entries.empty()) {
                ++matched_queries;matched_fasta<<'>'<<query_id<<'\n'<<q.sequence<<'\n';
                full_matched_queries+=std::any_of(q.entries.begin(),q.entries.end(),[&](const auto& e){return e.seed.begin==0 && e.seed.end==q.sequence.size();});
            }
            unresolved_queries+=q.status=="unresolved";
            std::set<std::string> context_rejected;
            std::map<std::string,std::vector<const Decision*>> pending;
            for(const auto& d:q.decisions) {
                decisions<<query_id<<'\t'<<d.seed<<'\t'<<d.begin<<'\t'<<d.end<<'\t'<<d.outcome<<'\t'<<d.reason<<'\t'<<d.left_stop<<'\t'<<d.right_stop<<'\t'<<d.expanded<<'\t'<<d.replacement<<'\n';
                if(d.outcome=="unresolved" && !q.viable_seeds.count(d.seed))pending[d.seed].push_back(&d);
                if(d.reason=="source_array_context")context_rejected.insert(d.seed);
                redundant_paths+=d.outcome=="redundant";
            }
            for(auto it=context_rejected.begin();it!=context_rejected.end();) {
                if(q.viable_seeds.count(*it) || pending.count(*it))it=context_rejected.erase(it);else ++it;
            }
            unresolved_seeds+=pending.size();rejected_context_seeds+=context_rejected.size();
            for(const auto& [seed,items]:pending) {
                const auto& first=*items.front();auto seq=q.sequence.substr(first.begin,first.end-first.begin);
                std::set<std::string> reasons;for(auto item:items)reasons.insert(item->reason);
                std::string joined;for(const auto& reason:reasons){if(!joined.empty())joined+=',';joined+=reason;}
                unresolved_fasta<<'>'<<seed<<'\n'<<seq<<'\n';
                unresolved_out<<seed<<'\t'<<query_id<<'\t'<<first.begin<<'\t'<<first.end<<'\t'<<seq<<'\t'<<joined<<'\n';
            }
            for(const auto& source:q.sources) {
                ++sources;query_out<<source.id<<'\t'<<source.array<<'\t'<<source.observation<<'\t'<<query_id<<'\t'<<source.spacer<<'\t'<<q.status<<'\t'<<q.eligible<<'\t'<<q.truncated<<'\t'<<q.entries.size()<<'\t'<<pending.size()<<'\n';
            }
            for(const auto& e:q.entries) {
                ++seeds;auto seq=q.sequence.substr(e.seed.begin,e.seed.end-e.seed.begin);
                seed_out<<e.id<<'\t'<<identifier(q.sequence,"query_")<<'\t'<<e.seed.begin<<'\t'<<e.seed.end<<'\t'<<seq.size()<<'\t'<<double(seq.size())/q.sequence.size()<<"\t1\t"<<seq<<'\t'<<q.sources.size()<<'\t'<<path_text(e.seed.path)<<'\n';
                seed_fasta<<'>'<<e.id<<'\n'<<seq<<'\n';size_t rank=0;
                for(const auto& c:e.paths) {
                    ++rank;auto key=canonical(c.sequence);auto id=identifier(key,"region_");
                    auto [it,inserted]=regions.emplace(id,std::make_pair(key,c.sequence));
                    if(it->second.first!=key)throw std::runtime_error("region identifier collision");
                    if(inserted) {
                        Path edges;for(size_t i=1;i<c.path.size();++i)edges.push_back(graph.edge(c.path[i-1],c.path[i]));
                        regions_out<<id<<'\t'<<c.sequence.size()<<'\t'<<c.sequence<<'\t'<<path_text(c.path)<<'\t'<<path_text(edges)<<"\tnot_assessed\n";
                        fasta<<'>'<<id<<'\n'<<c.sequence<<'\n';
                    }
                    bool flip=c.sequence!=it->second.second;
                    auto stable=stability(c.left.coverage,c.right.coverage);
                    const auto& lp=(flip?c.right:c.left).coverage.profile;const auto& rp=(flip?c.left:c.right).coverage.profile;
                    bool left_enough=lp.samples>=options.coverage_window,right_enough=rp.samples>=options.coverage_window;
                    selection<<query_id<<'\t'<<e.id<<'\t'<<id<<'\t'<<q.compared_paths<<'\t'<<q.compared_regions<<'\t'<<c.score<<'\t'
                        <<stable.score()<<'\t'<<stable.within<<'\t'<<stable.jumps<<'\t'<<stable.between<<'\t'<<stable.support<<'\t'
                        <<lp.samples<<'\t'<<rp.samples<<'\t'<<lp.level()<<'\t'<<rp.level()<<'\t'<<double(e.seed.end-e.seed.begin)/q.sequence.size()<<'\t'
                        <<(left_enough&&right_enough?"both_flanks":left_enough||right_enough?"one_flank":"insufficient")<<'\t';
                    if(q.has_runner)selection<<q.runner_score;else selection<<"NA";
                    selection<<'\t';if(q.has_runner && q.comparable_runner)selection<<c.score-q.runner_score;else selection<<"NA";
                    selection<<'\t'<<q.score_tie<<'\t'<<(c.search_limited||q.truncated)<<'\n';
                    for(const auto& source:q.sources) {
                        ++links;size_t qb=source.reverse?q.sequence.size()-e.seed.end:e.seed.begin;
                        size_t qe=source.reverse?q.sequence.size()-e.seed.begin:e.seed.end;
                        size_t sb=flip?c.sequence.size()-c.seed_end:c.seed_begin;
                        size_t se=flip?c.sequence.size()-c.seed_begin:c.seed_end;
                        auto pam=pam_status(q,e,c,source,rules);
                        bool l=c.external_left>0,r=c.external_right>0;
                        associations<<source.id<<'\t'<<source.array<<'\t'<<e.id<<'\t'<<id<<'\t'<<qb<<'\t'<<qe<<'\t'<<sb<<'\t'<<se<<'\t'<<((source.reverse^flip)?"-":"+")<<'\t'
                            <<(l&&r?"two_sided":"one_sided")<<'\t'<<(flip?c.right.stop:c.left.stop)<<'\t'<<(flip?c.left.stop:c.right.stop)<<'\t'
                            <<(flip?c.right.coverage:c.left.coverage).mean()<<'\t'<<(flip?c.left.coverage:c.right.coverage).mean()<<'\t'<<(flip?c.right.coverage:c.left.coverage).dispersion()<<'\t'<<(flip?c.left.coverage:c.right.coverage).dispersion()<<'\t'
                            <<c.left.coverage.scored+c.right.coverage.scored<<'\t'<<c.left.coverage.score()+c.right.coverage.score()<<'\t'<<(c.score-c.pam_bonus+(pam.first=="match"?options.pam_weight:0))<<'\t'<<rank<<'\t'<<(c.search_limited||q.truncated)<<'\t'<<c.expanded<<'\t'
                            <<(reads1.empty()?"not_provided":"exact_retained_candidates")<<'\t'<<c.supported_junctions<<'\t'<<c.junctions<<'\t'<<c.read_support<<'\t'<<c.pair_support<<'\t'<<c.unique_pairs<<'\t'<<pam.first<<'\t'<<pam.second<<'\t'<<(flip?c.external_right:c.external_left)<<'\t'<<(flip?c.external_left:c.external_right)<<'\t'<<(q.score_tie?"selected_score_tie":q.compared_regions>1?"selected_best_of_alternatives":(c.search_limited||q.truncated)?"search_limited":"single_candidate")<<'\t'<<(flip?c.right.context:c.left.context)<<'\t'<<(flip?c.left.context:c.right.context)<<'\n';
                    }
                }
            }
        }
        auto parameters=output(staging/"protospacer_parameters.json");
        std::string adapters_json="[";
        for(const auto& adapter:options.adapters) {
            if(adapters_json.size()>1)adapters_json+=",";
            adapters_json+="{\"sequence\":"+json_string(adapter.sequence)+",\"exact_begin\":"+std::to_string(adapter.exact_begin)+
                ",\"exact_end\":"+std::to_string(adapter.exact_end)+",\"max_mismatches\":"+std::to_string(adapter.max_mismatches)+"}";
        }adapters_json+="]";
        parameters<<"{\n  \"schema\": 3,\n  \"selection\": \"one_prediction_per_distinct_spacer\",\n  \"ranking\": \"log2_multiplicity_stability\",\n  \"graph\": "<<json_string(fs::absolute(graph_path).string())<<",\n  \"arrays\": "<<json_string(fs::absolute(arrays).string())
            <<",\n  \"max_length\": "<<options.max_bp<<",\n  \"beam_width\": "<<options.beam_width<<",\n  \"max_states_per_side\": "<<options.max_states
            <<",\n  \"max_paths_per_seed\": "<<options.max_paths<<",\n  \"max_seeds_per_query\": "<<options.max_seeds<<",\n  \"max_query_states\": "<<options.max_query_states<<",\n  \"min_match_bp\": "<<std::max(graph.k(),options.min_match)
            <<",\n  \"max_retained_bases\": "<<max_retained_bases
            <<",\n  \"coverage_model\": "<<json_string(options.coverage_model)<<",\n  \"coverage_window\": "<<options.coverage_window<<",\n  \"max_coverage_penalty\": "<<options.max_coverage_penalty
            <<",\n  \"array_repeat_identity\": "<<options.array_repeat_identity
            <<",\n  \"adapter_patterns\": "<<adapters_json<<",\n  \"adapter_match\": \"anchored_hamming_both_strands\""
            <<",\n  \"adapter_blocked_branches\": "<<adapter_branches<<",\n  \"adapter_diagnostic_rows\": "<<adapter_rows
            <<",\n  \"adapter_contaminated_queries\": "<<adapter_queries
            <<",\n  \"threads\": "<<threads<<",\n  \"reads1\": "<<json_string(reads1.string())<<",\n  \"reads2\": "<<json_string(reads2.string())
            <<",\n  \"min_insert\": "<<min_insert<<",\n  \"max_insert\": "<<max_insert<<",\n  \"max_read_index\": "<<max_index<<",\n  \"max_read_work\": "<<max_read_work
            <<",\n  \"read_weight\": "<<options.read_weight<<",\n  \"pair_weight\": "<<options.pair_weight<<",\n  \"pam_rules\": "<<json_string(pam_file.string())<<",\n  \"pam_weight\": "<<options.pam_weight
            <<",\n  \"source_observations\": "<<sources<<",\n  \"unique_queries\": "<<queries.size()<<",\n  \"matched_queries\": "<<matched_queries<<",\n  \"full_length_matched_queries\": "<<full_matched_queries
            <<",\n  \"partial_only_matched_queries\": "<<matched_queries-full_matched_queries<<",\n  \"unresolved_queries\": "<<unresolved_queries<<",\n  \"unresolved_seed_intervals\": "<<unresolved_seeds
            <<",\n  \"rejected_context_seed_intervals\": "<<rejected_context_seeds<<",\n  \"redundant_interval_paths\": "<<redundant_paths
            <<",\n  \"seed_candidates\": "<<seeds<<",\n  \"regions\": "<<regions.size()<<",\n  \"associations\": "<<links<<"\n}\n";
        query_out.close();seed_out.close();seed_fasta.close();regions_out.close();fasta.close();associations.close();parameters.close();
        matched_fasta.close();unresolved_fasta.close();unresolved_out.close();decisions.close();selection.close();adapter_stops.close();
        for(const auto& item:fs::directory_iterator(staging))fs::rename(item.path(),destination/item.path().filename());
        fs::remove(staging);
        std::cout<<"Matched "<<matched_queries<<" / "<<queries.size()<<" distinct spacers ("<<full_matched_queries<<" with full-length matches; "<<matched_queries-full_matched_queries<<" partial-only)\n"
            <<"Selected one prediction per matched spacer: "<<seeds<<" predictions, "<<regions.size()<<" unique regions, "<<links<<" source associations to "<<destination<<'\n'
            <<"Separate unresolved intervals: "<<unresolved_seeds<<"; source-array intervals rejected: "<<rejected_context_seeds<<"; redundant interval paths suppressed: "<<redundant_paths<<'\n';
        std::cout<<"Adapter branches stopped: "<<adapter_branches<<"; adapter-containing input spacers rejected: "<<adapter_queries<<'\n';
    }catch(...){std::error_code ec;fs::remove_all(staging,ec);throw;}
    return 0;
}
} // namespace mcaat::proto
