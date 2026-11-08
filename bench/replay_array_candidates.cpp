// Re-run final curation on already inferred candidates; NOT cycle/repeat inference.
#include "crispr_array/array_output.h"
#include <unordered_map>
using namespace mcaat::arrays;
using Fields=std::map<std::string,std::string>;
std::vector<std::string> split(const std::string& s,char delimiter) {
    std::vector<std::string> out;std::istringstream in(s);std::string token;
    while(std::getline(in,token,delimiter)) out.push_back(token);return out;
}
Path path(const std::string& s) {
    Path out;if(s=="-" || s.empty())return out;
    for(const auto& v:split(s,','))out.push_back(std::stoull(v));return out;
}
template<class F> void rows(const std::filesystem::path& file,F f) {
    std::ifstream in(file);if(!in)throw std::runtime_error("cannot read "+file.string());
    std::string line;std::getline(in,line);auto names=split(line,'\t');
    while(std::getline(in,line)) {
        if(line.empty())continue;auto values=split(line,'\t');
        if(values.size()!=names.size())throw std::runtime_error("column mismatch: "+file.string());
        Fields row;for(size_t i=0;i<names.size();++i)row[names[i]]=values[i];f(row);
    }
    if(in.bad())throw std::runtime_error("read failed: "+file.string());
}
struct ReportGraph {
    size_t k;
    std::unordered_map<Node,std::string> labels;
    std::map<std::pair<Node,Node>,Node> links;
    const std::string& label(Node n) {return labels.at(n);}
    Path edges(const Path& c) {
        Path result;for(size_t i=0;i<c.size();++i)result.push_back(links.at({c[i],c[(i+1)%c.size()]}));return result;
    }
    void add(const Observation& row,const Path& edge_ids) {
        const auto sequence=row.repeat+row.spacer,padded=sequence+sequence.substr(0,k-1);
        if(row.cycle.size()!=sequence.size() || edge_ids.size()!=row.cycle.size())throw std::runtime_error("cycle lengths disagree");
        for(size_t i=0;i<row.cycle.size();++i) {
            auto [it,inserted]=labels.emplace(row.cycle[i],padded.substr(i,k));
            if(!inserted && it->second!=padded.substr(i,k))throw std::runtime_error("conflicting node labels");
            auto [edge,fresh]=links.emplace(std::make_pair(row.cycle[i],row.cycle[(i+1)%row.cycle.size()]),edge_ids[i]);
            if(!fresh && edge->second!=edge_ids[i])throw std::runtime_error("conflicting edge IDs");
        }
    }
};
int main(int argc,char** argv) {
    try {
        if(argc<3)throw std::runtime_error("replay_array_candidates INPUT_REPORT_DIR OUTPUT_DIR [--array-filter-mode conservative|strict] [--cycle-min-length 47 --cycle-max-length 100]");
        namespace fs=std::filesystem;fs::path input=argv[1],output=argv[2];
        if(fs::weakly_canonical(input)==fs::weakly_canonical(output))throw std::runtime_error("replay output must differ from input");
        Options options;options.cycle_min=47;options.cycle_max=100;
        rows(input/"CRISPR_graph_summary.tsv",[&](const Fields& r){
            const auto& key=r.at("metric");const auto& value=r.at("value");
            if(key=="k")options.k=std::stoull(value);
            else if(key=="repeat_support")options.repeat_support=parse_repeat_support(value);
            else if(key=="repeat_min_length")options.repeat_min=std::stoull(value);
            else if(key=="repeat_max_length")options.repeat_max=std::stoull(value);
            else if(key=="spacer_min_length")options.spacer_min=std::stoull(value);
            else if(key=="spacer_max_length")options.spacer_max=std::stoull(value);
        });
        for(int i=3;i<argc;++i) {
            std::string flag=argv[i];if(++i>=argc)throw std::runtime_error("missing value: "+flag);
            if(flag=="--array-filter-mode")options.quality.strict=parse_filter_mode(argv[i]);
            else if(flag=="--cycle-min-length" || flag=="--cycle-max-length" ||
                    flag=="--min-repeat-len" || flag=="--max-repeat-len" ||
                    flag=="--min-spacer-len" || flag=="--max-spacer-len") {
                size_t used=0;auto n=std::stoull(argv[i],&used);
                if(used!=std::string(argv[i]).size() || std::string(argv[i]).find_first_not_of("0123456789")!=std::string::npos || !n || n>std::numeric_limits<int>::max())throw std::runtime_error("invalid length: "+flag);
                if(flag=="--cycle-min-length")options.cycle_min=n;
                else if(flag=="--cycle-max-length")options.cycle_max=n;
                else if(flag=="--min-repeat-len")options.repeat_min=n;
                else if(flag=="--max-repeat-len")options.repeat_max=n;
                else if(flag=="--min-spacer-len")options.spacer_min=n;
                else options.spacer_max=n;
            } else throw std::runtime_error("unknown option: "+flag);
        }
        options.validate();
        ReportGraph graph{options.k};std::vector<System> systems;std::map<std::string,size_t> ids;
        std::vector<size_t> expected;
        rows(input/"CRISPR_graph_systems.tsv",[&](const Fields& r){
            if(!ids.emplace(r.at("array_id"),systems.size()).second)throw std::runtime_error("duplicate array ID");
            System s;s.repeat_path=path(r.at("repeat_vertices"));s.sources=path(r.at("source_start_edges"));
            s.support_count=std::stoull(r.at("support_cycles"));s.candidate_count=std::stoull(r.at("candidate_cycles"));
            s.alternate_repeat_paths=std::stoull(r.at("alternative_repeat_paths"));systems.push_back(std::move(s));
            expected.push_back(std::stoull(r.at("assigned_cycles")));
        });
        size_t total=0;
        rows(input/"CRISPR_graph_paths.tsv",[&](const Fields& r){
            auto& s=systems.at(ids.at(r.at("array_id")));Observation row;
            row.cycle=path(r.at("cycle_vertices"));row.repeat=r.at("repeat");row.spacer=r.at("spacer");
            row.repeat_nodes=path(r.at("repeat_vertices")).size();row.supports_path=r.at("supports_repeat_path")=="1";
            if(row.repeat_nodes+options.k-1!=row.repeat.size())throw std::runtime_error("repeat length mismatch");
            graph.add(row,path(r.at("cycle_edges")));s.candidate_cycles.push_back(cycle_key(row.cycle));s.observations.push_back(std::move(row));++total;
        });
        for(size_t i=0;i<systems.size();++i) {
            size_t support=0;for(const auto& r:systems[i].observations)support+=r.supports_path;
            if(expected[i]!=systems[i].observations.size() || support!=systems[i].support_count)throw std::runtime_error("report counts disagree");
        }
        RepeatFilter<ReportGraph> filter(graph,options);filter.stats.input_cycles=total;filter.stats.input_groups=systems.size();
        auto started=std::chrono::steady_clock::now();
        auto accepted=filter.curate(std::move(systems));
        const double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
        write_array_outputs(graph,filter,options,accepted,output.string());
        auto scope=array_output(output/"REPLAY_SCOPE.txt");
        scope<<"Final curation of "<<expected.size()<<" previously retained candidates / "<<total<<" assigned observations.\n"
             <<"Not a replay of original cycle search or repeat inference. Rejected earlier groups and unresolved paths are absent.\n"
             <<"Original candidate-cycle denominators, including unresolved counts, are retained.\n"
             <<"candidate_id N refers to original input table row N; original array IDs:\n";
        for(const auto& [id,index]:ids)scope<<index+1<<'\t'<<id<<'\n';
        std::cout<<"curation_seconds="<<elapsed<<'\n';
    } catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
