#include "crispr_array/graph_repeat_filter.h"
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <unordered_map>
using namespace mcaat::arrays;
struct Graph {
    size_t k=23;
    std::unordered_map<std::string,Node> ids;
    std::vector<std::string> labels;
    std::vector<std::vector<Node>> out;
    Node node(const std::string& s) {
        auto it=ids.find(s);if(it!=ids.end())return it->second;
        Node id=labels.size();ids.emplace(s,id);labels.push_back(s);out.emplace_back();return id;
    }
    Path add_one(const std::string& period) {
        Path path;auto padded=period+period.substr(0,k-1);
        for(size_t i=0;i<period.size();++i)path.push_back(node(padded.substr(i,k)));
        for(size_t i=0;i<path.size();++i) {
            auto& edges=out[path[i]];Node to=path[(i+1)%path.size()];
            if(std::find(edges.begin(),edges.end(),to)==edges.end())edges.push_back(to);
        }
        return path;
    }
    Path add(const std::string& period) {
        auto p=add_one(period);add_one(reverse_complement(period));return p;
    }
    const std::string& label(Node n) const{return labels.at(n);}
    bool adjacent(Node a,Node b) const {
        auto& edges=out[a];return std::find(edges.begin(),edges.end(),b)!=edges.end();
    }
    Node reverse_node(Node n) const {
        auto it=ids.find(reverse_complement(label(n)));
        return it==ids.end()?missing_node:it->second;
    }
};
std::string dna(std::mt19937& rng,size_t n) {
    std::string result;while(n--)result+="ACGT"[rng()%4];return result;
}
void write_path(std::ostream& out,const Path& p) {for(Node n:p)out<<n<<',';out<<'\t';}
int main(int argc,char**argv) {
    try {
        if(argc<4)throw std::runtime_error("benchmark_repeat_filter FAMILIES GROUPS_PER_FAMILY OUTPUT [seed=42] [identical|mixed|variants|ambiguous] [support=0.85]");
        size_t families=std::stoull(argv[1]),copies=std::stoull(argv[2]);
        if(!families || !copies)throw std::runtime_error("families and groups must be positive");
        std::mt19937 rng(argc>4?std::stoul(argv[4]):42);
        const std::string mode=argc>5?argv[5]:"identical";
        if(mode!="identical" && mode!="mixed" && mode!="variants" && mode!="ambiguous")
            throw std::runtime_error("unknown fixture mode");
        const bool mixed=mode!="identical",variants=mode=="variants",ambiguous=mode=="ambiguous";
        Graph graph;Options options;options.k=23;
        if(argc>6)options.repeat_support=parse_repeat_support(argv[6]);
        if(ambiguous) {graph.k=9;options.k=9;options.repeat_min=12;options.spacer_min=12;options.spacer_max=80;}
        std::map<Node,std::vector<Path>> input;
        for(size_t family=0;family<families;++family) {
            auto repeat=dna(rng,ambiguous?20:30),other=dna(rng,20);
            std::vector<Path> cycles;
            for(size_t i=0;i<32;++i) {
                auto spacer=dna(rng,30);spacer.front()="ACGT"[i%4];spacer.back()="TGCA"[i%4];
                auto rep=repeat;
                if(variants && i>=24) {
                    const size_t pos=i%3==0?0:(i%3==1?rep.size()-1:14);
                    rep[pos]=rep[pos]=='A'?'C':'A';
                }
                if(ambiguous) {
                    // Two potentially conserved regions with overlapping but
                    // different support. Subgroups can select different roles.
                    if(i>=20)rep=dna(rng,20);
                    spacer=(i<12?dna(rng,20):other)+dna(rng,15);
                }
                cycles.push_back(graph.add(rep+spacer));
            }
            for(size_t j=0;j<copies;++j) {
                auto& group=input[family*copies+j];
                for(size_t i=0;i<cycles.size();++i) {
                    if(mixed && rng()%4<(mode=="mixed"?1:2))continue;
                    auto c=rotate(cycles[i],(j*7)%cycles[i].size());
                    if(mixed && j%2) {std::reverse(c.begin(),c.end());for(auto& n:c)n=graph.reverse_node(n);}
                    group.push_back(std::move(c));
                }
            }
        }
        RepeatFilter<Graph> filter(graph,options);
        auto start=std::chrono::steady_clock::now();
        auto systems=filter.run(input);
        double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        std::ofstream out(argv[3]);
        out.exceptions(std::ios::badbit|std::ios::failbit);
        for(const auto& s:systems) {
            write_path(out,s.repeat_path);out<<s.support_count<<'\t'<<s.alternate_repeat_paths<<'\t';write_path(out,s.sources);out<<'\n';
            for(const auto& c:s.candidate_cycles){out<<"C\t";write_path(out,c);out<<'\n';}
            for(const auto& row:s.observations){out<<"O\t";write_path(out,row.cycle);out<<row.repeat_nodes<<'\t'<<row.repeat<<'\t'<<row.spacer<<'\t'<<row.supports_path<<'\n';}
        }
        const auto& s=filter.stats;
        out<<"stats\t"<<s.input_cycles<<'\t'<<s.invalid_cycles<<'\t'<<s.duplicate_cycles<<'\t'<<s.input_groups<<'\t'<<s.groups_without_repeat<<'\t'<<s.duplicate_systems<<'\t'<<s.merged_systems<<'\n';
        std::cout<<std::setprecision(8)<<"groups="<<input.size()<<" cycles="<<s.input_cycles<<" systems="<<systems.size()<<" seconds="<<elapsed<<'\n';
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
