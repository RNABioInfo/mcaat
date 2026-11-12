#include "protospacer_detection/sdbg_extension_graph.h"
#include "complete_target_placements.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <omp.h>
using namespace mcaat::proto;
// BOSS colexicographic intervals let us prune whole ranges of k-mers without
// materializing the graph. Each leaf is checked against its actual node label.
class PlacementGraph:public ExtensionGraph {
    const SDBG& raw_;
public:
    explicit PlacementGraph(const SDBG& raw):ExtensionGraph(raw),raw_(raw){}
    std::vector<std::pair<Node,size_t>> find_masked(const std::string& pattern,size_t budget)const {
        // Read2SdbgS2 emits every edge together with its reverse complement.
        // For that unmodified double-stranded graph, search the orientation
        // with fewer trailing wildcards, then map the node labels back.
        auto search_pattern=pattern;bool mirrored=false;
        if(std::getenv("MCAAT_AUDIT_BOTH_STRANDS")) {
            auto leading=pattern.find_first_not_of('.'),last=pattern.find_last_not_of('.');
            size_t trailing=last==std::string::npos?pattern.size():pattern.size()-last-1;
            if(trailing>leading) {
                mirrored=true;std::reverse(search_pattern.begin(),search_pattern.end());
                for(auto& c:search_pattern)if(c!='.')c="TGCA"[std::string("ACGT").find(c)];
            }
        }
        std::map<Node,size_t> result;std::vector<uint8_t> bases(k());
        std::function<void(Node,Node,size_t,size_t)> visit;
        visit=[&](Node lo,Node hi,size_t depth,size_t errors) {
            if(lo>=hi)return;
            if(depth==k() || hi-lo<=32) {
                for(Node index=lo;index<hi;) {
                    Node v=raw_.GetLastIndex(index);raw_.GetLabel(v,bases.data());std::string word;
                    for(auto b:bases){if(b<1||b>4){word.clear();break;}word+="ACGT"[b-1];}
                    if(!word.empty()){size_t actual=0;for(size_t j=0;j<k();++j)actual+=search_pattern[j]!='.'&&search_pattern[j]!=word[j];auto mapped=find(word);if(mapped!=missing&&actual<=budget)result.emplace(mapped,actual);}
                    if(v<index)throw std::runtime_error("invalid BOSS range");index=v+1;
                }
                return;
            }
            const size_t at=k()-1-depth;
            auto bound=[&](unsigned base) {
                auto a=lo,b=hi;
                while(a<b){auto mid=a+(b-a)/2;raw_.GetLabel(mid,bases.data());if(bases[at]<base)a=mid+1;else b=mid;}
                return a;
            };
            for(unsigned base=1;base<=4;++base) {
                size_t cost=errors+(search_pattern[at]!='.'&&search_pattern[at]!="ACGT"[base-1]);
                if(cost>budget)continue;auto first=bound(base),last=bound(base+1);
                visit(first,last,depth+1,cost);
            }
        };
        visit(0,raw_.size(),0,0);
        if(mirrored) {
            std::map<Node,size_t> forward;
            for(auto [node,errors]:result) {
                auto original=find(reverse_complement(label(node)));
                if(original==missing)throw std::runtime_error("double-stranded lookup encountered an asymmetric vertex");
                forward.emplace(original,errors);
            }
            result=std::move(forward);
        }
        if(std::getenv("MCAAT_AUDIT_CHECK_INDEX")) {
            if(raw_.size()>100000)throw std::runtime_error("index brute check only for tiny fixtures");
            std::map<Node,size_t> expected;
            for(Node i=0;i<raw_.size();){auto v=raw_.GetLastIndex(i);raw_.GetLabel(v,bases.data());std::string word;
                for(auto b:bases){if(b<1||b>4){word.clear();break;}word+="ACGT"[b-1];}
                if(!word.empty()){auto mapped=find(word);size_t errors=0;for(size_t j=0;j<k();++j)errors+=pattern[j]!='.'&&pattern[j]!=word[j];if(mapped!=missing&&errors<=budget)expected.emplace(mapped,errors);}
                i=v+1;
            }
            if(result!=expected)throw std::runtime_error("BOSS masked index differs from full graph scan");
        }
        return {result.begin(),result.end()};
    }
};
int main(int argc,char** argv){try {
    if(argc!=3)throw std::runtime_error("usage: audit_complete_placements GRAPH_PREFIX DIR");
    std::string dir=argv[2];omp_set_num_threads(2);SDBG raw;raw.LoadFromFile(argv[1]);PlacementGraph graph(raw);
    std::ifstream input(dir+"/targets.tsv");std::ofstream output(dir+"/complete_placements.tsv");
    if(!input||!output)throw std::runtime_error("missing target table/output directory");
    output<<"seed_id\tside\tvertex\tobserved_target\tcontext\tflank_prefix\tquery_begin\tquery_end\tmismatches\n";
    std::string line;std::getline(input,line);std::vector<std::pair<std::string,std::string>> jobs;
    while(std::getline(input,line)) {
        if(!line.empty()&&line.back()=='\r')line.pop_back();
        std::istringstream row(line);std::string id,query;std::getline(row,id,'\t');std::getline(row,query,'\t');
        jobs.emplace_back(id,query);
    }
    std::exception_ptr error;omp_set_num_threads(4);
    #pragma omp parallel for schedule(dynamic,1)
    for(size_t job=0;job<jobs.size()*2;++job)try {
            const auto& [id,query]=jobs[job/2];bool back=job%2;std::ostringstream rows;
            auto found=audit::complete_target_placements(graph,back?reverse_complement(query):query);
            for(const auto& [key,p]:found.endpoints) {
                auto v=p.endpoint;
                auto outward=p.flank_prefix;if(back)for(auto& c:outward)c="TGCA"[std::string("ACGT").find(c)];
                auto node=back?graph.find(reverse_complement(graph.label(v))):v;
                if(node==missing)throw std::runtime_error("reverse-complement endpoint missing");
                rows<<id<<'\t'<<(back?"left":"right")<<'\t'<<node<<'\t'<<(back?reverse_complement(p.observed):p.observed)
                      <<'\t'<<(back?reverse_complement(p.context):p.context)<<'\t'<<outward<<'\t'<<(back?query.size()-p.end:p.begin)
                      <<'\t'<<(back?query.size()-p.begin:p.end)<<'\t'<<p.errors<<'\n';
            }
            #pragma omp critical
            {output<<rows.str()<<std::flush;std::cout<<id<<' '<<(back?"left":"right")<<": "<<found.endpoints.size()<<" endpoints, "<<found.lookups<<" lookups, "<<found.states<<" states"<<std::endl;}
    }catch(...){
        #pragma omp critical
        {if(!error)error=std::current_exception();}
    }
    if(error)std::rethrow_exception(error);
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
