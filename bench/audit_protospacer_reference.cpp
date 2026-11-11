#include "protospacer_detection/sdbg_extension_graph.h"
#include "reference_graph_oracle.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <omp.h>
using namespace mcaat::proto;
using Row=std::map<std::string,std::string>;
std::vector<std::string> split(std::string s){if(!s.empty()&&s.back()=='\r')s.pop_back();std::vector<std::string> r;size_t a=0;while(true){auto b=s.find('\t',a);r.push_back(s.substr(a,b-a));if(b==std::string::npos)return r;a=b+1;}}
std::vector<Row> table(const std::string& path){std::ifstream f(path);if(!f)throw std::runtime_error("cannot open "+path);std::string line;std::getline(f,line);auto keys=split(line);std::vector<Row> rows;while(std::getline(f,line)){auto v=split(line);if(v.size()!=keys.size())throw std::runtime_error("bad table");Row r;for(size_t i=0;i<keys.size();++i)r[keys[i]]=v[i];rows.push_back(std::move(r));}return rows;}
std::string choices(ExtensionGraph& g,Node v,bool back){std::ostringstream s;for(auto a:g.neighbors(v,back))s<<g.extension_base(a,back)<<':'<<a.count<<':'<<a.to<<',';return s.str();}
template<class V>std::string joined(const V& values){std::ostringstream s;for(size_t i=0;i<values.size();++i){if(i)s<<',';s<<values[i];}return s.str();}
int main(int argc,char** argv){try{
    if(argc<3)throw std::runtime_error("usage: audit_reference GRAPH_PREFIX DIR [max_states] [case_list]");
    std::string dir=argv[2];size_t limit=argc>3?std::stoull(argv[3]):100000;
    std::set<std::string> filter;if(argc>4){std::ifstream in(argv[4]);std::string s;while(in>>s)filter.insert(s);}
    auto cases=table(dir+"/reference_cases.tsv"),targets=table(dir+"/targets.tsv");
    std::map<std::string,std::string> sequences;for(auto r:targets)sequences[r["seed_id"]]=r["sequence"];
    omp_set_num_threads(4);std::cout<<"Loading existing graph\n"<<std::flush;SDBG raw;raw.LoadFromFile(argv[1]);ExtensionGraph g(raw);
    std::map<std::string,Path> seeds;for(auto [id,seq]:sequences)seeds[id]=g.trace(seq);
    std::ofstream exact(dir+"/exact_native_routes.tsv"),summary(dir+"/oracle_routes.tsv"),paths(dir+"/oracle_paths.tsv");
    exact<<"case_id\tseed_id\tside\treference\tseed_length\tseed_mapped\treference_seed_exact\treference_boundary_clipped\tavailable_reference_bp\tknown_reference_bp\texact_bp\tfirst_stop\tvertex\texpected_base\tavailable_edges\n";
    summary<<"case_id\tseed_id\tside\treference\tcheckpoint_ref_bp\tstatus\tquery_bp\tedit_cost\tidentity\tstates_discovered\tstates_popped\tlimited\texhausted\tlowest_unexplored_cost\n";
    paths<<"case_id\tseed_id\tside\tcheckpoint_ref_bp\tquery_bp\tsequence_outward\toperations\tvertices\tcounts\n";
    size_t done=0,totalstates=0;
    for(const auto& r:cases){if(!filter.empty()&&!filter.count(r.at("case_id")))continue;
        auto sid=r.at("seed_id");const auto& seed=seeds.at(sid);
        for(bool back:{true,false}){
            auto ref=r.at(back?"left_outward":"right_outward");size_t available=ref.size();auto unknown=ref.find_first_not_of("ACGT");if(unknown!=std::string::npos)ref.resize(unknown);
            Node start=seed.empty()?missing:back?seed.front():seed.back(),v=start;size_t n=0;
            if(v!=missing)while(n<ref.size()){
                Node next=missing;for(auto a:g.neighbors(v,back))if(g.extension_base(a,back)==ref[n]){next=a.to;break;}
                if(next==missing)break;v=next;++n;
            }
            std::string reason=start==missing?(sequences[sid].size()<g.k()?"input_shorter_than_k":"input_seed_unmapped"):n<ref.size()?"exact_edge_absent":unknown!=std::string::npos?"reference_ambiguous_base":available<750?"reference_end":"750_exact_bases";
            exact<<r.at("case_id")<<'\t'<<sid<<'\t'<<(back?"left":"right")<<'\t'<<r.at("reference")<<'\t'<<sequences[sid].size()<<'\t'<<!seed.empty()<<'\t'<<(sequences[sid]==r.at("reference_allele"))<<'\t'<<r.at("reference_boundary_clipped")<<'\t'<<available<<'\t'<<ref.size()<<'\t'<<n<<'\t'<<reason<<'\t'<<v<<'\t'<<(n<ref.size()?ref.substr(n,1):"-")<<'\t'<<(v!=missing?choices(g,v,back):"-")<<'\n';
            std::set<size_t> checkpoints;for(size_t t:{24,100,250,500,750})if(t<=ref.size())checkpoints.insert(t);if(!ref.empty())checkpoints.insert(ref.size());
            audit::Search result;if(start!=missing&&!checkpoints.empty())result=audit::align(g,start,ref,back,checkpoints,limit);
            totalstates+=result.discovered;
            if(checkpoints.empty())checkpoints.insert(0);
            for(auto t:checkpoints){auto it=result.found.find(t);bool found=it!=result.found.end();
                std::string status=found?"verified_graph_path":start==missing?reason:t==0?reason:result.lowest_unexplored_cost<=t?"search_limit_unresolved":"no_path_at_90pct_anchored_identity";
                summary<<r.at("case_id")<<'\t'<<sid<<'\t'<<(back?"left":"right")<<'\t'<<r.at("reference")<<'\t'<<t<<'\t'<<status<<'\t';
                if(found){auto& w=it->second;summary<<w.sequence.size()<<'\t'<<w.edits<<'\t'<<100.*w.matches/w.columns;
                    std::vector<double> counts;for(size_t i=1;i<w.nodes.size();++i)counts.push_back(back?g.count(w.nodes[i],w.nodes[i-1]):g.count(w.nodes[i-1],w.nodes[i]));
                    paths<<r.at("case_id")<<'\t'<<sid<<'\t'<<(back?"left":"right")<<'\t'<<t<<'\t'<<w.sequence.size()<<'\t'<<w.sequence<<'\t'<<w.operations<<'\t'<<joined(w.nodes)<<'\t'<<joined(counts)<<'\n';
                }else summary<<"0\tNA\tNA";
                summary<<'\t'<<result.discovered<<'\t'<<result.popped<<'\t'<<result.limited<<'\t'<<result.exhausted<<'\t'<<result.lowest_unexplored_cost<<'\n';
            }
        }
        if(++done%20==0){exact.flush();summary.flush();paths.flush();std::cout<<"Audited "<<done<<" reference locations; "<<totalstates<<" diagnostic states\n"<<std::flush;}
    }
    std::cout<<"Complete: "<<done<<" reference locations\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
