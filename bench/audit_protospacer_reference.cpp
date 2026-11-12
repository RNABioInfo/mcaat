#include "protospacer_detection/sdbg_extension_graph.h"
#include "reference_graph_oracle.h"
#include "reference_unitig_oracle.h"
#include "reference_rolling_oracle.h"
#include "reference_segment_bounds.h"
#include "short_seed_placements.h"
#include "reachable_membership.h"
#include "adapter_filtered_graph.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <future>
#include <mutex>
#include <omp.h>
using namespace mcaat::proto;
using Row=std::map<std::string,std::string>;
std::vector<std::string> split(std::string s){if(!s.empty()&&s.back()=='\r')s.pop_back();std::vector<std::string> r;size_t a=0;while(true){auto b=s.find('\t',a);r.push_back(s.substr(a,b-a));if(b==std::string::npos)return r;a=b+1;}}
std::vector<Row> table(const std::string& path){std::ifstream f(path);if(!f)throw std::runtime_error("cannot open "+path);std::string line;std::getline(f,line);auto keys=split(line);std::vector<Row> rows;while(std::getline(f,line)){auto v=split(line);if(v.size()!=keys.size())throw std::runtime_error("bad table");Row r;for(size_t i=0;i<keys.size();++i)r[keys[i]]=v[i];rows.push_back(std::move(r));}return rows;}
template<class Graph>std::string choices(Graph& g,Node v,bool back){std::ostringstream s;for(auto a:g.neighbors(v,back))s<<g.extension_base(a,back)<<':'<<a.count<<':'<<a.to<<',';return s.str();}
template<class V>std::string joined(const V& values){std::ostringstream s;for(size_t i=0;i<values.size();++i){if(i)s<<',';s<<values[i];}return s.str();}
int main(int argc,char** argv){try{
    if(argc<3)throw std::runtime_error("usage: audit_reference GRAPH_PREFIX DIR [max_states] [case_list]");
    std::string dir=argv[2];size_t limit=argc>3?std::stoull(argv[3]):100000;
    bool any_length=std::getenv("MCAAT_AUDIT_ANY_LENGTH")!=nullptr;
    bool use_reachability=std::getenv("MCAAT_AUDIT_REACHABLE_BOUNDS")!=nullptr;
    bool prefer_progress=std::getenv("MCAAT_AUDIT_PREFER_PROGRESS")!=nullptr;
    bool filter_adapters=std::getenv("MCAAT_AUDIT_FILTER_ADAPTERS")!=nullptr;
    bool strong_bounds=std::getenv("MCAAT_AUDIT_STRONG_BOUNDS")!=nullptr;
    bool unitigs=std::getenv("MCAAT_AUDIT_UNITIGS")!=nullptr;
    bool exact_only=std::getenv("MCAAT_AUDIT_EXACT_ONLY")!=nullptr;
    bool rolling=std::getenv("MCAAT_AUDIT_ROLLING")!=nullptr;
    std::set<std::string> filter;if(argc>4){std::ifstream in(argv[4]);std::string s;while(in>>s)filter.insert(s);}
    auto cases=table(dir+"/reference_cases.tsv"),targets=table(dir+"/targets.tsv");
    std::map<std::string,std::string> sequences;for(auto r:targets)sequences[r["seed_id"]]=r["sequence"];
    std::map<std::pair<std::string,bool>,std::vector<FlankContext>> source_models;
    std::map<std::string,std::vector<FlankContext>> cycle_models;
    const char* cycles_path=std::getenv("MCAAT_AUDIT_CYCLES");
    if(cycles_path)for(const auto& row:table(cycles_path))cycle_models[row.at("seed_id")].push_back({"","",row.at("cycle")});
    const char* contexts_path=std::getenv("MCAAT_AUDIT_CONTEXTS");
    if(contexts_path) {
        for(const auto& row:table(contexts_path))for(bool back:{true,false}) {
            auto id=row.at("seed_id");auto begin=std::stoull(row.at("begin")),end=std::stoull(row.at("end"));
            if(row.at("spacer").substr(begin,end-begin)!=sequences.at(id))throw std::runtime_error("wrong source context for input");
            auto models=flank_contexts(row.at("spacer"),begin,end,{row.at("repeat")},back);
            auto& into=source_models[{id,back}];into.insert(into.end(),models.begin(),models.end());
        }
        for(const auto& item:sequences)for(bool back:{true,false})if(!source_models.count({item.first,back}))throw std::runtime_error("missing source context");
    }
    int threads=std::getenv("MCAAT_AUDIT_THREADS")?std::stoi(std::getenv("MCAAT_AUDIT_THREADS")):4;
    if(threads<1 || threads>4)throw std::runtime_error("audit threads must be 1..4");
    omp_set_num_threads(threads);std::cout<<"Loading existing graph\n"<<std::flush;SDBG raw;raw.LoadFromFile(argv[1]);ExtensionGraph g(raw);
    std::map<std::string,Path> seeds;for(auto [id,seq]:sequences)seeds[id]=g.trace(seq);
    std::map<std::pair<std::string,bool>,Path> endpoints;
    for(auto [id,seq]:sequences)for(bool back:{true,false}) {
        endpoints[{id,back}]={};
        if(seq.size()<g.k())endpoints[{id,back}]=audit::short_seed_endpoints(g,seq,back);
        else if(!seeds[id].empty())endpoints[{id,back}]={back?seeds[id].front():seeds[id].back()};
    }
    std::ofstream exact_file(dir+"/exact_native_routes.tsv"),summary_file(dir+"/oracle_routes.tsv"),paths_file(dir+"/oracle_paths.tsv");
    exact_file<<"case_id\tseed_id\tside\treference\tseed_length\tseed_mapped\treference_seed_exact\treference_boundary_clipped\tavailable_reference_bp\tknown_reference_bp\texact_bp\tfirst_stop\tvertex\texpected_base\tavailable_edges\tplacement_count\texact_start_vertex\n";
    summary_file<<"case_id\tseed_id\tside\treference\tcheckpoint_ref_bp\tstatus\tquery_bp\tedit_cost\tidentity\tstates_discovered\tstates_popped\tlimited\texhausted\tlowest_unexplored_cost\n";
    paths_file<<"case_id\tseed_id\tside\tcheckpoint_ref_bp\tquery_bp\tsequence_outward\toperations\tvertices\tcounts\n";
    size_t done=0,totalstates=0;
    std::vector<std::map<std::string,size_t>> bound_caches(4);
    using Filtered=audit::AdapterFilteredGraph<ExtensionGraph>;
    std::vector<std::unordered_map<Node,bool>> adapter_caches(4);
    std::vector<std::map<std::pair<std::string,bool>,std::unique_ptr<Filtered>>> filtered_graphs(4);
    using Reach=audit::ReachableMembership<Filtered>;
    std::vector<std::map<std::pair<std::string,bool>,std::unique_ptr<Reach>>> reachability(4);
    std::ofstream bounds_file(dir+"/segment_bounds.tsv");bounds_file<<"case_id\tside\treference_bp\tglobal_cost_bound\tblock_costs\n";
    std::ofstream strong_file(dir+"/strong_bounds.tsv");strong_file<<"case_id\tside\treference_bp\tglobal_cost_bound\tblock_costs\n";
    std::exception_ptr error;
    std::mutex search_mutex;
    std::map<std::string,std::shared_future<audit::Search>> searches;
    #pragma omp parallel for schedule(dynamic,1)
    for(size_t case_index=0;case_index<cases.size();++case_index){try{
        const auto& r=cases[case_index];if(!filter.empty()&&!filter.count(r.at("case_id")))continue;
        auto& bound_cache=bound_caches.at(omp_get_thread_num());
        std::ostringstream exact,summary,paths,bound_report,strong_report;size_t case_states=0;
        auto sid=r.at("seed_id");
        for(bool back:{true,false}){
            auto ref=r.at(back?"left_outward":"right_outward");size_t available=ref.size();auto unknown=ref.find_first_not_of("ACGT");if(unknown!=std::string::npos)ref.resize(unknown);
            const auto& roots=endpoints.at({sid,back});
            auto& filtered=filtered_graphs.at(omp_get_thread_num())[{sid,back}];
            if(!filtered)filtered=std::make_unique<Filtered>(g,filter_adapters,roots,sequences.at(sid).size(),back,adapter_caches.at(omp_get_thread_num()));
            auto& g=*filtered;
            Reach* reachable=nullptr;
            if(use_reachability) {
                auto& ptr=reachability.at(omp_get_thread_num())[{sid,back}];
                if(!ptr)ptr=std::make_unique<Reach>(g,roots,back);
                reachable=ptr.get();
            }
            Node start=roots.empty()?missing:roots.front(),v=start;size_t n=0;
            for(auto root:roots) {
                auto current=root;size_t count=0;
                while(count<ref.size()){
                    Node next=missing;for(auto a:g.neighbors(current,back))if(g.extension_base(a,back)==ref[count]){next=a.to;break;}
                    if(next==missing)break;current=next;++count;
                }
                if(count>n){n=count;v=current;start=root;}
            }
            std::string reason=start==missing?(sequences[sid].size()<g.k()?"input_shorter_than_k":"input_seed_unmapped"):n<ref.size()?"exact_edge_absent":unknown!=std::string::npos?"reference_ambiguous_base":available<750?"reference_end":"750_exact_bases";
            exact<<r.at("case_id")<<'\t'<<sid<<'\t'<<(back?"left":"right")<<'\t'<<r.at("reference")<<'\t'<<sequences[sid].size()<<'\t'<<!roots.empty()<<'\t'<<(sequences[sid]==r.at("reference_allele"))<<'\t'<<r.at("reference_boundary_clipped")<<'\t'<<available<<'\t'<<ref.size()<<'\t'<<n<<'\t'<<reason<<'\t'<<v<<'\t'<<(n<ref.size()?ref.substr(n,1):"-")<<'\t'<<(v!=missing?choices(g,v,back):"-")<<'\t'<<roots.size()<<'\t'<<start<<'\n';
            std::set<size_t> checkpoints;for(size_t t:{24,100,250,500,750})if(t<=ref.size())checkpoints.insert(t);if(!ref.empty())checkpoints.insert(ref.size());
            if(any_length){checkpoints.clear();for(size_t t=100;t<=ref.size();++t)checkpoints.insert(t);}
            audit::SegmentBounds bounds;bounds.block=g.k();
            for(size_t i=0;!exact_only&&i+g.k()<=ref.size();i+=g.k()) {
                auto word=ref.substr(i,g.k());if(back)std::reverse(word.begin(),word.end());
                auto key=(use_reachability||filter_adapters)?sid+std::to_string(back)+word:word;
                auto cached=bound_cache.find(key);
                if(cached==bound_cache.end())cached=bound_cache.emplace(key,reachable?audit::segment_lower_bound(*reachable,word):audit::segment_lower_bound(g,word)).first;
                bounds.costs.push_back(cached->second);
            }
            bounds.prepare();if(!exact_only)bound_report<<r.at("case_id")<<'\t'<<(back?"left":"right")<<'\t'<<ref.size()<<'\t'<<bounds(0,ref.size())<<'\t'<<joined(bounds.costs)<<std::endl;
            audit::SegmentBounds strong;strong.block=g.k()+2;
            if(strong_bounds) {
                for(size_t i=0;i+strong.block<=ref.size();i+=strong.block) {
                    auto word=ref.substr(i,strong.block);if(back)std::reverse(word.begin(),word.end());
                    auto key="strong"+sid+std::to_string(back)+word;auto cached=bound_cache.find(key);
                    if(cached==bound_cache.end())cached=bound_cache.emplace(key,audit::two_edit_segment_lower_bound(g,word)).first;
                    strong.costs.push_back(cached->second);
                }
                strong.prepare();strong_report<<r.at("case_id")<<'\t'<<(back?"left":"right")<<'\t'<<ref.size()<<'\t'<<strong(0,ref.size())<<'\t'<<joined(strong.costs)<<std::endl;
            }
            auto remaining=[&](size_t pos,size_t target){return strong_bounds?std::max(bounds(pos,target),strong(pos,target)):bounds(pos,target);};
            audit::PrefixContext prefix(contexts_path?source_models.at({sid,back}):std::vector<FlankContext>{});
            audit::CycleContext cycle(g.k(),cycles_path?cycle_models.at(sid):std::vector<FlankContext>{});
            audit::Search result;if(!exact_only&&start!=missing&&!checkpoints.empty()) {
                std::string key=sid+"|"+std::to_string(back)+"|"+ref;
                std::shared_ptr<std::promise<audit::Search>> promise;std::shared_future<audit::Search> future;
                {
                    std::lock_guard<std::mutex> lock(search_mutex);
                    auto it=searches.find(key);
                    if(it==searches.end()) {
                        promise=std::make_shared<std::promise<audit::Search>>();future=promise->get_future().share();searches.emplace(key,future);
                    }else future=it->second;
                }
                if(promise)try{promise->set_value(rolling?
                    audit::align_rolling(g,start,ref,back,checkpoints,limit,remaining,roots,any_length,contexts_path?&prefix:nullptr,cycles_path?&cycle:nullptr):unitigs?
                    audit::align_unitigs(g,start,ref,back,checkpoints,limit,remaining,roots,any_length,prefer_progress,contexts_path?&prefix:nullptr,cycles_path?&cycle:nullptr):
                    audit::align(g,start,ref,back,checkpoints,limit,remaining,roots,any_length,prefer_progress,contexts_path?&prefix:nullptr,cycles_path?&cycle:nullptr));}
                    catch(...){promise->set_exception(std::current_exception());}
                result=future.get();
            }
            case_states+=result.discovered;
            if(any_length){checkpoints.clear();checkpoints.insert(result.found.empty()?ref.size():result.found.begin()->first);}
            if(checkpoints.empty())checkpoints.insert(0);
            for(auto t:checkpoints){auto it=result.found.find(t);bool found=it!=result.found.end();
                std::string status=exact_only?"approximate_search_not_requested":found?"verified_graph_path":start==missing?reason:((any_length&&ref.size()<100)||t==0)?reason:result.lowest_unexplored_cost<=t?"search_limit_unresolved":"no_path_at_90pct_anchored_identity";
                summary<<r.at("case_id")<<'\t'<<sid<<'\t'<<(back?"left":"right")<<'\t'<<r.at("reference")<<'\t'<<t<<'\t'<<status<<'\t';
                if(found){auto& w=it->second;summary<<w.sequence.size()<<'\t'<<w.edits<<'\t'<<100.*w.matches/w.columns;
                    std::vector<double> counts;for(size_t i=1;i<w.nodes.size();++i)counts.push_back(back?g.count(w.nodes[i],w.nodes[i-1]):g.count(w.nodes[i-1],w.nodes[i]));
                    paths<<r.at("case_id")<<'\t'<<sid<<'\t'<<(back?"left":"right")<<'\t'<<t<<'\t'<<w.sequence.size()<<'\t'<<w.sequence<<'\t'<<w.operations<<'\t'<<joined(w.nodes)<<'\t'<<joined(counts)<<'\n';
                }else summary<<"0\tNA\tNA";
                summary<<'\t'<<result.discovered<<'\t'<<result.popped<<'\t'<<result.limited<<'\t'<<result.exhausted<<'\t'<<result.lowest_unexplored_cost<<'\n';
            }
        }
        #pragma omp critical
        {
            exact_file<<exact.str()<<std::flush;summary_file<<summary.str()<<std::flush;paths_file<<paths.str()<<std::flush;bounds_file<<bound_report.str()<<std::flush;strong_file<<strong_report.str()<<std::flush;
            totalstates+=case_states;
            if(++done%5==0)std::cout<<"Audited "<<done<<" reference locations; "<<totalstates<<" diagnostic states\n"<<std::flush;
        }
    }catch(...){
        #pragma omp critical
        {if(!error)error=std::current_exception();}
    }
    }
    if(error)std::rethrow_exception(error);
    std::cout<<"Complete: "<<done<<" reference locations\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
