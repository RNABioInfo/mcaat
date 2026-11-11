#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <cstdint>
template<class S>void audit_side(const char*,bool,const S&,const std::string&);
template<class S,class A>void audit_arc(const char*,bool,const S&,const A&,double,const std::string&);
template<class C>void audit_candidate(const char*,const C&,size_t);
// The exact production context loader and extension implementation.
#include "../src/protospacer_detection/extension_module.cpp"
using namespace mcaat::proto;
std::ofstream events,candidates,sides;std::string current_id;bool logged_back=false;
std::map<std::pair<std::string,bool>,std::set<size_t>> watched_lengths;
uint64_t hash_path(const Path& p){uint64_t h=14695981039346656037ULL;for(auto v:p){h^=v;h*=1099511628211ULL;}return h;}
template<class S>void event_prefix(const char* kind,bool back,const S& s){logged_back=back;events<<current_id<<'\t'<<kind<<'\t'<<(back?"left":"right")<<'\t'<<s.nodes.size()-1<<'\t'<<hash_path(s.nodes)<<'\t'<<s.nodes.back()<<'\t';}
template<class C>auto broad_median(const C& c,int)->decltype(c.sampled_recent,double()){return median(c.sampled_recent);}
template<class C>double broad_median(const C& c,long){return c.local_mean();}
template<class S>void event_suffix(const S& s,double penalty,const std::string& why){auto& p=s.coverage.profile;events<<s.coverage.local_mean()<<'\t'<<penalty<<'\t'<<p.loss()<<'\t'<<s.beam_profile.loss()<<'\t'<<why<<'\t'<<s.context<<'\t'<<p.samples<<'\t'<<broad_median(s.coverage,0)<<'\t';
 auto it=watched_lengths.find({current_id,logged_back});uint64_t h=14695981039346656037ULL;
 if(it!=watched_lengths.end())for(size_t i=0;i<s.nodes.size();++i){h^=s.nodes[i];h*=1099511628211ULL;if(it->second.count(i))events<<i<<':'<<h<<',';}
 events<<'\n';}
template<class S>void audit_side(const char* kind,bool back,const S& s,const std::string& why){if(s.nodes.empty())return;event_prefix(kind,back,s);events<<"-\t-\t-\t";event_suffix(s,0,why);if(std::string(kind)=="side_retained")sides<<current_id<<'\t'<<(back?"left":"right")<<'\t'<<s.stop<<'\t'<<s.coverage.profile.loss()<<'\t'<<path_text(s.nodes)<<'\n';}
template<class S,class A>void audit_arc(const char* kind,bool back,const S& s,const A& a,double penalty,const std::string& why){if(s.nodes.empty())return;event_prefix(kind,back,s);events<<a.edge<<'\t'<<a.to<<'\t'<<a.count<<'\t';event_suffix(s,penalty,why);}
template<class C>void audit_candidate(const char* stage,const C& c,size_t rank){auto st=stability(c.left.coverage,c.right.coverage);candidates<<current_id<<'\t'<<stage<<'\t'<<rank<<'\t'<<c.score<<'\t'<<c.sequence<<'\t'<<c.seed_begin<<'\t'<<c.seed_end<<'\t'<<c.left.stop<<'\t'<<c.right.stop<<'\t'<<st.within<<'\t'<<st.jumps<<'\t'<<st.between<<'\t'<<st.weak<<'\t'<<st.support<<'\n';}
int main(int argc,char** argv){try{
 if(argc!=4)throw std::runtime_error("usage: audit_losses GRAPH_PREFIX ARRAYS_DIR AUDIT_DIR");std::string dir=argv[3];Options options;
 events.open(dir+"/production_events.tsv");candidates.open(dir+"/production_candidates.tsv");sides.open(dir+"/production_sides.tsv");
 events<<std::setprecision(17);candidates<<std::setprecision(17);
 events<<"seed_id\tevent\tside\tdepth\tpath_hash\tvertex\tedge\ttarget\tcount\tlocal_median\tpenalty\tside_loss\tbeam_loss\treason\tcontext\tsamples\tbroad_median\tprefix_hashes\n";
 candidates<<"seed_id\tstage\trank\tscore\tsequence\tseed_begin\tseed_end\tleft_stop\tright_stop\twithin\tjumps\tbetween\tweak\tsupport\n";
 sides<<"seed_id\tside\tstop\tloss\tvertices\n";
 omp_set_num_threads(4);std::cout<<"Loading graph\n"<<std::flush;SDBG raw;raw.LoadFromFile(argv[1]);ExtensionGraph graph(raw);auto queries=load_queries(argv[2],graph.k(),options);
 std::ofstream selected(dir+"/production_selected.tsv"),contexts(dir+"/production_contexts.tsv"),cycles(dir+"/production_source_cycles.tsv");
 selected<<"seed_id\tstatus\treason\tsequence\tseed_begin\tseed_end\tscore\texpanded\n";contexts<<"seed_id\tspacer\trepeat\tbegin\tend\n";cycles<<"seed_id\tcycle\n";
 if(fs::exists(dir+"/oracle_paths.tsv"))table(dir+"/oracle_paths.tsv",[&](const Fields& r){watched_lengths[{r.at("seed_id"),r.at("side")=="left"}].insert(std::stoull(r.at("query_bp")));});
 if(fs::exists(dir+"/watch_selected.tsv")) {
  std::ofstream out(dir+"/watched_selected_paths.tsv");
  out<<"case_id\tseed_id\tside\tquery_bp\tsequence_outward\tvertices\n";
  table(dir+"/watch_selected.tsv",[&](const Fields& r){
   if(r.at("sequence").empty())return;
   auto full=graph.trace(r.at("sequence"));if(full.empty())throw std::runtime_error("baseline prediction no longer maps to graph");
   auto lo=std::stoull(r.at("seed_begin")),hi=std::stoull(r.at("seed_end"));
   for(bool back:{true,false}) {
    Path path=back?Path(full.begin(),full.begin()+lo+1):Path(full.begin()+hi-graph.k(),full.end());
    auto seq=back?r.at("sequence").substr(0,lo):r.at("sequence").substr(hi);
    if(back){std::reverse(path.begin(),path.end());std::reverse(seq.begin(),seq.end());}
    watched_lengths[{r.at("seed_id"),back}].insert(seq.size());
    out<<"baseline_"<<r.at("seed_id")<<'\t'<<r.at("seed_id")<<'\t'<<(back?"left":"right")<<'\t'<<seq.size()<<'\t'<<seq<<'\t'<<path_text(path)<<'\n';
   }
  });
 }
 size_t done=0;
 table(dir+"/targets.tsv",[&](const Fields& row){
  current_id=row.at("seed_id");auto seq=row.at("sequence");std::vector<TargetContext> models;std::set<std::string> source_cycles;
  std::set<std::tuple<std::string,std::string,size_t,size_t>> unique;
  for(const auto& q:queries)for(bool flip:{false,true}){
   auto full=flip?reverse_complement(q.sequence):q.sequence;
   for(size_t pos=full.find(seq);pos!=std::string::npos;pos=full.find(seq,pos+1)){
    for(auto c:q.source_cycles)source_cycles.insert(flip?reverse_complement(c):c);
    for(const auto& s:q.sources){bool source_flip=s.reverse!=flip;auto repeat=source_flip?reverse_complement(s.repeat):s.repeat;
     if(unique.emplace(full,repeat,pos,pos+seq.size()).second)models.push_back({full,repeat,pos,pos+seq.size()});
    }
   }
  }
  for(const auto& m:models)contexts<<current_id<<'\t'<<m.spacer<<'\t'<<m.repeat<<'\t'<<m.begin<<'\t'<<m.end<<'\n';
  for(const auto& c:source_cycles)cycles<<current_id<<'\t'<<c<<'\n';
  auto e=extend_retained_target(graph,seq,models,options,&source_cycles);
  selected<<current_id<<'\t'<<e.status<<'\t'<<e.reason<<'\t'<<(e.has_extension?e.selected.sequence:"")<<'\t'<<e.selected.seed_begin<<'\t'<<e.selected.seed_end<<'\t'<<e.selected.score<<'\t'<<e.expanded<<'\n';
  if(++done%10==0){selected.flush();events.flush();candidates.flush();std::cout<<"Recorded decisions for "<<done<<"/81 inputs\n"<<std::flush;}
 });
 std::cout<<"Complete: "<<done<<" inputs\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
