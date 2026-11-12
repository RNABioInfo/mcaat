// Offline experiment: extend specified observed graph alleles independently.
// The main MCAAT executable does not read this diagnostic fixture.
#include "protospacer_detection/sdbg_extension_graph.h"
#include "protospacer_detection/retained_target.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <omp.h>
using namespace mcaat::proto;
using Row=std::map<std::string,std::string>;
std::vector<std::string> split(const std::string& s){std::vector<std::string> r;std::istringstream in(s);std::string v;while(std::getline(in,v,'\t'))r.push_back(v);return r;}
std::vector<Row> table(const std::string& path){std::ifstream in(path);if(!in)throw std::runtime_error("cannot open "+path);std::string line;std::getline(in,line);auto keys=split(line);std::vector<Row> rows;while(std::getline(in,line)){auto values=split(line);if(values.size()!=keys.size())throw std::runtime_error("invalid table");Row r;for(size_t i=0;i<keys.size();++i)r[keys[i]]=values[i];rows.push_back(r);}return rows;}
int main(int argc,char** argv){try {
    if(argc!=3)throw std::runtime_error("usage: audit_variant_extensions GRAPH_PREFIX DIR");
    std::string dir=argv[2];omp_set_num_threads(2);SDBG raw;raw.LoadFromFile(argv[1]);ExtensionGraph graph(raw);
    std::map<std::string,std::vector<TargetContext>> contexts;std::map<std::string,std::set<std::string>> cycles;
    for(const auto& r:table(dir+"/contexts.tsv"))contexts[r.at("seed_id")].push_back({r.at("spacer"),r.at("repeat"),std::stoull(r.at("begin")),std::stoull(r.at("end"))});
    for(const auto& r:table(dir+"/cycles.tsv"))cycles[r.at("seed_id")].insert(r.at("cycle"));
    std::ofstream out(dir+"/isolated_extensions.tsv");out<<"seed_id\tbeam\tobserved_core\textended\tseed_begin\tseed_end\tlength\tscore\texpanded\tlimited\tleft_stop\tright_stop\tsequence\n";
    for(const auto& r:table(dir+"/cores.tsv"))for(size_t n:{16,32,64}) {
        auto id=r.at("seed_id"),sequence=r.at("observed_core");auto offset=std::stoull(r.at("query_offset")),length=std::stoull(r.at("query_length"));
        auto models=contexts.at(id);for(auto& c:models){c.begin-=offset;c.end+=sequence.size()-offset-length;}
        Options options;options.beam_width=n;options.max_flank_bp-=sequence.size()-length;
        auto result=extend_retained_target(graph,sequence,models,options,&cycles.at(id),false);
        auto& c=result.selected;out<<id<<'\t'<<n<<'\t'<<sequence<<'\t'<<result.has_extension<<'\t'<<c.seed_begin+offset<<'\t'<<c.seed_begin+offset+length<<'\t'<<c.sequence.size()<<'\t'<<c.score<<'\t'<<result.expanded<<'\t'<<result.search_limited<<'\t'<<c.left.stop<<'\t'<<c.right.stop<<'\t'<<c.sequence<<'\n';out.flush();
        std::cout<<id<<" N="<<n<<": "<<c.sequence.size()<<" bases"<<std::endl;
    }
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
