#include "protospacer_detection/sdbg_extension_graph.h"
#include "protospacer_detection/target_placement.h"
#include "short_seed_placements.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <omp.h>
using namespace mcaat::proto;
int main(int argc,char** argv){try {
    if(argc!=3)throw std::runtime_error("usage: audit_target_placements GRAPH_PREFIX DIR");
    std::string dir=argv[2];omp_set_num_threads(2);SDBG raw;raw.LoadFromFile(argv[1]);ExtensionGraph graph(raw);
    std::ifstream input(dir+"/targets.tsv");std::ofstream output(dir+"/placements.tsv");
    if(!input||!output)throw std::runtime_error("missing target table/output directory");
    output<<"seed_id\tquery\tobserved_core\tquery_offset\tquery_begin\tquery_end\tidentity\tcoverage\twork\tlimited\n";
    std::string line;std::getline(input,line);
    while(std::getline(input,line)) {
        std::istringstream row(line);std::string id,query;std::getline(row,id,'\t');std::getline(row,query,'\t');
        if(query.size()<graph.k()) {
            bool capped=false;auto indexed=graph.find_suffix(query,100000,capped);
            auto brute=audit::short_seed_endpoints(graph,query,false);std::sort(brute.begin(),brute.end());
            if(capped||indexed!=brute)throw std::runtime_error("suffix index disagrees with exhaustive padding: "+id);
        }
        auto found=place_target(graph,query,0);
        for(const auto& p:found.values)output<<id<<'\t'<<query<<'\t'<<p.sequence<<'\t'<<p.query_offset<<'\t'<<p.query_begin<<'\t'<<p.query_end<<'\t'<<p.identity<<'\t'<<p.coverage<<'\t'<<found.work<<'\t'<<found.limited<<'\n';
        output.flush();std::cout<<id<<": "<<found.values.size()<<" placements, "<<found.work<<" work, limited="<<found.limited<<std::endl;
    }
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
