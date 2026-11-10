// Read-only per-k-mer trace for diagnosing an existing graph. No reference is
// used by the extension algorithm; this utility only reports graph evidence.
#include "protospacer_detection/sdbg_extension_graph.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <omp.h>
using namespace mcaat::proto;
int main(int argc,char** argv) {
    if(argc!=3){std::cerr<<"usage: trace_protospacer_graph GRAPH_PREFIX QUERIES.fasta\n";return 2;}
    try {
        omp_set_num_threads(4);SDBG raw;raw.LoadFromFile(argv[1]);ExtensionGraph graph(raw);
        std::ifstream in(argv[2]);if(!in)throw std::runtime_error("cannot open queries");
        std::cout<<"query_id\toffset\tkmer\tvertex\tincoming\toutgoing\tprevious_edge_count\n";
        auto arcs=[&](Node node,bool reverse) {
            std::ostringstream out;bool first=true;
            if(node!=missing)for(auto a:graph.neighbors(node,reverse)) {
                if(!first)out<<',';first=false;
                out<<graph.extension_base(a,reverse)<<':'<<a.count<<':'<<a.to;
            }
            return out.str();
        };
        auto trace=[&](const std::string& id,const std::string& seq) {
            Node previous=missing;
            for(size_t i=0;i+graph.k()<=seq.size();++i) {
                Node node=missing;double count=-1;
                if(previous!=missing)for(auto a:graph.neighbors(previous,false))
                    if(graph.extension_base(a,false)==seq[i+graph.k()-1]){node=a.to;count=a.count;break;}
                if(node==missing)node=graph.find(seq.substr(i,graph.k()));
                std::cout<<id<<'\t'<<i<<'\t'<<seq.substr(i,graph.k())<<'\t';
                if(node!=missing)std::cout<<node;
                std::cout<<'\t'<<arcs(node,true)<<'\t'<<arcs(node,false)<<'\t';
                if(count>=0)std::cout<<count;
                std::cout<<'\n';previous=node;
            }
        };
        std::string line,id,sequence;
        while(std::getline(in,line)) {
            if(line.empty())continue;
            if(line[0]=='>'){if(!id.empty())trace(id,sequence);id=line.substr(1,line.find(' ')-1);sequence.clear();}
            else {if(line.back()=='\r')line.pop_back();sequence+=line;}
        }
        if(!id.empty())trace(id,sequence);
        if(!in.eof())throw std::runtime_error("failed reading queries");
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
