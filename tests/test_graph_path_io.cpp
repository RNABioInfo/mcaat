#include "crispr_array/post_processing.h"
#include "core/array_to_nodes.h"
#include <iostream>
#include <filesystem>
#include <omp.h>
using namespace mcaat::arrays;
size_t checks=0;
void check(bool ok,const std::string& msg){++checks;if(!ok)throw std::runtime_error(msg);}
std::string contents(const std::filesystem::path& path) {
    std::ifstream in(path);return {std::istreambuf_iterator<char>(in),std::istreambuf_iterator<char>()};
}
void write(const std::filesystem::path& path,const std::string& text){std::ofstream out(path);out<<text;}
int main(int argc,char**argv) {
    try {
        if(argc!=5)throw std::runtime_error("test_graph_path_io GRAPH CYCLES OUTPUT SUPPORT");
        omp_set_num_threads(1);
        SDBG graph;graph.LoadFromFile(argv[1]);
        Settings settings;settings.sdbg=&graph;settings.output_folder=argv[3];
        settings.dna_sequence_settings.repeat_support=parse_repeat_support(argv[4]);
        std::filesystem::create_directories(settings.output_folder);
        PostProcessor processor(settings);processor.run_analysis(argv[2]);
        auto arrays=BuildFilteredArraysFromDir(settings.output_folder,graph);
        check(!arrays.empty(),"fixture has arrays");
        SdbgRepeatGraph g(graph);
        size_t spacers=0,nodes=0;
        for(auto& a:arrays) {
            std::string repeat;
            for(auto e:a.repeat_path) {
                check(graph.IsValidEdge(e),"downstream repeat edges valid");
                auto s=g.label(e);
                if(repeat.empty()) repeat=s;else repeat+=s.back();
                ++nodes;
            }
            check(repeat==a.repeat_sequence,"downstream repeat path spells observed repeat");
            for(auto& p:a.spacer_node_paths) {
                ++spacers;nodes+=p.size();
                for(size_t i=0;i<p.size();++i) {
                    check(graph.IsValidEdge(p[i]),"downstream spacer edges valid");
                    if(i)check(graph.Forward(p[i-1])==uint64_t(graph.GetLastIndex(p[i])),"spacer traversal continuity");
                }
            }
        }
        auto dir=std::filesystem::path(settings.output_folder);
        const auto path_table=contents(dir/"CRISPR_graph_paths.tsv");
        auto before=contents(dir/"CRISPR_Arrays_1.txt");
        write(dir/"CRISPR_Arrays_99.txt","stale");write(dir/"user_notes.txt","keep");
        processor.run_analysis(argv[2]);
        check(!std::filesystem::exists(dir/"CRISPR_Arrays_99.txt"),"stale numbered output cleared");
        check(contents(dir/"user_notes.txt")=="keep","unrelated output preserved");
        check(contents(dir/"CRISPR_graph_paths.tsv")==path_table && contents(dir/"CRISPR_Arrays_1.txt")==before,"rerun deterministic");
        auto corrupt=path_table;
        size_t pos=corrupt.find('\n')+1;
        pos=corrupt.find('\t',pos)+1;pos=corrupt.find('\t',pos)+1;
        corrupt[pos]=corrupt[pos]=='A'?'C':'A';
        write(dir/"bad_paths.tsv",corrupt);
        bool failed=false;try{BuildFilteredArraysFromGraphPaths((dir/"bad_paths.tsv").string(),graph);}catch(const std::runtime_error&){failed=true;}
        check(failed,"mismatched graph/sequence sidecar rejected");
        std::filesystem::remove(dir/"bad_paths.tsv");std::filesystem::remove(dir/"user_notes.txt");
        // Edge-slot validity is distinct from k-mer vertex identity.
        bool tested=false;
        for(Node e=0;e<graph.size() && !tested;++e) {
            if(!graph.IsValidEdge(e)) continue;
            Node vertex=graph.GetLastIndex(e);
            if(vertex==e || !graph.IsValidEdge(vertex)) continue;
            auto label=g.label(vertex);auto target=graph.Forward(e);
            if(graph.Forward(vertex)==target)continue;
            graph.SetInvalidEdge(vertex);
            check(g.find(label)==vertex && g.edge(vertex,target)==e,"invalid last edge slot does not erase source vertex");
            graph.SetValidEdge(vertex);tested=true;
        }
        check(tested,"branched source-vertex fixture exercised");
        Settings empty=settings;empty.output_folder=(dir/"empty").string();
        std::filesystem::create_directories(empty.output_folder);
        write(std::filesystem::path(empty.output_folder)/"cycles.txt",contents(argv[2]));
        PostProcessor(empty).run_analysis();
        check(BuildFilteredArraysFromDir(empty.output_folder,graph).empty(),"empty current run cannot reuse stale cycles.txt");
        auto config=dir/"test-settings.txt";write(config,"repeat_support=0.60\n");
        Settings configured;configured.LoadFromFile(config.string());
        check(configured.dna_sequence_settings.repeat_support==.6,"settings-file fraction applied");
        write(config,"repeat_support=nan\n");failed=false;
        try{configured.LoadFromFile(config.string());}catch(const std::invalid_argument&){failed=true;}
        check(failed,"invalid settings-file fraction rejected");std::filesystem::remove(config);
        std::cout<<checks<<" SDBG/output checks passed: "<<arrays.size()<<" arrays, "<<spacers<<" spacers, "<<nodes<<" downstream path nodes\n";
    }catch(const std::exception&e){std::cerr<<"FAILED: "<<e.what()<<'\n';return 1;}
}
