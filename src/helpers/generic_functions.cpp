#include "generic_functions.h"
using namespace std;


/// @brief Helper Function gets all outgoing edges of a node
/// @param uint64_t node
string FetchNodeLabel(SDBG& succinct_de_bruijn_graph,auto node){
    std::string label;            
    uint8_t seq[succinct_de_bruijn_graph.k()];
    uint32_t t = succinct_de_bruijn_graph.GetLabel(node, seq);
    for (int i = succinct_de_bruijn_graph.k() - 1; i >= 0; --i) label.append(1, "ACGT"[seq[i] - 1]);
    reverse(label.begin(), label.end());
    return label;
}

/// @brief Helper Function gets all outgoing edges of a node
/// @param uint64_t node
/// @return Outgoing edges of a node
set<uint64_t> GetOutgoingNodes(SDBG& succinct_de_bruijn_graph, uint64_t node){
    if (succinct_de_bruijn_graph.EdgeOutdegree(node) == 0)
        return {};
            
    uint64_t *outgoings = new uint64_t[succinct_de_bruijn_graph.EdgeOutdegree(node)];
    succinct_de_bruijn_graph.OutgoingEdges(node,outgoings);
    set<uint64_t> outgoings_set(outgoings, outgoings + succinct_de_bruijn_graph.EdgeOutdegree(node));
    return outgoings_set;
};
        
/// @brief Helper Function gets all incoming edges of a node
/// @param uint64_t node
/// @return Incoming edges of a node
set<uint64_t> GetIncomingNodes(SDBG& succinct_de_bruijn_graph,uint64_t node){
    if (succinct_de_bruijn_graph.EdgeIndegree(node) == 0)
        return {};
            
    uint64_t *incomings = new uint64_t[succinct_de_bruijn_graph.EdgeIndegree(node)];
    succinct_de_bruijn_graph.IncomingEdges(node,incomings);
    set<uint64_t> incomings_set(incomings, incomings + succinct_de_bruijn_graph.EdgeIndegree(node));
    return incomings_set;
};

/*
void GenericWrite(vector<uint64_t> path, uint64_t startnode, string type, string genome_id, 
    SDBG& succinct_de_bruijn_graph, char file_mode,string folder_path="proof_of_concept/data/")
{
    path.push_back(startnode);
    ofstream generic_file;
    switch(file_mode){
        case('p'):
            generic_file.open(folder_path+genome_id+"/cycles"+type+"/str_paths.txt",std::ios_base::app);
            for(int j=0; j<path.size(); j++)
                generic_file << FetchNodeLabel(succinct_de_bruijn_graph,path[j]) << " ";
            
        case ('i'):
            generic_file.open(folder_path+genome_id+"/cycles"+type+"/id_paths.txt",std::ios_base::app);
            for(int j=0; j<path.size(); j++)
                generic_file << path[j] << " ";
                break;
            
        case ('m'):
            generic_file.open(folder_path+genome_id+"/cycles"+type+"/edge_multiplicities.txt",std::ios_base::app);
            for(int j=0; j<path.size(); j++)
                generic_file << succinct_de_bruijn_graph.EdgeMultiplicity(path[j]) << " ";
            break;
        }
    
    generic_file << "\n";
    generic_file.close();
    return;
}
void WritePath(vector<uint64_t> path, uint64_t startnode, string modes, string type, string genome_id, SDBG& succinct_de_bruijn_graph,string folder_path="proof_of_concept/data/")
{            
    uint8_t N = modes.length();
    for (int i = 0; i < N; i++)
        GenericWrite(path, startnode, type, genome_id, succinct_de_bruijn_graph, modes[i], folder_path);
    return;
};

void WriteBlankLine(string type,string genome_id,string folder_path="proof_of_concept/data/"){
    
    ofstream id_file;
    id_file.open(folder_path+genome_id+"/cycles"+type+"/id_paths.txt",std::ios_base::app);
    id_file << "\n";
    id_file.close();

    ofstream str_file;
    str_file.open(folder_path+genome_id+"/cycles"+type+"/str_paths.txt",std::ios_base::app);
    str_file << "\n";
    str_file.close();
    return;
}*/
