#ifndef path_writer_h
#define path_writer_h

#include <iostream>
#include <fstream>
#include <string>
#include <../src/sdbg/sdbg.h>
#include "helpers/generic_functions.cpp"

using namespace std;

class PathWriter{
    private:
        SDBG& succinct_de_bruijn_graph;
        const string folder_path = "proof_of_concept/data/";
        string genome_id;
        string type;
    public:
        PathWriter(string mode, SDBG& succinct_de_bruijn_graph, vector<uint64_t> path,string genome_id, string type): 
            succinct_de_bruijn_graph(succinct_de_bruijn_graph), genome_id(genome_id), type(type){
            
            if(mode=="p"){
                WritePathString(path, path[0]);
            }
            if (mode=="i"){
                WritePathIDs(path, path[0]);
            }
            if (mode=="pi"){
                WritePathString(path, path[0]);
                WritePathIDs(path, path[0]);
            }
            if (mode=="im"){
                WriteEdgeMultiplicities(path, path[0]);
                WritePathIDs(path, path[0]);
            }
            if (mode=="a"){
                WritePathString(path, path[0]);
                WritePathIDs(path, path[0]);
                WriteEdgeMultiplicities(path, path[0]);
            }
           
        }
         ~PathWriter(){};
        
        void WritePathString(vector<uint64_t> path, uint64_t startnode)
        {            
            ofstream cycle_report_file;
            string string_filename = this->folder_path+this->genome_id+"/cycles/str_paths.txt";
            if (type=="fasta")
                string_filename = this->folder_path+this->genome_id+"/cycles_genome/str_paths.txt";
            cycle_report_file.open(string_filename,std::ios_base::app);
            path.push_back(startnode);
            
            for(int j=0; j<path.size(); j++)
                cycle_report_file << FetchNodeLabel(this->succinct_de_bruijn_graph,path[j]) << " ";
            
            cycle_report_file << "\n";
            cycle_report_file.close();
    };
    
    void WritePathIDs(vector<uint64_t> path, uint64_t startnode)
    {            
        ofstream path_report_file;
        string id_filename = this->folder_path+this->genome_id+"/cycles/id_paths.txt";
        if (type=="fasta")
                id_filename = this->folder_path+this->genome_id+"/cycles_genome/id_paths.txt";
        path_report_file.open(id_filename,std::ios_base::app);
        path.push_back(startnode);

        for(int j=0; j<path.size(); j++)
            path_report_file << path[j] << " ";
        
        path_report_file << "\n";
        path_report_file.close();
    };
    
    void WriteEdgeMultiplicities(vector<uint64_t> path, uint64_t startnode)
    {            
        ofstream path_report_file;
        string id_filename = this->folder_path+this->genome_id+"/cycles/multiplicity_distribution.txt";
        
        path_report_file.open(id_filename,std::ios_base::app);
        path.push_back(startnode);

        for(int j=0; j<path.size(); j++)
            path_report_file << this->succinct_de_bruijn_graph.EdgeMultiplicity(path[j]) << " ";
        
        path_report_file << "\n";
        path_report_file.close();
    };
};

#endif