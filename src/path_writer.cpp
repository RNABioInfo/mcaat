#include "path_writer.h"

PathWriter::PathWriter(string mode, SDBG& sdbg, vector<uint64_t> path, string genome_id, string type, vector<bool>& visited) 
    : sdbg(sdbg), genome_id(genome_id), type(type), visited(visited) {
    this->CollectPathsIntoStringStreams(path);
}

PathWriter::~PathWriter() {
    // Destructor
    
}

string PathWriter::FetchNodeLabel(size_t node) {
    std::string label;            
    uint8_t seq[sdbg.k()];
    uint32_t t = sdbg.GetLabel(node, seq);
    for (int i = sdbg.k() - 1; i >= 0; --i) label.append(1, "ACGT"[seq[i] - 1]);
    reverse(label.begin(), label.end());
    return label;
}

string PathWriter::FetchFirstNodeLabel(size_t node) {
    std::string label;            
    uint8_t seq[sdbg.k()];
    uint32_t t = sdbg.GetLabel(node, seq);
    for (int i = sdbg.k() - 1; i >= 0; --i) label.append(1, "ACGT"[seq[i] - 1]);
    reverse(label.begin(), label.end());
    return label;
}

//@brief Write the path nodes ids, their labels and the multiplicities to the respective files
void PathWriter::CollectPathsIntoStringStreams(vector<uint64_t> path) {
    // Reserve space and add the start node to the path
    path.reserve(path.size() + 1);
    path.push_back(path[0]);

    // Create string streams for path, id, and multiplicity
    stringstream ss_path, ss_id, ss_multiplicity;

    // Iterate over the path and populate the string streams
    for (uint64_t node : path) {
        ss_path << FetchNodeLabel(node) << " ";
        ss_id << node << " ";
        ss_multiplicity << this->sdbg.EdgeMultiplicity(node) << " ";
        this->visited[node] = true;
    }

    // Add new lines to the string streams
    ss_path << endl;
    ss_id << endl;
    ss_multiplicity << endl;
    ofstream cycle_report_file;
    ofstream path_report_file;
    ofstream multiplicity_report_file;

    string cycle_file_name = this->folder_path + this->genome_id + "/cycles/str_paths.txt";
    string id_file_name = this->folder_path + this->genome_id + "/cycles/id_paths.txt";
    string multiplicity_file_name = this->folder_path + this->genome_id + "/cycles/multiplicity_distribution.txt";
    
    multiplicity_report_file.open(multiplicity_file_name, std::ios_base::app);
    cycle_report_file.open(cycle_file_name, std::ios_base::app);
    path_report_file.open(id_file_name, std::ios_base::app);
    // Write the string streams to the respective files
    cycle_report_file << ss_path.str();
    path_report_file << ss_id.str();
    multiplicity_report_file << ss_multiplicity.str();

    // Close the files
    cycle_report_file.close();
    path_report_file.close();
    multiplicity_report_file.close();
}