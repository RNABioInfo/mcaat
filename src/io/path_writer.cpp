#include "io/path_writer.h"
#ifdef DEBUG
PathWriter::PathWriter(string mode, SDBG& sdbg, vector<uint64_t> path, string genome_id, string type) 
    : sdbg(sdbg), genome_id(genome_id), type(type) {
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




//@brief Write the path nodes ids, their labels and the multiplicities to the respective files
void PathWriter::CollectPathsIntoStringStreams(vector<uint64_t> path) {
    // Create the folder and get its name

    // Create string streams for path and id
    stringstream ss_path, ss_id;

    // Iterate over the path and populate the string streams
    for (uint64_t node : path) {
        ss_path << FetchNodeLabel(node) << " ";
        ss_id << node << " ";
    }

    // Add new lines to the string streams
    ss_path << endl;
    ss_id << endl;

    // Open files in the created folder
    ofstream cycles_file(this->genome_id + "/cycles.txt", std::ios_base::app);
    ofstream labels_file(this->genome_id + "/labels.txt", std::ios_base::app);

    // Write to files
    cycles_file << ss_path.str();
    labels_file << ss_id.str();

    // Close files
    cycles_file.close();
    labels_file.close();
}


#endif