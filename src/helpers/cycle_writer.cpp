#include "cycle_writer.h"

// generate documentation with: doxygen Doxyfile

CycleWriter::CycleWriter(const std::string& filename, SDBG& sdbg) :
    sequences_output_file(filename, std::ios::out | std::ios::app),
    ids_output_file(filename + "_ids", std::ios::out | std::ios::app),
    sdbg(sdbg)
{
    if (!sequences_output_file.is_open() or !sequences_output_file.is_open()) {
        throw std::invalid_argument("Could not open file " + filename);
    }
}

// destructor
CycleWriter::~CycleWriter() {
    sequences_output_file.close();
    ids_output_file.close();
}

void CycleWriter::write_path(std::vector<int> &path) {
    
    int k = sdbg.k();

    for(int j=0; j<path.size(); j++){
        
        std::string label;
        
        uint8_t seq[k];
        uint32_t t = sdbg.GetLabel(path[j], seq);

        for (int i = k - 1; i >= 0; --i) 
            label.append(1, "ACGT"[seq[i] - 1]);

        reverse(label.begin(), label.end());
        this->sequences_output_file << label << " ";
        this->ids_output_file << path[j] << " ";
    }
    
    this->sequences_output_file << "\n";
    this->ids_output_file << "\n";
}