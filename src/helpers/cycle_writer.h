#ifndef OUTPUT_HANDLER_H
#define OUTPUT_HANDLER_H

#include <iostream>
#include <fstream>
#include <string>
#include <../src/sdbg/sdbg.h>

class CycleWriter {
public:

    CycleWriter(const std::string& filename, SDBG& sdbg);
    ~CycleWriter();
    void write_path(std::vector<int> &path);
    
private:

    SDBG& sdbg;
    std::ofstream sequences_output_file; 
    std::ofstream ids_output_file;
   
};

#endif // OUTPUT_HANDLER_H
