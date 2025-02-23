#ifndef path_writer_h
#define path_writer_h

#include <iostream>
#include <fstream>
#include <string>
#include <sdbg/sdbg.h>
#include <vector>
#include <sstream>
#include <algorithm>
using namespace std;

class PathWriter{
    private:
        SDBG& sdbg;
        const string folder_path = "somefolder";
        string genome_id;
        string type;
        vector<bool>& visited;
    public:
        PathWriter(string mode, SDBG& sdbg, vector<uint64_t> path,string genome_id, string type,  vector<bool>& visited);
        ~PathWriter();
        string FetchNodeLabel(size_t node);
        string FetchFirstNodeLabel(size_t node);
        void CollectPathsIntoStringStreams(vector<uint64_t> path);
};

#endif
