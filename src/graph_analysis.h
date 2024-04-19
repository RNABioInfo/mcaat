#ifndef GRAPH_ANALYSIS_H
#define GRAPH_ANALYSIS_H


#include <../src/sdbg/sdbg.h>
#include <vector>
#include <map>
#include <set>
#include <array>
#include <stack>
#include <omp.h>

using namespace std;

class GraphAnalysis {
    
    private:
        SDBG& sdbg;
    public:
        GraphAnalysis(SDBG& sdbg): sdbg(sdbg) {}
        
        void Analysis(){
            for (uint64_t i = 0; i < sdbg.size(); i++)
               if (i<=25)
                   cout << "i :" <<i << " reverse complement: "<< sdbg.EdgeReverseComplement(i) << endl;
            
        }

};

#endif 