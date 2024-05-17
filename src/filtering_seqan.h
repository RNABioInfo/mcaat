#include <seqan/align.h>
#include "edlib.h"
using std::istringstream;
using namespace seqan;
using namespace std;



typedef String<char> TSequence;
typedef vector<TSequence> TCycleStream;
typedef Align<TSequence, ArrayGaps> TAlign;
typedef Row<TAlign>::Type TRow;
typedef Iterator<TRow>::Type TRowIterator;
using IteratorFikrators = Iterator<TRow>::Type;

TCycleStream AlignCycles(TCycleStream cycles, int longest_sequence_pos) {
    TCycleStream keep_cycles_list;
    TSequence longest_cycle = cycles[longest_sequence_pos];


    for (int i = 0; i < cycles.size(); i++) {
        TAlign align;
        resize(rows(align), 2);
        assignSource(row(align, 0), longest_cycle);
        assignSource(row(align, 1), cycles[i]);

        int score = globalAlignment(align, Score<int, Simple>(1, -1, -1));
        /*
        TRow & row1 = row(align, 0);
        TRow & row2 = row(align, 1);
        // Initialize the row iterators.
        TRowIterator itRow1 = begin(row1);
        TRowIterator itEndRow1 = end(row1);
        TRowIterator itRow2 = begin(row2);
        int gapCount = 0;
        for (; itRow1 != itEndRow1; ++itRow1, ++itRow2)
        {
            if (isGap(itRow1))
            {
                gapCount += countGaps(itRow1);  // Count the number of consecutive gaps from the current position in row1.
                itRow1 += countGaps(itRow1);    // Jump to next position to check for gaps.
                itRow2 += countGaps(itRow1);    // Jump to next position to check for gaps.
            }
            if (isGap(itRow2))
            {
                gapCount += countGaps(itRow2);  // Count the number of consecutive gaps from the current position in row2.
                itRow1 += countGaps(itRow2);    // Jump to next position to check for gaps.
                itRow2 += countGaps(itRow2);    // Jump to next position to check for gaps.
            }
        }
        int editDistance = abs(score) - gapCount;
        */
        auto length_of_short = length(cycles[i]);

        double percentage = score/double(length_of_short);
        // Print the result.
        if (percentage > 0.10) keep_cycles_list.push_back(cycles[i]);
        
    }
    if (keep_cycles_list.size() != 0) keep_cycles_list.push_back(longest_cycle);
    return keep_cycles_list;
}

void WriteCycles(TCycleStream cycles) {
    //write cycles to file
    ofstream log_file;
    string log_filename = "log.txt";
    log_file.open(log_filename,std::ios_base::app);
    for (int i = 0; i < cycles.size(); i++) {
        log_file << cycles[i] << endl;
    }
    log_file.close();

}
void TestEdlib(){
    EdlibAlignResult result = edlibAlign("GTTTCCATTCAATTAGTTTCCCCAGCGAGTGGGGACACTGAACCAGGAAAACCAATGGGCTACTTTAAAGATTCTGAGTTTCCATTCAATTAGTTTCCCC", 100, "GTTTCCATTCAATTAGTTTCCCCAGCGAGTGGGGACTCCCTCCTTCAGGAGGAGATCCTTCGGGCGTTCGCCTGCCTCCCAGAGTTTCCATTCAATTAGTTTCCCC", 122, edlibDefaultAlignConfig());
    if (result.status == EDLIB_STATUS_OK) {
        printf("edit_distance = %d\n", result.editDistance);
    }
    edlibFreeAlignResult(result);
}
int FilterOut(string filename) {
    std::ifstream file(filename);
    bool initial_line = true;
    TCycleStream cycles;
    vector<size_t> lengths;
    for (std::string line; std::getline(file, line);) {
        if(initial_line) {
            initial_line = false;
            continue;
        }
        if(line.size() < 40){
            if(cycles.size()==1){
                WriteCycles(cycles);
                cycles.clear();
                lengths.clear();
                continue;
            }
            
            int longest_sequence_pos = max_element(lengths.begin(), lengths.end()) - lengths.begin();
            TCycleStream kept_cycles = AlignCycles(cycles,longest_sequence_pos);

            if(kept_cycles.size()>0) WriteCycles(kept_cycles);
            cycles.clear();
            lengths.clear();
            continue;
        }
        else{
            lengths.push_back(line.size());
            cycles.push_back(line);
        }
    }
    return 0;


}
