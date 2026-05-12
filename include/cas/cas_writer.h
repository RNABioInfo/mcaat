#ifndef CAS_WRITER_H
#define CAS_WRITER_H

#include "cas/cas_workflow.h"
#include <string>
#include <vector>
#include <fstream>

/**
 * CasWriter — writes Cas gene detection results to numbered output files.
 *
 * Output format per system:
 *
 *   ================================================================================
 *   Array  N   subtype  I-E   span  4832 bp   direction  downstream
 *   Repeat     ATCGATCGATCGATCGATCGATCG
 *
 *     cas8e    score  0.82   length  144 aa   offset   127 bp
 *     MKLQELIAKNDEIQLAQRFVEDLKSKDQHPQLPSSFLDDLFKAKEAGKPIVEQVFKQMKQILQMKLQ
 *     ELIAKNDEIQLAQRFVEDLKSKDQHPQLPSSFLDDLFKAKEAGKPIVEQVFKQMKQILQMKLAAKQM
 *
 *     cas7     score  0.71   length  298 aa   offset  1843 bp
 *     MRILAQVDPQKFST...
 *
 * Files are split at max_lines_per_file and named CAS_Systems_1.txt, CAS_Systems_2.txt, ...
 */
class CasWriter {
public:
    static constexpr int DEFAULT_MAX_LINES = 20000;
    static constexpr int AA_LINE_WIDTH     = 70;

    explicit CasWriter(const std::string& output_dir,
                       int max_lines_per_file = DEFAULT_MAX_LINES);
    ~CasWriter();

    /**
     * Write all results to numbered files in output_dir.
     * results = vector of (repeat_sequence, [upstream_cassette, downstream_cassette])
     * Returns total number of systems (cassettes with at least one gene) written.
     */
    int Write(const std::vector<std::pair<std::string,
                  std::vector<CasCassette>>>& results);

private:
    void WriteFileHeader();
    void WriteSystemBlock(const std::string& repeat_seq,
                          const CasCassette& cassette,
                          int array_number);
    void WriteAASequence(const std::string& aa);

    // Emits text and tracks line count; rolls to next file when limit reached.
    void Emit(const std::string& line);
    void EmitBlank();
    void OpenNextFile();
    void CloseFile();
    std::string CurrentFilePath() const;

    std::string output_dir_;
    int max_lines_per_file_;
    int file_index_     = 1;
    int line_count_     = 0;   // lines in current file (including header)
    int systems_written_= 0;
    int total_arrays_   = 0;
    int total_genes_    = 0;
    std::ofstream out_;
    std::string generated_ts_;  // filled once at construction
};

#endif // CAS_WRITER_H
