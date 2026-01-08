# mCAAT - Metagenomic CRISPR Array Analysis Tool

## New Components (January 2026)

### AminoAcidator
Beam search-based amino acid translator for de Bruijn graph traversal.

**Files:**
- `include/amino_acidator.h`
- `src/amino_acidator.cpp`
- `src/main_test_acidator.cpp`

**Compile & Run:**
```bash
g++ -std=c++17 -O3 -march=native -fopenmp -g \
  -I./include -I./libs/megahit/src -I./libs/kseqpp/include \
  -DXXH_INLINE_ALL -ftemplate-depth=3000 -Wall -Wno-unused-function \
  -fprefetch-loop-arrays -funroll-loops \
  src/main_test_acidator.cpp src/amino_acidator.cpp \
  libs/megahit/src/sdbg/sdbg_meta.cpp \
  libs/megahit/src/sdbg/sdbg_raw_content.cpp \
  libs/megahit/src/sdbg/sdbg_writer.cpp \
  libs/megahit/src/utils/options_description.cpp \
  -lz -lpthread -o test_acidator

./test_acidator
```

**Features:**
- Beam search traversal through succinct de Bruijn graph
- Converts DNA triplets to amino acids using standard genetic code
- Avoids revisiting nodes (loop detection)
- Returns paths with amino acid sequences, node paths, and scores

---

### Profile (HMMER3 HMM Parser)
Reads and stores HMMER3 profile HMM files for sequence alignment scoring.

**Files:**
- `include/profile.h`
- `src/profile.cpp`
- `src/test_profile.cpp`

**Compile & Run:**
```bash
g++ -std=c++17 -O3 -g -I./include \
  src/test_profile.cpp src/profile.cpp \
  -o test_profile

./test_profile
# or with custom HMM file:
./test_profile path/to/profile.hmm
```

**Features:**
- Parses HMMER3 format HMM files
- Stores match/insert emission scores (20 amino acids)
- Stores transition probabilities (M->M, M->I, M->D, I->M, I->I, D->M, D->D)
- Extracts consensus sequence
- Query emission and transition scores by position

**Example HMM:** `hmm_test.hmm` (COG1518, 328 states)

---

### Simulated Read Generator
Python script to generate FASTQ reads with flanking sequences.

**File:** `build/generate_sim_reads.py`

**Usage:**
```bash
python3 generate_sim_reads.py input.fasta [left_flank_bp] [right_flank_bp] [coverage] [read_length]

# Default: 1M bp flanks, 30x coverage, 150 bp reads
python3 generate_sim_reads.py the_sequence.fasta 1000000 1000000 30 150
```

**Output:** `sim_reads.fastq` (uppercase DNA sequences, Phred 30-40 quality)

---

## Key Data Structures

**AminoAcidPathInfo:**
- `amino_acids`: vector of amino acid strings
- `node_path`: vector of node IDs traversed
- `scores`: per-node scores
- `total_score`: cumulative path score
- `dna_sequence`: accumulated DNA sequence

**ProfileState:**
- `position`: match state position
- `match_emissions`: 20 amino acid emission scores
- `insert_emissions`: 20 amino acid insertion scores
- `transitions`: 7 transition probabilities
- `consensus`: consensus amino acid character

---

## Testing

All test files compile independently with minimal dependencies (see commands above).
