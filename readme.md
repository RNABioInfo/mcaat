# MCAAT-CAS: Cas Gene Detector

HMM-guided beam search for detecting Cas genes in de Bruijn graphs using HMMER3 profiles.

## Algorithm

**Beam Search + Viterbi Alignment:**
- Traverses de Bruijn graph, translates DNA → amino acids
- Scores paths with HMMER3 Viterbi alignment at each step
- HMMER-compatible log-odds scoring with null model correction
- Beam width: 10 paths, max depth: HMM_length × 3

## Compile

```bash
# Cas gene detector
g++ -std=c++17 -O3 -march=native -fopenmp -I./include -I./libs/megahit/src -I./libs/kseqpp/include -o test_cas_gene_detector src/test_cas_gene_detector.cpp src/cas_gene_detector.cpp src/profile.cpp libs/megahit/src/sdbg/sdbg_meta.cpp libs/megahit/src/sdbg/sdbg_raw_content.cpp libs/megahit/src/sdbg/sdbg_writer.cpp libs/megahit/src/utils/options_description.cpp

# Graph builder
g++ -std=c++17 -O3 -march=native -fopenmp -I./include -I./libs/megahit/src -I./libs/kseqpp/include -o build_graph_cli src/build_graph_cli.cpp src/sdbg_build.cpp libs/megahit/src/sdbg/sdbg_meta.cpp libs/megahit/src/sdbg/sdbg_raw_content.cpp libs/megahit/src/sdbg/sdbg_writer.cpp libs/megahit/src/sorting/kmer_counter.cpp libs/megahit/src/sorting/read_to_sdbg_s1.cpp libs/megahit/src/sorting/read_to_sdbg_s2.cpp libs/megahit/src/sorting/seq_to_sdbg.cpp libs/megahit/src/utils/options_description.cpp libs/megahit/src/sorting/base_engine.cpp libs/megahit/src/sorting/kmsort_selector.cpp libs/megahit/src/sequence/io/fastx_reader.cpp libs/megahit/src/sequence/io/sequence_lib.cpp libs/megahit/src/sequence/io/paired_fastx_reader.cpp -lz

# K-mer finder
g++ -std=c++17 -O3 -march=native -fopenmp -I./include -I./libs/megahit/src -I./libs/kseqpp/include -o find_kmer_id src/find_kmer_id.cpp libs/megahit/src/sdbg/sdbg_meta.cpp libs/megahit/src/sdbg/sdbg_raw_content.cpp libs/megahit/src/sdbg/sdbg_writer.cpp libs/megahit/src/utils/options_description.cpp
```

## Usage

```bash
# 1. Build graph from FASTQ
./build_graph_cli reads.fastq output_dir

# 2. Find k-mer node ID
./find_kmer_id output_dir/graph/graph ATCGATCGATCGATCGATCGATC

# 3. Run Cas gene detector
./test_cas_gene_detector output_dir/graph/graph profile.hmm <node_id>
```

## Validation (704 Cas Profiles)

```bash
./validate.sh
```

**Results:**
- Profiles tested: 704
- Average score difference vs HMMER: **2.16%**
- Scoring method: HMMER3 log-odds with null model correction

## Files

**Core:**
- `src/cas_gene_detector.cpp` - Beam search implementation
- `src/profile.cpp` - HMMER3 parser & Viterbi alignment
- `src/test_cas_gene_detector.cpp` - Main entry point

**Tools:**
- `src/build_graph_cli.cpp` - Graph builder
- `src/find_kmer_id.cpp` - K-mer lookup
- `validate.sh` - Validation script

**HMM Scores:**
- **HMMER3 File Format**: Stores positive values representing `-log(P)`
- **Internal Storage**: Negated to log probabilities (all negative)
- **Viterbi DP**: Pure log-probability addition (sums remain negative)
- **Null Model Correction**: Subtracts background model score
  - `log_odds = viterbi_score - null_score`
  - Can be positive for good matches!
- **Bit Score Conversion**: `bit_score = log_odds / log(2)`
  - **Positive scores** indicate good HMM matches
  - **Higher is better** (e.g., 100.2 bits > 50 bits)
- **Score Range Examples**:
  - Poor match: -50 to 0 bits
  - Weak match: 0 to 30 bits
  - Good match: 30 to 80 bits
  - Strong match: 80+ bits (our result: 100.215 bits!)

**All Transition Types Used:**
- M→M, M→I, M→D (from Match state)
- I→M, I→I (from Insert state)
- D→M, D→D (from Delete state)

---

## Documentation

### Technical Documentation
- **HMM_VITERBI_BEAM_SEARCH.md** - Complete technical documentation
  - Architecture and algorithm descriptions
  - Detailed scoring scheme explanations
  - HMMER validation and comparison
  - Usage examples and performance analysis

### Visualizations
See **DIAGRAMS_README.md** for all visualization files:

**GraphViz Diagrams** (paste into https://dreampuf.github.io/GraphvizOnline/):
- `algorithm_flowchart.dot` - Main beam search algorithm
- `viterbi_flowchart.dot` - Detailed Viterbi DP with null model
- `scoring_pipeline.dot` - Complete scoring transformation pipeline

**PlantUML Diagrams** (paste into https://www.planttext.com/):
- `sequence_diagram.puml` - Component interaction sequence
- `class_diagram.puml` - Complete class architecture

---

## Testing

All test files compile independently with minimal dependencies (see commands above).
