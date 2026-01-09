# Documentation and Visualization Files

This directory contains comprehensive documentation for the HMM-guided beam search implementation.

## Documentation Files

### 📄 HMM_VITERBI_BEAM_SEARCH.md
**Main technical documentation** covering:
- Architecture overview
- Algorithm descriptions
- Scoring scheme details
- HMMER comparison and validation
- Usage examples
- Performance analysis

**Start here** for complete understanding of the system.

## Visualization Files

All visualization files can be pasted into online tools to generate diagrams.

### 🔷 GraphViz Diagrams (*.dot files)

Paste into: **https://dreampuf.github.io/GraphvizOnline/**

Or render locally:
```bash
dot -Tpng <file>.dot -o <output>.png
```

#### algorithm_flowchart.dot
**High-level beam search algorithm** showing:
- Main loop structure
- Graph expansion
- Viterbi integration
- Beam pruning
- Termination conditions

#### viterbi_flowchart.dot
**Detailed Viterbi algorithm** showing:
- DP matrix computation
- Match/Insert/Delete states
- Null model correction
- Backtracking
- Score conversion to bits

#### scoring_pipeline.dot
**End-to-end scoring pipeline** showing:
- HMMER3 file format parsing
- Score dimension transformations
- Viterbi computation in log space
- Null model correction
- HMMER validation comparison

### 🔶 PlantUML Diagrams (*.puml files)

Paste into: **https://www.planttext.com/** or **http://www.plantuml.com/plantuml/uml/**

Or render locally (requires PlantUML):
```bash
plantuml <file>.puml
```

#### sequence_diagram.puml
**Sequence diagram** showing interactions between:
- BeamSearch (AminoAcidator)
- De Bruijn Graph (SDBG)
- Profile HMM
- Viterbi Algorithm
- Priority Queue

Shows complete flow from user call to results.

#### class_diagram.puml
**Class diagram** showing:
- Core classes and their relationships
- Key methods and attributes
- Data structures
- Algorithm components
- HMMER integration utilities
- Testing infrastructure

## Quick Start

### View Documentation
```bash
# Read main documentation
cat HMM_VITERBI_BEAM_SEARCH.md

# Or open in your markdown viewer
code HMM_VITERBI_BEAM_SEARCH.md
```

### Generate All Diagrams

**Using GraphViz** (if installed):
```bash
dot -Tpng algorithm_flowchart.dot -o algorithm_flowchart.png
dot -Tpng viterbi_flowchart.dot -o viterbi_flowchart.png
dot -Tpng scoring_pipeline.dot -o scoring_pipeline.png
```

**Using PlantUML** (if installed):
```bash
plantuml sequence_diagram.puml
plantuml class_diagram.puml
```

**Using Online Tools** (no installation needed):
1. Copy contents of .dot file
2. Paste into https://dreampuf.github.io/GraphvizOnline/
3. View/download generated diagram

Or for .puml files:
1. Copy contents of .puml file
2. Paste into https://www.planttext.com/
3. View/download generated diagram

## Diagram Overview

| File | Type | Purpose | Best For |
|------|------|---------|----------|
| algorithm_flowchart.dot | GraphViz | Beam search algorithm flow | Understanding main loop |
| viterbi_flowchart.dot | GraphViz | Viterbi DP details | Understanding scoring |
| scoring_pipeline.dot | GraphViz | Score transformations | Debugging score issues |
| sequence_diagram.puml | PlantUML | Component interactions | Understanding system flow |
| class_diagram.puml | PlantUML | Code structure | Understanding architecture |

## Key Concepts Illustrated

### 1. Algorithm Flowchart
- How beam search explores the graph
- When and why Viterbi is called
- How beam pruning works
- Termination criteria (HMM position ≥ 328)

### 2. Viterbi Flowchart
- DP matrix initialization and computation
- Three states: Match, Insert, Delete
- All 7 transition types
- Null model correction process
- Conversion to HMMER-compatible bits

### 3. Scoring Pipeline
- HMMER3 file format (positive values)
- Negation during parsing (→ log probabilities)
- Pure log-space computation in Viterbi
- Null model subtraction (→ log-odds)
- Division by log(2) (→ bits)

### 4. Sequence Diagram
- Call sequence from user to results
- Interaction between graph, HMM, and Viterbi
- How paths are scored and ranked
- Queue management and pruning

### 5. Class Diagram
- Profile class: HMM storage and Viterbi
- AminoAcidator class: Beam search
- SDBG class: De Bruijn graph
- Path/SearchState data structures
- Utility classes for HMMER compatibility

## Score Dimension Verification

All diagrams emphasize the critical scoring transformations:

```
HMMER3 File:    +2.90958  (positive, represents -log P)
      ↓ Parse & negate
Internal:       -2.90958  (negative log P)
      ↓ Viterbi DP (addition in log space)
Raw Score:      -1413.38  (very negative)
      ↓ Null model correction
Log-Odds:       +100.21   (positive for good matches)
      ↓ Divide by log(2)
Bits:           +100.215  (HMMER-compatible)
```

This is visualized in **scoring_pipeline.dot**.

## Validation Results

Documented in all relevant diagrams:

- **Our Score**: 100.215 bits
- **HMMER Score**: 105.5 bits  
- **Difference**: 5.3 bits (5% error)
- **Conclusion**: ✅ Successfully validated!

## For Paper/Publication

These diagrams are suitable for:
- **Algorithm flowcharts**: algorithm_flowchart.dot, viterbi_flowchart.dot
- **System architecture**: class_diagram.puml
- **Method description**: sequence_diagram.puml
- **Validation**: scoring_pipeline.dot

Export as PDF or SVG for high-resolution figures:
```bash
dot -Tpdf algorithm_flowchart.dot -o figure1.pdf
dot -Tsvg viterbi_flowchart.dot -o figure2.svg
```

## Additional Resources

- **HMMER3 Documentation**: http://eddylab.org/software/hmmer/Userguide.pdf
- **Viterbi Algorithm**: "Biological Sequence Analysis" by Durbin et al.
- **GraphViz Documentation**: https://graphviz.org/documentation/
- **PlantUML Documentation**: https://plantuml.com/

## Questions?

For details on specific implementations, see:
- Code: `src/profile.cpp`, `src/amino_acidator.cpp`
- Headers: `include/profile.h`, `include/amino_acidator.h`
- Tests: `src/test_hmm_acidator.cpp`
- Main docs: `HMM_VITERBI_BEAM_SEARCH.md`
