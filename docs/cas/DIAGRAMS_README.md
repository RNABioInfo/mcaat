# Documentation & Diagrams

## Documentation

### HMM_VITERBI_BEAM_SEARCH.md
Technical documentation: architecture, algorithms, scoring, validation, usage.

## Diagrams

Paste into online tools to generate visualizations.

### GraphViz (*.dot) → https://dreampuf.github.io/GraphvizOnline/

**algorithm_flowchart.dot**: Beam search main loop
**viterbi_flowchart.dot**: Viterbi DP details
**scoring_pipeline.dot**: Score transformations

Render locally:
```bash
dot -Tpng <file>.dot -o <output>.png
```

### PlantUML (*.puml) → https://www.planttext.com/

**sequence_diagram.puml**: Component interactions
**class_diagram.puml**: Code structure

Render locally:
```bash
plantuml <file>.puml
```

## Quick Reference

| File | Type | Purpose |
|------|------|---------|
| algorithm_flowchart.dot | GraphViz | Main loop |
| viterbi_flowchart.dot | GraphViz | Scoring |
| scoring_pipeline.dot | GraphViz | Transformations |
| sequence_diagram.puml | PlantUML | System flow |
| class_diagram.puml | PlantUML | Architecture |

## Score Transformations

```
HMMER3:    +2.90958  (-log P)
Parse:     -2.90958  (log P)
Viterbi:   -1413.38  (raw)
Null:      -1513.59
Log-Odds:  +100.21
Bits:      +100.215
```

Visualized in **scoring_pipeline.dot**.

## Validation

- CasGeneDetector: 105.259 bits
- HMMER: 105.5 bits
- Difference: 0.23%
- 704 profiles: 2.16% average

## For Publications

Export as PDF/SVG:
```bash
dot -Tpdf algorithm_flowchart.dot -o figure1.pdf
dot -Tsvg viterbi_flowchart.dot -o figure2.svg
```

## Resources

- HMMER3: http://eddylab.org/software/hmmer/Userguide.pdf
- GraphViz: https://graphviz.org/documentation/
- PlantUML: https://plantuml.com/

**Implementation files**:
- Code: [src/profile.cpp](src/profile.cpp), [src/cas_gene_detector.cpp](src/cas_gene_detector.cpp)
- Headers: [include/profile.h](include/profile.h), [include/cas_gene_detector.h](include/cas_gene_detector.h)
- Tests: [src/test_cas_gene_detector.cpp](src/test_cas_gene_detector.cpp)
