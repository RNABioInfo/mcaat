# HMM Scoring: Mathematical Framework

## HMMER-Compatible Log-Odds Scoring

### Formula
```
Score = log₂ [P(sequence | HMM) / P(sequence | null)]
      = [log P(sequence | HMM) - log P(sequence | null)] / log(2)
```

### Components

**HMM Score** - Viterbi algorithm:
- Transition probabilities: M→M, M→I, M→D, I→M, I→I, D→M, D→D
- Emission probabilities: match/insert states
- Special states: N, B, E, C (Plan7 architecture)

**Null Model** - background:
- Background AA frequencies from HMM COMPO line
- Σᵢ log P(aaᵢ | background)

**Log-Odds** - likelihood ratio:
- HMM score - null score
- Length-independent statistic

### Interpretation

| Score | Meaning |
|-------|---------|
| > 0 bits | HMM family (homolog) |
| = 0 bits | Equally likely |
| < 0 bits | Random background |

**N bits** = 2^N times more likely from HMM than random.

### Validation

| Metric | CasGeneDetector | HMMER | Difference |
|--------|----------------|-------|------------|
| Score (bits) | 105.259 | 105.5 | 0.23% |
| Sequence | 362 AA | 362 AA | - |
| Coverage | 328/328 | 296/328 | Global vs local |

**704 profiles tested**: 2.16% average difference

### Stop Codon Handling

**Emission**: Return -∞ for stop codon (`*`)
- Blocks paths with stop codons
- Biologically correct (no internal stops)

**Null model**: Skip stop codons (emission filter)

### Mathematical Properties

**Length independence**:
```
Score(seq₁ ∪ seq₂) ≈ Score(seq₁) + Score(seq₂)
```

**Monotonicity**: Higher score = stronger homology

**Statistical validity**: Likelihood ratio test

### Implementation

**File**: [src/profile.cpp](src/profile.cpp) (lines 540-555)
**Algorithm**: HMMER3 Plan7 local alignment
**States**: N, B, M, I, D, E, C (J disabled)

## References

1. Eddy SR. (2011) Accelerated Profile HMM Searches. *PLoS Comput Biol* 7(10): e1002195.
2. Durbin R, et al. (1998) *Biological Sequence Analysis.* Cambridge University Press.
3. HMMER User Guide. http://eddylab.org/software/hmmer/Userguide.pdf
