# HMM Scoring: Mathematical Framework

## Implementation: HMMER-Compatible Log-Odds Scoring

Based on the HMMER3 algorithm.

### Formula
Log-odds scoring (standard HMM approach):
```
Score = log₂ [P(sequence | HMM) / P(sequence | null)]
      = [log P(sequence | HMM) - log P(sequence | null)] / log(2)
```

### Components

**HMM Score** (numerator) - Viterbi algorithm:
- Sum of transition probabilities (M→M, M→I, etc.)
- Sum of emission probabilities (match/insert states)
- Special state transitions (N, B, E, C for local alignment)

**Null Model Score** (denominator) - background model:
- Background amino acid frequencies from HMM file (COMPO line)
- Σᵢ log P(aaᵢ | background)

**Log-Odds** - likelihood ratio test:
- Subtract null from HMM score before converting to bits
- Result is length-independent statistical measure

### Interpretation

| Score Range | Meaning |
|-------------|---------|
| > 0 bits    | More likely from HMM family than random (homolog) |
| = 0 bits    | Equally likely HMM or random |
| < 0 bits    | More likely random than HMM family |

**Bit score N** means sequence is 2^N times more likely to come from HMM family than random background.

### Example Results

**Good sequence** (362 AA, no stop codons):
- Our score: 105.259 bits
- HMMER score: 105.5 bits
- Difference: 0.23% (excellent agreement)
- Interpretation: 2^105.259 ≈ 10^31 times more likely to be homolog

**Sequence with stop codons**:
- Score: -286 bits
- Interpretation: Strongly rejected (not from HMM family)

### Stop Codon Handling

**Emission scores**: Return -INFINITY for stop codon (`*`)
- Blocks paths with stop codons in DP
- Results in very large negative scores
- Biologically correct: genes shouldn't have internal stop codons

**Null model**: Skip stop codons (use default background frequency)
- Stop codons filtered at emission level
- Should never contribute to null score

### Mathematical Properties

**Length independence**:
```
Score(seq₁ ∪ seq₂) ≈ Score(seq₁) + Score(seq₂)
```
Each position contributes independent evidence.

**Monotonicity**:
Higher score = stronger homology evidence
Can rank sequences or paths by score.

**Statistical validity**:
- Likelihood ratio test statistic
- Score > 0 provides evidence for homology
- Magnitude indicates strength of evidence

### Comparison to Alternative Approaches

**Raw log probability** (tested but rejected):
- Formula: Score = log P(seq | HMM)
- Problems:
  - Always negative (probabilities < 1)
  - More negative for longer sequences
  - No reference point for "good" vs "bad"
  - Can't compare to HMMER
  
**MegaGTA approach** (A* graph search):
- Formula: Score = log P(path | HMM) + exit_penalty
- Context: Ranking variable-length paths in de Bruijn graph
- Not applicable to our Viterbi alignment use case

### Implementation Details

**File**: src/profile.cpp, lines 540-555
**Algorithm**: HMMER3 local alignment - Plan7 architecture
**States**: N, B, M, I, D, E, C, J (J disabled for single-domain)
**Special transitions**: Loaded from HMM file (xsc table)

### Validation

✓ Matches HMMER within 0.3% on test sequences
✓ Stop codons properly filtered
✓ Scores biologically interpretable
✓ Length-independent comparisons

## References

1. Eddy SR. (2011) Accelerated Profile HMM Searches. *PLoS Comput Biol* 7(10): e1002195. https://doi.org/10.1371/journal.pcbi.1002195

2. Durbin R, Eddy SR, Krogh A, Mitchison G. (1998) *Biological Sequence Analysis: Probabilistic Models of Proteins and Nucleic Acids.* Cambridge University Press.

3. HMMER User Guide. http://eddylab.org/software/hmmer/Userguide.pdf
