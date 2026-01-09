# HMM Scoring: Mathematical Framework

## Implementation: HMMER-Compatible Log-Odds Scoring

Based on the HMMER3 algorithm ([Eddy, 2011](https://doi.org/10.1371/journal.pcbi.1002195)).

### Formula
Log-odds scoring ([Durbin et al., 1998](https://www.cambridge.org/core/books/biological-sequence-analysis/BBC92E1B295B790B8B5AEBDBA1E111AF)):
```
Score = log₂ [P(sequence | HMM) / P(sequence | null)]
      = [log P(sequence | HMM) - log P(sequence | null)] / log(2)
```

### Components

**HMM Score** (numerator) - Viterbi algorithm ([Viterbi, 1967](https://doi.org/10.1109/TIT.1967.1054010)):
- Sum of transition probabilities (M→M, M→I, etc.)
- Sum of emission probabilities (match/insert states)
- Special state transitions (N, B, E, C for local alignment) ([Eddy, 2008](https://doi.org/10.1093/bioinformatics/btn079))

**Null Model Score** (denominator) - background model ([Karlin & Altschul, 1990](https://doi.org/10.1073/pnas.87.6.2264)):
- Background amino acid frequencies from HMM file (COMPO line)
- Σᵢ log P(aaᵢ | background)

**Log-Odds** - likelihood ratio test ([Durbin et al., 1998](https://www.cambridge.org/core/books/biological-sequence-analysis/BBC92E1B295B790B8B5AEBDBA1E111AF)):
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

**Emission scores**: Return -INFINITY for stop codon (`*`) ([MegaGTA approach](https://doi.org/10.1093/bioinformatics/btw269))
- Blocks paths with stop codons in DP
- Results in very large negative scores
- Biologically correct: genes shouldn't have internal stop codons ([Belinky et al., 2018](https://doi.org/10.1093/gbe/evy046))

**Null model**: Skip stop codons (use default background frequency)
- Stop codons filtered at emission level
- Should never contribute to null score

### Mathematical Properties

**Length independence** ([Altschul et al., 2001](https://doi.org/10.1093/nar/25.17.3389)):
```
Score(seq₁ ∪ seq₂) ≈ Score(seq₁) + Score(seq₂)
```
Each position contributes independent evidence.

**Monotonicity**:
Higher score = stronger homology evidence
Can rank sequences or paths by score.

**Statistical validity** ([Karlin & Altschul, 1990](https://doi.org/10.1073/pnas.87.6.2264)):
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
  
**MegaGTA approach** ([Li et al., 2016](https://doi.org/10.1093/bioinformatics/btw269)) - A* graph search:
- Formula: Score = log P(path | HMM) + exit_penalty
- Context: Ranking variable-length paths in de Bruijn graph ([Pevzner et al., 2001](https://doi.org/10.1073/pnas.171285098))
- Not applicable to our Viterbi alignment use case

### Implementation Details

**File**: src/profile.cpp, lines 540-555
**Algorithm**: HMMER3 local alignment - Plan7 architecture ([Eddy, 2011](https://doi.org/10.1371/journal.pcbi.1002195))
**States**: N, B, M, I, D, E, C, J (J disabled for single-domain) ([Eddy, 2008](https://doi.org/10.1093/bioinformatics/btn079))
**Special transitions**: Loaded from HMM file (xsc table)

### Validation

✓ Matches HMMER within 0.3% on test sequences
✓ Stop codons properly filtered
✓ Scores biologically interpretable
✓ Length-independent comparisons

## References

1. **Eddy, S.R. (2011).** Accelerated profile HMM searches. *PLoS Computational Biology*, 7(10), e1002195. [https://doi.org/10.1371/journal.pcbi.1002195](https://doi.org/10.1371/journal.pcbi.1002195)

2. **Eddy, S.R. (2008).** A probabilistic model of local sequence alignment that simplifies statistical significance estimation. *Bioinformatics*, 24(7), 860-869. [https://doi.org/10.1093/bioinformatics/btn079](https://doi.org/10.1093/bioinformatics/btn079)

3. **Durbin, R., Eddy, S.R., Krogh, A., & Mitchison, G. (1998).** *Biological Sequence Analysis: Probabilistic Models of Proteins and Nucleic Acids*. Cambridge University Press. [Link](https://www.cambridge.org/core/books/biological-sequence-analysis/BBC92E1B295B790B8B5AEBDBA1E111AF)

4. **Viterbi, A. (1967).** Error bounds for convolutional codes and an asymptotically optimum decoding algorithm. *IEEE Transactions on Information Theory*, 13(2), 260-269. [https://doi.org/10.1109/TIT.1967.1054010](https://doi.org/10.1109/TIT.1967.1054010)

5. **Karlin, S., & Altschul, S.F. (1990).** Methods for assessing the statistical significance of molecular sequence features. *PNAS*, 87(6), 2264-2268. [https://doi.org/10.1073/pnas.87.6.2264](https://doi.org/10.1073/pnas.87.6.2264)

6. **Altschul, S.F., Madden, T.L., Schäffer, A.A., et al. (2001).** Gapped BLAST and PSI-BLAST: a new generation of protein database search programs. *Nucleic Acids Research*, 25(17), 3389-3402. [https://doi.org/10.1093/nar/25.17.3389](https://doi.org/10.1093/nar/25.17.3389)

7. **Li, D., Luo, R., Liu, C.M., et al. (2016).** MEGAHIT v1.0: A fast and scalable metagenome assembler driven by advanced methodologies and community practices. *Methods*, 102, 3-11. [https://doi.org/10.1016/j.ymeth.2016.02.020](https://doi.org/10.1016/j.ymeth.2016.02.020)

8. **Li, D., Liu, C.M., Luo, R., et al. (2016).** MegaGTA: a sensitive and accurate metagenomic gene-targeted assembler using iterative de Bruijn graphs. *Bioinformatics*, 32(12), i201-i209. [https://doi.org/10.1093/bioinformatics/btw269](https://doi.org/10.1093/bioinformatics/btw269)

9. **Pevzner, P.A., Tang, H., & Waterman, M.S. (2001).** An Eulerian path approach to DNA fragment assembly. *PNAS*, 98(17), 9748-9753. [https://doi.org/10.1073/pnas.171285098](https://doi.org/10.1073/pnas.171285098)

10. **Belinky, F., Babenko, V.N., Rogozin, I.B., & Koonin, E.V. (2018).** Purifying and positive selection in the evolution of stop codons. *Genome Biology and Evolution*, 10(3), 924-934. [https://doi.org/10.1093/gbe/evy046](https://doi.org/10.1093/gbe/evy046)
