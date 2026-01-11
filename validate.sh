#!/bin/bash

echo "Profile,MCAAT_Score,HMMER_Score,Difference_Percent" > validation_results.txt

for HMM in profiles/*.hmm; do
    PROFILE=$(basename "$HMM" .hmm)
    DATASET="profile_datasets/$PROFILE"
    
    if [ ! -d "$DATASET" ]; then
        continue
    fi
    
    FASTQ=$(find "$DATASET" -name "*.fastq" -type f | head -n 1)
    if [ -z "$FASTQ" ]; then
        continue
    fi
    
    KMER=$(basename "$FASTQ" .fastq)
    AA_FILE="$DATASET/${PROFILE}_consensus.fasta"
    
    if [ ! -f "$AA_FILE" ]; then
        continue
    fi
    
    echo "=========================================="
    echo "Profile: $PROFILE"
    echo "=========================================="
    
    # Step 1: Build graph
    ./build_graph_cli "$FASTQ" "$DATASET" > /dev/null 2>&1
    
    # Step 2: Find k-mer node
    NODE=$(./find_kmer_id "$DATASET/graph/graph" "$KMER" 2>&1 | grep "Node ID:" | sed 's/.*Node ID: //' | awk '{print $1}')
    
    if [ -z "$NODE" ]; then
        echo "ERROR: K-mer not found"
        rm -rf "$DATASET/graph" "$DATASET/cycles"
        continue
    fi
    
    # Step 3: Beam search
    MCAAT=$(./test_cas_gene_detector "$DATASET/graph/graph" "$HMM" "$NODE" 2>&1 | grep "Viterbi score:" | head -n 1 | awk '{print $3}')
    
    if [ -z "$MCAAT" ]; then
        echo "ERROR: Beam search failed"
        rm -rf "$DATASET/graph" "$DATASET/cycles"
        continue
    fi
    
    # Step 4: HMMER
    HMMER=$(hmmsearch --max "$HMM" "$AA_FILE" 2>&1 | grep -E "^\s+[0-9]" | grep "${PROFILE}_consensus" | awk '{print $2}')
    
    if [ -z "$HMMER" ]; then
        echo "ERROR: HMMER failed"
        rm -rf "$DATASET/graph" "$DATASET/cycles"
        continue
    fi
    
    # Step 5: Calculate difference
    DIFF=$(awk -v m="$MCAAT" -v h="$HMMER" 'BEGIN {printf "%.2f", ((m-h)/h)*100}')
    DIFF_ABS=$(echo "$DIFF" | tr -d '-')
    
    echo "MCAAT:  $MCAAT bits"
    echo "HMMER:  $HMMER bits"
    echo "Diff:   $DIFF_ABS%"
    echo ""
    
    echo "$PROFILE,$MCAAT,$HMMER,$DIFF_ABS" >> validation_results.txt
    
    # Step 6: Clean up
    rm -rf "$DATASET/graph" "$DATASET/cycles"
done

echo "=========================================="
echo "DONE - Results in validation_results.txt"
echo "=========================================="
