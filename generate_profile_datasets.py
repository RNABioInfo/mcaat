#!/usr/bin/env python3
"""
Generate test datasets from HMM profiles.
For each profile:
1. Extract consensus sequence
2. Translate AA to DNA (random codon choice)
3. Generate simulated reads with flanking regions
"""

import os
import sys
import random
import subprocess
from pathlib import Path

# Standard genetic code (codon → amino acid)
CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

# Reverse lookup: AA → list of codons
AA_TO_CODONS = {}
for codon, aa in CODON_TABLE.items():
    if aa != '*':  # Skip stop codons
        if aa not in AA_TO_CODONS:
            AA_TO_CODONS[aa] = []
        AA_TO_CODONS[aa].append(codon)


def parse_hmm_consensus(hmm_file):
    """Extract consensus sequence using hmmemit from HMMER"""
    try:
        # Use hmmemit to generate consensus sequence
        result = subprocess.run(
            ['hmmemit', '-c', hmm_file],
            capture_output=True,
            text=True,
            check=True
        )
        
        # Parse FASTA output
        lines = result.stdout.strip().split('\n')
        sequence = ''
        for line in lines:
            if not line.startswith('>'):
                sequence += line.strip()
        
        return sequence
        
    except subprocess.CalledProcessError as e:
        print(f"ERROR running hmmemit: {e}")
        print(f"stderr: {e.stderr}")
        return ''
    except FileNotFoundError:
        print("ERROR: hmmemit not found. Please install HMMER3.")
        print("  Ubuntu/Debian: sudo apt install hmmer")
        print("  macOS: brew install hmmer")
        return ''


def get_profile_name(hmm_file):
    """Get profile name from filename"""
    return Path(hmm_file).stem


def translate_aa_to_dna(aa_sequence):
    """Translate amino acid sequence to DNA using random codon choices"""
    dna = []
    for aa in aa_sequence:
        if aa in AA_TO_CODONS:
            # Randomly choose one codon for this amino acid
            codon = random.choice(AA_TO_CODONS[aa])
            dna.append(codon)
        else:
            print(f"Warning: Unknown amino acid '{aa}', skipping")
    return ''.join(dna)


def generate_random_dna(length):
    """Generate random DNA sequence"""
    bases = ['A', 'C', 'G', 'T']
    return ''.join(random.choice(bases) for _ in range(length))


def generate_reads_for_sequence(dna_sequence, output_dir, kmer_name, coverage=30, read_length=150):
    """Generate simulated reads using generate_sim_reads.py approach"""
    
    # Add flanking regions
    left_flank = generate_random_dna(100000)
    right_flank = generate_random_dna(100000)
    full_sequence = left_flank + dna_sequence + right_flank
    
    # Write temporary FASTA file
    temp_fasta = output_dir / "temp_sequence.fasta"
    with open(temp_fasta, 'w') as f:
        f.write(f">sequence\n")
        # Write in 80-character lines
        for i in range(0, len(full_sequence), 80):
            f.write(full_sequence[i:i+80] + '\n')
    
    # Generate reads
    fastq_file = output_dir / f"{kmer_name}.fastq"
    seq_length = len(full_sequence)
    num_reads = int((coverage * seq_length) / read_length)
    
    print(f"  Generating {num_reads} reads (coverage={coverage}x)...")
    
    with open(fastq_file, 'w') as out:
        for i in range(num_reads):
            # Random start position
            if seq_length > read_length:
                start = random.randint(0, seq_length - read_length)
            else:
                start = 0
            
            # Extract read
            read_seq = full_sequence[start:start + read_length]
            
            # Random quality scores (Phred 30-40)
            quality = ''.join(chr(random.randint(63, 73)) for _ in range(len(read_seq)))
            
            # Write FASTQ format
            out.write(f"@read_{i+1}_{start}\n")
            out.write(f"{read_seq}\n")
            out.write("+\n")
            out.write(f"{quality}\n")
    
    # Clean up temporary file
    temp_fasta.unlink()
    
    actual_coverage = (num_reads * read_length) / seq_length
    print(f"  Reads written to {fastq_file.name}")
    print(f"  Actual coverage: {actual_coverage:.2f}x")


def process_profile(hmm_file, output_base_dir, coverage=30, read_length=150):
    """Process a single HMM profile"""
    
    print(f"\n{'='*60}")
    print(f"Processing: {hmm_file}")
    print('='*60)
    
    # Get profile name
    profile_name = get_profile_name(hmm_file)
    print(f"Profile name: {profile_name}")
    
    # Extract consensus sequence
    consensus_aa = parse_hmm_consensus(hmm_file)
    if not consensus_aa:
        print(f"ERROR: Could not extract consensus from {hmm_file}")
        return
    
    print(f"Consensus AA sequence: {len(consensus_aa)} residues")
    print(f"First 50 AA: {consensus_aa[:50]}")
    
    # Translate to DNA
    dna_sequence = translate_aa_to_dna(consensus_aa)
    print(f"DNA sequence: {len(dna_sequence)} nucleotides")
    print(f"First 69 nt: {dna_sequence[:69]}")
    
    # Get first 23 nucleotides for filename
    kmer_23 = dna_sequence[:23] if len(dna_sequence) >= 23 else dna_sequence
    print(f"23-mer for filename: {kmer_23}")
    
    # Create output directory
    output_dir = output_base_dir / profile_name
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Save AA FASTA
    aa_fasta = output_dir / f"{profile_name}_consensus.fasta"
    with open(aa_fasta, 'w') as f:
        f.write(f">{profile_name}_consensus\n")
        for i in range(0, len(consensus_aa), 80):
            f.write(consensus_aa[i:i+80] + '\n')
    print(f"AA sequence saved to: {aa_fasta.name}")
    
    # Save DNA FASTA
    dna_fasta = output_dir / f"{profile_name}_dna.fasta"
    with open(dna_fasta, 'w') as f:
        f.write(f">{profile_name}_dna\n")
        for i in range(0, len(dna_sequence), 80):
            f.write(dna_sequence[i:i+80] + '\n')
    print(f"DNA sequence saved to: {dna_fasta.name}")
    
    # Generate reads
    generate_reads_for_sequence(dna_sequence, output_dir, kmer_23, coverage, read_length)
    
    print(f"✓ Completed: {profile_name}")


def main():
    if len(sys.argv) < 2:
        print("Usage: python generate_profile_datasets.py <profiles_dir> [output_dir] [coverage] [read_length]")
        print("\nGenerates test datasets from HMM profiles:")
        print("  1. Extract consensus sequence from each profile")
        print("  2. Translate AA → DNA (random codon choices)")
        print("  3. Generate simulated reads with 100kb flanking regions")
        print("\nArguments:")
        print("  profiles_dir : Directory containing .hmm files")
        print("  output_dir   : Output directory (default: profile_datasets)")
        print("  coverage     : Read coverage depth (default: 30)")
        print("  read_length  : Read length in bp (default: 150)")
        sys.exit(1)
    
    profiles_dir = Path(sys.argv[1])
    output_dir = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("profile_datasets")
    coverage = int(sys.argv[3]) if len(sys.argv) > 3 else 30
    read_length = int(sys.argv[4]) if len(sys.argv) > 4 else 150
    
    # Find all HMM files
    hmm_files = list(profiles_dir.glob("*.hmm"))
    if not hmm_files:
        print(f"ERROR: No .hmm files found in {profiles_dir}")
        sys.exit(1)
    
    print(f"Found {len(hmm_files)} HMM profiles")
    print(f"Output directory: {output_dir}")
    print(f"Coverage: {coverage}x")
    print(f"Read length: {read_length} bp")
    
    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Process each profile
    for hmm_file in sorted(hmm_files):
        try:
            process_profile(hmm_file, output_dir, coverage, read_length)
        except Exception as e:
            print(f"ERROR processing {hmm_file}: {e}")
            import traceback
            traceback.print_exc()
    
    print(f"\n{'='*60}")
    print(f"✓ All profiles processed!")
    print(f"Output directory: {output_dir}")
    print('='*60)


if __name__ == "__main__":
    main()
