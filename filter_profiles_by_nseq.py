#!/usr/bin/env python3
"""
Filter HMM profiles by NSEQ (number of sequences used to build the profile).
Removes profiles where NSEQ < threshold.

Usage:
    python filter_profiles_by_nseq.py profiles_dir/ output_dir/ --min-nseq 5
    python filter_profiles_by_nseq.py profiles_dir/ output_dir/ --min-nseq 10 --report
"""

import os
import sys
import argparse
import shutil
from pathlib import Path


def parse_hmm_header(filepath):
    """Parse HMM file header and extract metadata."""
    metadata = {
        'name': None,
        'acc': None,
        'desc': None,
        'leng': None,
        'nseq': None,
        'effn': None,
    }
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            
            if line.startswith('HMMER'):
                continue
            if line.startswith('NAME'):
                metadata['name'] = line.split(None, 1)[1] if len(line.split(None, 1)) > 1 else None
            elif line.startswith('ACC'):
                metadata['acc'] = line.split(None, 1)[1] if len(line.split(None, 1)) > 1 else None
            elif line.startswith('DESC'):
                metadata['desc'] = line.split(None, 1)[1] if len(line.split(None, 1)) > 1 else None
            elif line.startswith('LENG'):
                try:
                    metadata['leng'] = int(line.split()[1])
                except (IndexError, ValueError):
                    pass
            elif line.startswith('NSEQ'):
                try:
                    metadata['nseq'] = int(line.split()[1])
                except (IndexError, ValueError):
                    pass
            elif line.startswith('EFFN'):
                try:
                    metadata['effn'] = float(line.split()[1])
                except (IndexError, ValueError):
                    pass
            elif line.startswith('HMM'):
                # Reached the model section, stop parsing header
                break
    
    return metadata


def filter_profiles(input_dir, output_dir, min_nseq, report=False):
    """Filter profiles by NSEQ threshold."""
    
    input_path = Path(input_dir)
    output_path = Path(output_dir)
    
    if not input_path.exists():
        print(f"Error: Input directory '{input_dir}' does not exist.")
        sys.exit(1)
    
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Find all .hmm files
    hmm_files = list(input_path.glob('*.hmm'))
    
    if not hmm_files:
        print(f"Error: No .hmm files found in '{input_dir}'")
        sys.exit(1)
    
    print(f"Found {len(hmm_files)} HMM profiles")
    print(f"Filtering by NSEQ >= {min_nseq}")
    print("-" * 60)
    
    kept = []
    removed = []
    no_nseq = []
    
    for hmm_file in sorted(hmm_files):
        metadata = parse_hmm_header(hmm_file)
        
        if metadata['nseq'] is None:
            no_nseq.append((hmm_file.name, metadata))
            # Keep files without NSEQ (can't filter)
            shutil.copy(hmm_file, output_path / hmm_file.name)
            kept.append((hmm_file.name, metadata))
        elif metadata['nseq'] >= min_nseq:
            shutil.copy(hmm_file, output_path / hmm_file.name)
            kept.append((hmm_file.name, metadata))
        else:
            removed.append((hmm_file.name, metadata))
    
    # Summary
    print(f"\nResults:")
    print(f"  Kept:    {len(kept)} profiles")
    print(f"  Removed: {len(removed)} profiles")
    print(f"  No NSEQ: {len(no_nseq)} profiles (kept by default)")
    
    if report:
        print("\n" + "=" * 60)
        print("REMOVED PROFILES (NSEQ < {})".format(min_nseq))
        print("=" * 60)
        for filename, meta in sorted(removed, key=lambda x: x[1]['nseq'] or 0):
            print(f"  {filename}")
            print(f"    NAME: {meta['name']}")
            print(f"    NSEQ: {meta['nseq']}")
            print(f"    LENG: {meta['leng']}")
            print()
        
        print("\n" + "=" * 60)
        print("KEPT PROFILES (NSEQ >= {})".format(min_nseq))
        print("=" * 60)
        for filename, meta in sorted(kept, key=lambda x: -(x[1]['nseq'] or 0)):
            nseq_str = str(meta['nseq']) if meta['nseq'] is not None else "N/A"
            print(f"  {filename}: NSEQ={nseq_str}, LENG={meta['leng']}")
        
        if no_nseq:
            print("\n" + "=" * 60)
            print("PROFILES WITHOUT NSEQ (kept)")
            print("=" * 60)
            for filename, meta in no_nseq:
                print(f"  {filename}: LENG={meta['leng']}")
    
    # Write summary file
    summary_file = output_path / "filter_summary.txt"
    with open(summary_file, 'w') as f:
        f.write(f"HMM Profile Filter Summary\n")
        f.write(f"=" * 60 + "\n")
        f.write(f"Input directory:  {input_dir}\n")
        f.write(f"Output directory: {output_dir}\n")
        f.write(f"Minimum NSEQ:     {min_nseq}\n")
        f.write(f"\n")
        f.write(f"Total profiles:   {len(hmm_files)}\n")
        f.write(f"Kept:             {len(kept)}\n")
        f.write(f"Removed:          {len(removed)}\n")
        f.write(f"No NSEQ:          {len(no_nseq)}\n")
        f.write(f"\n")
        f.write(f"Removed profiles:\n")
        for filename, meta in sorted(removed, key=lambda x: x[1]['nseq'] or 0):
            f.write(f"  {filename} (NSEQ={meta['nseq']})\n")
    
    print(f"\nSummary written to: {summary_file}")
    print(f"Filtered profiles in: {output_dir}")
    
    return kept, removed


def main():
    parser = argparse.ArgumentParser(
        description='Filter HMM profiles by NSEQ (number of sequences used to build profile)'
    )
    parser.add_argument('input_dir', help='Input directory containing .hmm files')
    parser.add_argument('output_dir', help='Output directory for filtered profiles')
    parser.add_argument('--min-nseq', type=int, default=5,
                        help='Minimum NSEQ threshold (default: 5)')
    parser.add_argument('--report', action='store_true',
                        help='Print detailed report of kept/removed profiles')
    
    args = parser.parse_args()
    
    filter_profiles(args.input_dir, args.output_dir, args.min_nseq, args.report)


if __name__ == '__main__':
    main()
