#!/usr/bin/env python3
"""
Generate buckets.h and updated subtype_lengths.txt from HMM profile files.
Reads LENG values from all .hmm files in profiles/ directory.
"""

import os
import re
import math
from pathlib import Path
from collections import defaultdict

def extract_leng_from_hmm(filepath):
    """Extract LENG value from HMM file."""
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('LENG'):
                match = re.search(r'LENG\s+(\d+)', line)
                if match:
                    return int(match.group(1))
    return None

def calculate_ranges(leng):
    """Calculate Min/Max AA and BP based on LENG value."""
    min_aa = int(leng * 0.85)  # -15%
    max_aa = int(leng * 1.25)  # +25%
    min_bp = min_aa * 3
    max_bp = max_aa * 3
    return min_aa, max_aa, min_bp, max_bp

def assign_to_bucket(leng):
    """Assign gene to size bucket based on LENG value (smaller buckets)."""
    if leng <= 75:
        return 1, "TINY", 30, 75, 90, 225
    elif leng <= 100:
        return 2, "VERY-SMALL", 75, 100, 225, 300
    elif leng <= 150:
        return 3, "SMALL", 100, 150, 300, 450
    elif leng <= 200:
        return 4, "SMALL-MEDIUM", 150, 200, 450, 600
    elif leng <= 250:
        return 5, "MEDIUM-SMALL", 200, 250, 600, 750
    elif leng <= 300:
        return 6, "MEDIUM", 250, 300, 750, 900
    elif leng <= 400:
        return 7, "MEDIUM-LARGE", 300, 400, 900, 1200
    elif leng <= 500:
        return 8, "LARGE-MEDIUM", 400, 500, 1200, 1500
    elif leng <= 650:
        return 9, "LARGE", 500, 650, 1500, 1950
    elif leng <= 800:
        return 10, "LARGE-XL", 650, 800, 1950, 2400
    elif leng <= 1000:
        return 11, "VERY-LARGE", 800, 1000, 2400, 3000
    elif leng <= 1200:
        return 12, "HUGE", 1000, 1200, 3000, 3600
    elif leng <= 1400:
        return 13, "VERY-HUGE", 1200, 1400, 3600, 4200
    else:
        return 14, "LARGEST", 1400, 2100, 4200, 6300

def main():
    profiles_dir = Path(__file__).parent / 'profiles'
    
    # Collect all HMM file data
    hmm_data = []
    
    for hmm_file in sorted(profiles_dir.glob('*.hmm')):
        leng = extract_leng_from_hmm(hmm_file)
        if leng is not None:
            min_aa, max_aa, min_bp, max_bp = calculate_ranges(leng)
            bucket_num, bucket_name, _, _, _, _ = assign_to_bucket(leng)
            
            hmm_data.append({
                'filename': hmm_file.name,
                'leng': leng,
                'min_aa': min_aa,
                'max_aa': max_aa,
                'min_bp': min_bp,
                'max_bp': max_bp,
                'bucket_num': bucket_num,
                'bucket_name': bucket_name
            })
    
    # Sort by LENG value
    hmm_data.sort(key=lambda x: x['leng'])
    
    # Group by buckets
    buckets = defaultdict(list)
    for item in hmm_data:
        buckets[item['bucket_num']].append(item)
    
    # Generate subtype_lengths.txt
    with open('docs/hmm_profile_lengths.txt', 'w') as f:
        f.write("HMM Profile Size Ranges (Auto-generated from profiles/*.hmm)\n")
        f.write("=" * 80 + "\n\n")
        f.write("Size range: -15% to +25% of LENG value\n")
        f.write("(accounts for natural variation and indels)\n\n")
        f.write("-" * 80 + "\n")
        f.write("ALL PROFILES TABLE\n")
        f.write("-" * 80 + "\n")
        f.write("Filename,LENG,Min_AA,Max_AA,Min_BP,Max_BP\n")
        
        for item in hmm_data:
            f.write(f"{item['filename']},{item['leng']},{item['min_aa']},"
                   f"{item['max_aa']},{item['min_bp']},{item['max_bp']}\n")
        
        f.write("\n\n")
        f.write("=" * 80 + "\n")
        f.write("SIZE BUCKETS FOR FAST FILTERING\n")
        f.write("=" * 80 + "\n\n")
        
        bucket_names = {
            1: (30, 75, 90, 225, "TINY"),
            2: (75, 100, 225, 300, "VERY-SMALL"),
            3: (100, 150, 300, 450, "SMALL"),
            4: (150, 200, 450, 600, "SMALL-MEDIUM"),
            5: (200, 250, 600, 750, "MEDIUM-SMALL"),
            6: (250, 300, 750, 900, "MEDIUM"),
            7: (300, 400, 900, 1200, "MEDIUM-LARGE"),
            8: (400, 500, 1200, 1500, "LARGE-MEDIUM"),
            9: (500, 650, 1500, 1950, "LARGE"),
            10: (650, 800, 1950, 2400, "LARGE-XL"),
            11: (800, 1000, 2400, 3000, "VERY-LARGE"),
            12: (1000, 1200, 3000, 3600, "HUGE"),
            13: (1200, 1400, 3600, 4200, "VERY-HUGE"),
            14: (1400, 2100, 4200, 6300, "LARGEST")
        }
        
        for bucket_num in sorted(buckets.keys()):
            min_aa, max_aa, min_bp, max_bp, name = bucket_names[bucket_num]
            f.write(f"Bucket {bucket_num}: {min_aa}-{max_aa} aa ({min_bp}-{max_bp} bp) — {name}\n")
            f.write(f"  Count: {len(buckets[bucket_num])} profiles\n")
            
            # Show first 10 files as examples
            examples = buckets[bucket_num][:10]
            f.write("  Examples: " + ", ".join(item['filename'] for item in examples))
            if len(buckets[bucket_num]) > 10:
                f.write(f", ... (+{len(buckets[bucket_num]) - 10} more)")
            f.write("\n\n")
        
        f.write("\n")
        f.write("=" * 80 + "\n")
        f.write("SUMMARY STATISTICS\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"Total profiles: {len(hmm_data)}\n")
        f.write(f"Smallest LENG: {hmm_data[0]['leng']} aa ({hmm_data[0]['filename']})\n")
        f.write(f"Largest LENG: {hmm_data[-1]['leng']} aa ({hmm_data[-1]['filename']})\n")
        f.write(f"\nOverall range to consider: {hmm_data[0]['min_aa']}-{hmm_data[-1]['max_aa']} aa ")
        f.write(f"({hmm_data[0]['min_bp']}-{hmm_data[-1]['max_bp']} bp)\n")
    
    # Generate buckets.h
    with open('include/buckets.h', 'w') as f:
        f.write("""/**
 * @file buckets.h
 * @brief HMM profile size buckets (auto-generated from profiles/*.hmm)
 * 
 * This file is AUTO-GENERATED by generate_buckets.py
 * DO NOT EDIT MANUALLY - your changes will be overwritten
 * 
 * Contains size ranges for all HMM profiles based on LENG values.
 * Size range: -15% to +25% of LENG value
 */

#ifndef INCLUDE_BUCKETS_H_
#define INCLUDE_BUCKETS_H_

#include <string>
#include <vector>
#include <cstdint>

namespace HMMProfiles {

/**
 * @brief HMM profile size information
 */
struct ProfileSize {
    std::string filename;
    int leng;       // LENG value from HMM file
    int min_aa;     // Minimum amino acids (-15%)
    int max_aa;     // Maximum amino acids (+25%)
    int min_bp;     // Minimum base pairs (min_aa * 3)
    int max_bp;     // Maximum base pairs (max_aa * 3)
};

/**
 * @brief Size bucket for fast filtering
 */
struct SizeBucket {
    int bucket_num;
    std::string name;
    int min_aa;
    int max_aa;
    int min_bp;
    int max_bp;
    std::vector<std::string> profiles;
};

""")
        
        # Write ALL_PROFILES array
        f.write("/**\n * @brief All HMM profiles with size ranges\n */\n")
        f.write("static const std::vector<ProfileSize> ALL_PROFILES = {\n")
        for i, item in enumerate(hmm_data):
            comma = "," if i < len(hmm_data) - 1 else ""
            f.write(f'    {{"{item["filename"]}", {item["leng"]}, {item["min_aa"]}, '
                   f'{item["max_aa"]}, {item["min_bp"]}, {item["max_bp"]}}}{comma}\n')
        f.write("};\n\n")
        
        # Write SIZE_BUCKETS array
        f.write("/**\n * @brief Size buckets for fast filtering\n */\n")
        f.write("static const std::vector<SizeBucket> SIZE_BUCKETS = {\n")
        
        for i, bucket_num in enumerate(sorted(buckets.keys())):
            min_aa, max_aa, min_bp, max_bp, name = bucket_names[bucket_num]
            profiles = [f'"{item["filename"]}"' for item in buckets[bucket_num]]
            
            comma = "," if i < len(buckets) - 1 else ""
            f.write(f'    {{{bucket_num}, "{name}", {min_aa}, {max_aa}, {min_bp}, {max_bp},\n')
            f.write(f'        {{{", ".join(profiles)}}}}}{comma}\n')
        
        f.write("};\n\n")
        
        # Write constants
        f.write("/**\n * @brief Overall size range constants\n */\n")
        f.write(f"constexpr int MIN_PROFILE_AA = {hmm_data[0]['min_aa']};\n")
        f.write(f"constexpr int MAX_PROFILE_AA = {hmm_data[-1]['max_aa']};\n")
        f.write(f"constexpr int MIN_PROFILE_BP = {hmm_data[0]['min_bp']};\n")
        f.write(f"constexpr int MAX_PROFILE_BP = {hmm_data[-1]['max_bp']};\n")
        f.write(f"constexpr int TOTAL_PROFILES = {len(hmm_data)};\n\n")
        
        f.write("}  // namespace HMMProfiles\n\n")
        f.write("#endif  // INCLUDE_BUCKETS_H_\n")
    
    print(f"✓ Processed {len(hmm_data)} HMM profiles")
    print(f"✓ Generated docs/hmm_profile_lengths.txt")
    print(f"✓ Generated include/buckets.h")
    print(f"\nSize range: {hmm_data[0]['min_aa']}-{hmm_data[-1]['max_aa']} aa "
          f"({hmm_data[0]['min_bp']}-{hmm_data[-1]['max_bp']} bp)")

if __name__ == '__main__':
    main()
