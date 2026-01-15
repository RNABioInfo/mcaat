#!/usr/bin/env python3
"""
Update buckets.h to remove filtered profiles.
Reads filter_summary.txt and removes corresponding profiles from buckets.h
"""

import re
from pathlib import Path

def extract_removed_profiles(filter_summary_path):
    """
    Extract list of removed profile filenames from filter_summary.txt
    
    Returns:
        set: Set of removed profile filenames
    """
    removed_profiles = set()
    
    with open(filter_summary_path, 'r') as f:
        in_removed_section = False
        for line in f:
            line = line.strip()
            
            # Check if we've entered the removed profiles section
            if line == "Removed profiles:":
                in_removed_section = True
                continue
            
            # Extract profile name from lines like "  Cas10_2_CAS-III.hmm (NSEQ=1)"
            if in_removed_section and line:
                match = re.match(r'\s*(\S+\.hmm)\s+\(NSEQ=\d+\)', line)
                if match:
                    removed_profiles.add(match.group(1))
    
    return removed_profiles

def update_buckets_h(buckets_path, removed_profiles):
    """
    Update buckets.h by removing profiles in the removed_profiles set
    
    Args:
        buckets_path: Path to buckets.h file
        removed_profiles: Set of profile filenames to remove
    """
    with open(buckets_path, 'r') as f:
        lines = f.readlines()
    
    new_lines = []
    removed_count = 0
    kept_count = 0
    in_all_profiles = False
    
    for i, line in enumerate(lines):
        # Detect if we're in the ALL_PROFILES vector
        if 'static const std::vector<ProfileSize> ALL_PROFILES = {' in line:
            in_all_profiles = True
            new_lines.append(line)
            continue
        
        # Detect end of ALL_PROFILES
        if in_all_profiles and line.strip() == '};':
            in_all_profiles = False
            new_lines.append(line)
            continue
        
        # Process profile entries
        if in_all_profiles:
            # Check if this line contains a profile entry
            match = re.search(r'{"([^"]+\.hmm)"', line)
            if match:
                profile_name = match.group(1)
                if profile_name in removed_profiles:
                    # Skip this line (don't add to new_lines)
                    removed_count += 1
                    print(f"Removing: {profile_name}")
                    continue
                else:
                    kept_count += 1
                    new_lines.append(line)
            else:
                # Not a profile entry, keep the line
                new_lines.append(line)
        else:
            # Not in ALL_PROFILES section
            # Update TOTAL_PROFILES constant
            if 'constexpr int TOTAL_PROFILES' in line:
                line = f'constexpr int TOTAL_PROFILES = {kept_count};\n'
            new_lines.append(line)
    
    # Write updated file
    with open(buckets_path, 'w') as f:
        f.writelines(new_lines)
    
    print(f"\nSummary:")
    print(f"  Profiles kept:    {kept_count}")
    print(f"  Profiles removed: {removed_count}")
    print(f"  Updated TOTAL_PROFILES constant to: {kept_count}")
    print(f"\nUpdated {buckets_path}")

def main():
    script_dir = Path(__file__).parent
    filter_summary_path = script_dir / 'filtered_profiles' / 'filter_summary.txt'
    buckets_path = script_dir / 'include' / 'buckets.h'
    
    print(f"Reading filter summary from: {filter_summary_path}")
    removed_profiles = extract_removed_profiles(filter_summary_path)
    print(f"Found {len(removed_profiles)} profiles to remove")
    
    print(f"\nUpdating buckets.h at: {buckets_path}")
    update_buckets_h(buckets_path, removed_profiles)

if __name__ == '__main__':
    main()
