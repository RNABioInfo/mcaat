/**
 * @file buckets.h
 * @brief HMM profile size buckets
 *
 * Contains size ranges for all HMM profiles based on LENG values.
 * Size range: -15% to +25% of LENG value
 */

#ifndef INCLUDE_BUCKETS_H_
#define INCLUDE_BUCKETS_H_

#include <string>
#include <vector>
#include <cstdint>
#include <fstream>
#include <filesystem>
#include <algorithm>
#include "profile.h"

namespace HMMProfiles {

struct ProfileSize {
    std::string filename;
    int leng;       // LENG value from HMM file
    int min_aa;     // Minimum amino acids (-15%)
    int max_aa;     // Maximum amino acids (+25%)
    int min_bp;     // Minimum base pairs (min_aa * 3)
    int max_bp;     // Maximum base pairs (max_aa * 3)
};

/**
 * @brief Load all HMM profiles from a directory
 *
 * Scans directory for .hmm files, reads LENG value from each,
 * calculates min/max bounds, and returns sorted list.
 *
 * @param profiles_dir Path to directory containing .hmm files
 * @param filter_pattern Optional: only load files matching this substring (e.g., "_0_" for primary variants)
 * @return Vector of ProfileSize structs, sorted by LENG (ascending)
 */
inline std::vector<ProfileSize> LoadProfilesFromDirectory(const std::string& profiles_dir,
                                                          const std::string& filter_pattern = "") {
    namespace fs = std::filesystem;
    std::vector<ProfileSize> profiles;

    if (!fs::exists(profiles_dir) || !fs::is_directory(profiles_dir)) {
        return profiles;  // Empty vector
    }

    // Scan directory for .hmm files
    for (const auto& entry : fs::directory_iterator(profiles_dir)) {
        if (!entry.is_regular_file()) continue;

        std::string filename = entry.path().filename().string();
        if (filename.size() < 4 || filename.substr(filename.size() - 4) != ".hmm") {
            continue;  // Not a .hmm file
        }

        // Apply filter if specified
        if (!filter_pattern.empty() && filename.find(filter_pattern) == std::string::npos) {
            continue;  // Doesn't match filter
        }

        // Fast read: only extract LENG line, don't parse entire HMM
        std::ifstream file(entry.path().string());
        if (!file.is_open()) continue;

        int leng = -1;
        std::string line;
        while (std::getline(file, line)) {
            // LENG line format: "LENG  123"
            if (line.substr(0, 4) == "LENG") {
                size_t pos = line.find_first_of("0123456789");
                if (pos != std::string::npos) {
                    leng = std::stoi(line.substr(pos));
                    break;  // Found it, stop reading
                }
            }
        }
        file.close();

        if (leng <= 0) continue;  // Invalid or not found

        // Calculate bounds: -15% to +25%
        ProfileSize ps;
        ps.filename = filename;
        ps.leng = leng;
        ps.min_aa = static_cast<int>(leng * 0.85);  // -15%
        ps.max_aa = static_cast<int>(leng * 1.25);  // +25%
        ps.min_bp = ps.min_aa * 3;
        ps.max_bp = ps.max_aa * 3;

        profiles.push_back(ps);
    }

    // Sort by LENG (ascending)
    std::sort(profiles.begin(), profiles.end(),
              [](const ProfileSize& a, const ProfileSize& b) {
                  return a.leng < b.leng;
              });

    return profiles;
}

}  // namespace HMMProfiles

#endif  // INCLUDE_BUCKETS_H_
