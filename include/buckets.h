/**
 * @file buckets.h
 * @brief HMM profile size buckets - dynamically loaded from HMM files
 */

#ifndef INCLUDE_BUCKETS_H_
#define INCLUDE_BUCKETS_H_

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <iostream>
#include <algorithm>

namespace HMMProfiles {

struct ProfileSize {
    std::string filename;
    int leng;       // LENG value from HMM file
    int min_aa;     // Minimum amino acids (-15%)
    int max_aa;     // Maximum amino acids (+25%)
    int min_bp;     // Minimum base pairs (min_aa * 3)
    int max_bp;     // Maximum base pairs (max_aa * 3)
};

// Parse LENG from an HMM file
inline int ParseLengFromHMM(const std::string& filepath) {
    std::ifstream file(filepath);
    if (!file.is_open()) return 0;
    
    std::string line;
    while (std::getline(file, line)) {
        if (line.substr(0, 4) == "LENG") {
            std::istringstream iss(line);
            std::string key;
            int leng;
            iss >> key >> leng;
            return leng;
        }
        // Stop at HMM line (header section ends)
        if (line.substr(0, 3) == "HMM") break;
    }
    return 0;
}

// Load all profiles from a directory
inline std::vector<ProfileSize> LoadProfiles(const std::string& profiles_dir) {
    std::vector<ProfileSize> profiles;
    
    // Check if directory exists
    if (!std::filesystem::exists(profiles_dir)) {
        std::cerr << "ERROR: Profiles directory does not exist: " << profiles_dir << std::endl;
        return profiles;
    }
    
    if (!std::filesystem::is_directory(profiles_dir)) {
        std::cerr << "ERROR: Not a directory: " << profiles_dir << std::endl;
        return profiles;
    }
    
    try {
        int file_count = 0;
        for (const auto& entry : std::filesystem::directory_iterator(profiles_dir)) {
            file_count++;
            std::string ext = entry.path().extension().string();
            if (ext == ".hmm") {
                std::string filename = entry.path().filename().string();
                std::string filepath = entry.path().string();
                
                int leng = ParseLengFromHMM(filepath);
                if (leng > 0) {
                    ProfileSize ps;
                    ps.filename = filename;
                    ps.leng = leng;
                    ps.min_aa = static_cast<int>(leng * 0.85);
                    ps.max_aa = static_cast<int>(leng * 1.25);
                    ps.min_bp = ps.min_aa * 3;
                    ps.max_bp = ps.max_aa * 3;
                    profiles.push_back(ps);
                } else {
                    std::cerr << "WARNING: Could not parse LENG from: " << filepath << std::endl;
                }
            }
        }
        if (profiles.empty() && file_count > 0) {
            std::cerr << "WARNING: Found " << file_count << " files but no .hmm files" << std::endl;
        }
    } catch (const std::filesystem::filesystem_error& e) {
        std::cerr << "ERROR: Failed to iterate directory: " << e.what() << std::endl;
        return profiles;
    }
    
    // Sort by LENG
    std::sort(profiles.begin(), profiles.end(), 
              [](const ProfileSize& a, const ProfileSize& b) { return a.leng < b.leng; });
    
    return profiles;
}

}  // namespace HMMProfiles

#endif  // INCLUDE_BUCKETS_H_
