#ifndef SETTINGS_H
#define SETTINGS_H

#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <thread>
#include <filesystem>
#include <chrono>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <stdexcept>
#include <cctype>
#include "sdbg/sdbg.h"
#include <fstream>

using namespace std;
namespace fs = std::filesystem;

struct Settings {
    std::string input_files; // Path to the input files
    double ram = 0.0; // Maximum RAM usage in gigabytes
    size_t threads = 0; // Number of threads to use
    std::string output_folder; // Output directory path. Empty until chosen by CLI or settings file.
    std::string graph_folder; // Folder for graph data
    std::string cycles_folder; // Folder for cycle data
    std::string output_file; // Path to the main output file
    std::string benchmark_file;
    // Sdbg

    struct CycleFinderSettings {
        uint64_t threshold_multiplicity = 20; // Minimum multiplicity threshold
        bool low_abundance = true; // Flag to enable low abundance mode
        int cycle_max_length = 77; // Maximum length of a cycle
        int cycle_min_length = 27; // Minimum length of a cycle
    } cycle_finder_settings;
    struct DNASequenceSettings {
        int spacer_min_length = 23; // Minimum length of a spacer
        int spacer_max_length = 50; // Maximum length of a spacer
        int repeat_min_length = 23; // Minimum length of a repeat
        int repeat_max_length = 50; // Maximum length of a repeat
    } dna_sequence_settings;

    SDBG* sdbg = nullptr; // Pointer to the SDBG graph

    Settings() {
        // Defaults are intentionally empty so parse_arguments can apply
        // a timestamp-based default only if neither CLI nor settings file specify the output folder.
        output_folder = "";
        graph_folder = "";
        cycles_folder = "";
        output_file = "";
    }

    // Generate timestamp in YYYY-MM-DD_HH-MM-SS format
    string get_timestamp() {
        auto now = chrono::system_clock::now();
        auto in_time_t = chrono::system_clock::to_time_t(now);
        stringstream ss;
        ss << put_time(localtime(&in_time_t), "%Y-%m-%d_%H-%M-%S");
        return ss.str();
    }

    // Get maximum available threads
    size_t max_threads() const {
        return thread::hardware_concurrency();
    }

    // Validate settings and return status map
    map<string, pair<bool, string>> validate_settings() const {
        map<string, pair<bool, string>> settings_message_map;
        
        // Check input files
        bool input_valid = !input_files.empty();
        settings_message_map["Input Files"] = make_pair(input_valid, 
            input_valid ? input_files + " exist(s)" : "No input files specified");

        // Check RAM (using reasonable default max of 128GB if not specified)
        string ram_str = to_string(ram).substr(0, to_string(ram).find(".") + 3);
        bool ram_valid = ram > 1.0;
        settings_message_map["RAM"] = make_pair(ram_valid, 
            ram_valid ? ram_str + " GB" : 
            "Value " + ram_str + " GB is invalid (must be greater than 1 GB)");

        // Check threads
        size_t max_t = max_threads();
        bool threads_valid = threads > 0 && threads <= max_t;
        settings_message_map["Threads"] = make_pair(threads_valid, 
            threads_valid ? to_string(threads) + " thread(s)" : 
            "Value " + to_string(threads) + " is invalid (must be between 1 and " + to_string(max_t) + ")");

        // Check folder paths
        bool output_valid = !output_folder.empty();
        settings_message_map["Output Folder"] = make_pair(output_valid, 
            output_valid ? output_folder : "Invalid output folder");

        return settings_message_map;
    }

    // Print settings summary and return erroneous properties
    string print_settings() {
        string erroneous_properties;
        auto settings_map = validate_settings();

        for (const auto& [key, value] : settings_map) {
            if (value.first) {
                cout << "[✔] " << key << ": " << value.second << endl;
            } else {
                erroneous_properties += key + " ";
                cout << "[✗] " << key << ": " << value.second << endl;
            }
        }
        return erroneous_properties;
    }

    // Read settings values from a simple key=value file. Lines starting with # or // are ignored.
    // Keys mirror struct property names; examples:
    // input_files=/path/a.fa /path/b.fa
    // ram=4G
    // threads=4
    // cycle_max_length=77
    // cycle_min_length=27
    // threshold_multiplicity=20
    // low_abundance=true
    bool LoadFromFile(const std::string& path) {
        std::ifstream file(path);
        if (!file.is_open()) {
            std::cerr << "Could not open settings file: " << path << std::endl;
            return false;
        }
        std::string line;
        auto trim = [](std::string s) {
            // trim in place
            const char* ws = " \t\n\r\f\v";
            s.erase(0, s.find_first_not_of(ws));
            s.erase(s.find_last_not_of(ws) + 1);
            return s;
        };
        while (std::getline(file, line)) {
            // Remove comments
            size_t posc = line.find('#');
            if (posc != std::string::npos) line = line.substr(0, posc);
            size_t pos2 = line.find("//");
            if (pos2 != std::string::npos) line = line.substr(0, pos2);
            string s = trim(line);
            if (s.empty()) continue;
            size_t eq = s.find('=');
            if (eq == std::string::npos) continue;
            string key = trim(s.substr(0, eq));
            string val = trim(s.substr(eq + 1));

            // Interpret known keys
            if (key == "input_files") {
                // allow space/comma/semicolon-separated file list; normalize to single-space separated
                vector<string> tokens;
                string cur;
                for (char c : val) {
                    if (c == ',' || c == ';') c = ' ';
                    if (!isspace(static_cast<unsigned char>(c))) {
                        cur.push_back(c);
                    } else {
                        if (!cur.empty()) {
                            tokens.push_back(cur);
                            cur.clear();
                        }
                    }
                }
                if (!cur.empty()) tokens.push_back(cur);
                // join with single space so older code (SDBGBuild) gets consistent formatting
                this->input_files.clear();
                for (size_t i = 0; i < tokens.size(); ++i) {
                    this->input_files += tokens[i];
                    if (i + 1 < tokens.size()) this->input_files += " ";
                }
            } else if (key == "ram") {
                // reuse same parsing as CLI: accept B/K/M/G suffix
                try {
                    double value = 0.0;
                    char unit = 'G';
                    size_t p = val.find_first_not_of("0123456789.");
                    if (p != std::string::npos) {
                        value = stod(val.substr(0, p));
                        unit = toupper(val[p]);
                    } else {
                        value = stod(val);
                    }
                    switch (unit) {
                        case 'B': this->ram = value / (1024.0 * 1024.0 * 1024.0); break;
                        case 'K': this->ram = value / (1024.0 * 1024.0); break;
                        case 'M': this->ram = value / 1024.0; break;
                        case 'G': this->ram = value; break;
                        default: throw runtime_error("Invalid RAM unit in settings file: " + val);
                    }
                } catch (...) {
                    std::cerr << "Warning: could not parse RAM value '" << val << "' in settings file" << std::endl;
                }
            } else if (key == "threads") {
                try { this->threads = stoul(val); } catch (...) { }
            } else if (key == "output_folder") { this->output_folder = val; }
            else if (key == "graph_folder") { this->graph_folder = val; }
            else if (key == "cycles_folder") { this->cycles_folder = val; }
            else if (key == "output_file") { this->output_file = val; }
            else if (key == "cycle_max_length") { this->cycle_finder_settings.cycle_max_length = stoi(val); }
            else if (key == "cycle_min_length") { this->cycle_finder_settings.cycle_min_length = stoi(val); }
            else if (key == "threshold_multiplicity") { this->cycle_finder_settings.threshold_multiplicity = stoull(val); }
            else if (key == "low_abundance") {
                std::transform(val.begin(), val.end(), val.begin(), ::tolower);
                this->cycle_finder_settings.low_abundance = (val == "true" || val == "1" || val == "yes");
            }
            else if (key == "spacer_min_length") { this->dna_sequence_settings.spacer_min_length = stoi(val); }
            else if (key == "spacer_max_length") { this->dna_sequence_settings.spacer_max_length = stoi(val); }
            else if (key == "repeat_min_length") { this->dna_sequence_settings.repeat_min_length = stoi(val); }
            else if (key == "repeat_max_length") { this->dna_sequence_settings.repeat_max_length = stoi(val); }
            // unknown keys are ignored for forward-compatibility
        }
        file.close();
        return true;
    }
};

#endif // SETTINGS_H