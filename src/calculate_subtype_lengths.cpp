#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <algorithm>
#include <filesystem>
#include <unistd.h>
#include "../include/profile.h"

namespace fs = std::filesystem;

// Helper to suppress stdout temporarily
class StdoutSuppressor {
    int saved_stdout;
    int saved_stderr;
public:
    StdoutSuppressor() {
        saved_stdout = dup(STDOUT_FILENO);
        saved_stderr = dup(STDERR_FILENO);
        freopen("/dev/null", "w", stdout);
        freopen("/dev/null", "w", stderr);
    }
    ~StdoutSuppressor() {
        dup2(saved_stdout, STDOUT_FILENO);
        dup2(saved_stderr, STDERR_FILENO);
        close(saved_stdout);
        close(saved_stderr);
    }
};

struct SubtypeGenes {
    std::string class_name;
    std::string type;
    std::string subtype;
    std::vector<std::string> genes;
};

// Map gene names to profile file patterns
std::map<std::string, std::vector<std::string>> gene_patterns = {
    {"Cas1", {"Cas1_"}},
    {"Cas2", {"Cas2_"}},
    {"Cas3", {"Cas3_", "Cas3HD_"}},
    {"Cas3HD", {"Cas3HD_"}},
    {"Cas4", {"Cas4_"}},
    {"Cas5a", {"Cas5_", "Cas5a_"}},
    {"Cas5b", {"Cas5_", "Cas5b_"}},
    {"Cas5c", {"Cas5_", "Cas5c_"}},
    {"Cas5d", {"Cas5d_"}},
    {"Cas5e", {"Cas5e_"}},
    {"Cas5f", {"Cas5f_"}},
    {"Cas5u", {"Cas5u_"}},
    {"Cas6", {"Cas6_"}},
    {"Cas6e", {"Cas6_", "Cas6e_"}},
    {"Cas6f", {"Cas6_", "Cas6f_"}},
    {"Cas7a", {"Cas7_", "Cas7a_"}},
    {"Cas7b", {"Cas7_", "Cas7b_"}},
    {"Cas7c", {"Cas7_", "Cas7c_"}},
    {"Cas7d", {"Cas7d_"}},
    {"Cas7e", {"Cas7e_"}},
    {"Cas7f", {"Cas7f_"}},
    {"Cas7u", {"Cas7u_"}},
    {"Cas8a", {"Cas8a_"}},
    {"Cas8b", {"Cas8b_"}},
    {"Cas8c", {"Cas8c_"}},
    {"Cas8e", {"Cas8e_"}},
    {"Cas8f", {"Cas8f_"}},
    {"Cas8u", {"Cas8u_"}},
    {"Cas9", {"Cas9_"}},
    {"Cas10", {"Cas10_"}},
    {"Cas10d", {"Cas10_", "Cas10d_"}},
    {"Cas11b", {"Cas11_", "Cas11b_"}},
    {"Cas11d", {"Cas11_", "Cas11d_"}},
    {"Cas11d2", {"Cas11_", "Cas11d2_"}},
    {"Cas11e", {"Cas11e_"}},
    {"Cas12a", {"Cas12a_", "Cpf1_"}},
    {"Cas12b1", {"Cas12b_"}},
    {"Cas12b2", {"Cas12b_"}},
    {"Cas12c", {"Cas12c_"}},
    {"Cas12d", {"Cas12d_"}},
    {"Cas12e", {"Cas12e_"}},
    {"Cas12f1", {"Cas12f_"}},
    {"Cas12f2", {"Cas12f_"}},
    {"Cas12f3", {"Cas12f_"}},
    {"Cas12g", {"Cas12g_"}},
    {"Cas12h", {"Cas12h_"}},
    {"Cas12i", {"Cas12i_"}},
    {"Cas12j", {"Cas12j_"}},
    {"Cas12k", {"Cas12k_"}},
    {"Cas12l", {"Cas12l_"}},
    {"Cas12m", {"Cas12m_"}},
    {"Cas13a", {"Cas13a_"}},
    {"Cas13b", {"Cas13b_"}},
    {"Cas13c", {"Cas13c_"}},
    {"Cas13d", {"Cas13d_"}},
    {"Csa5", {"Csa5_"}},
    {"Csm2", {"Csm2_"}},
    {"Csm3", {"Csm3_"}},
    {"Csm4", {"Csm4_"}},
    {"Csm5", {"Csm5_"}},
    {"Cmr1", {"Cmr1_"}},
    {"Cmr3", {"Cmr3_"}},
    {"Cmr4", {"Cmr4_"}},
    {"Cmr5", {"Cmr5_"}},
    {"Cmr6", {"Cmr6_"}},
    {"Cmr7", {"Cmr7_"}},
    {"Cmr8", {"Cmr8_"}},
    {"Csf1", {"Csf1_"}},
    {"Csf2", {"Csf2_"}},
    {"Csf3", {"Csf3_"}},
    {"Csf4", {"Csf4_"}},
    {"Csn2", {"Csn2_"}},
    {"Cas3-Cas2", {"Cas3-Cas2_"}},
    {"Csx10", {"Csx10_"}},
    {"Csx19", {"Csx19_"}},
    {"gRAMP", {"gRAMP_"}},
    {"SSgr11", {"SSgr11_"}},
    {"TniQ", {"TniQ_"}},
    {"RecD", {"RecD_"}},
};

std::vector<SubtypeGenes> parseSubtypes() {
    std::vector<SubtypeGenes> subtypes = {
        {"Class 1", "Type I", "I-A", {"Cas3", "Cas3HD", "Cas5a", "Cas6", "Cas7a", "Cas8a", "Csa5", "Cas1", "Cas2", "Cas4"}},
        {"Class 1", "Type I", "I-B", {"Cas3", "Cas5b", "Cas6", "Cas7b", "Cas8b", "Cas11b", "Cas1", "Cas2", "Cas4"}},
        {"Class 1", "Type I", "I-C", {"Cas3", "Cas5c", "Cas7c", "Cas8c", "Cas1", "Cas2", "Cas4"}},
        {"Class 1", "Type I", "I-D", {"Cas3", "Cas5d", "Cas6", "Cas7d", "Cas10d", "Cas11d", "Cas11d2", "Cas1", "Cas2", "Cas4"}},
        {"Class 1", "Type I", "I-E", {"Cas3", "Cas3HD", "Cas5e", "Cas6e", "Cas7e", "Cas8e", "Cas11e", "Cas1", "Cas2"}},
        {"Class 1", "Type I", "I-F", {"Cas5f", "Cas6f", "Cas7f", "Cas8f", "Cas3-Cas2", "Cas1"}},
        {"Class 1", "Type I", "I-F_T", {"Cas5f", "Cas6f", "Cas7f", "Cas8f", "TniQ"}},
        {"Class 1", "Type I", "I-G", {"Cas3", "Cas5u", "Cas7u", "Cas8u", "Cas1", "Cas2", "Cas4"}},
        
        {"Class 1", "Type III", "III-A", {"Cas10", "Csm2", "Csm3", "Csm4", "Csm5", "Cas1", "Cas2"}},
        {"Class 1", "Type III", "III-AM", {"Cas10", "Csm2", "Csm3", "Csm4", "Csm5", "Cas1", "Cas2"}},
        {"Class 1", "Type III", "III-B", {"Cas10", "Cmr1", "Cmr3", "Cmr4", "Cmr5", "Cmr6", "Cmr7", "Cmr8", "Cas1", "Cas2"}},
        {"Class 1", "Type III", "III-C", {"Cas10", "Cmr1", "Cmr3", "Cmr4", "Cmr5", "Cmr6"}},
        {"Class 1", "Type III", "III-D", {"Cas10", "Csm2", "Csm3", "Csm5", "Csx10", "Csx19", "Cas1", "Cas2"}},
        {"Class 1", "Type III", "III-E", {"gRAMP"}},
        {"Class 1", "Type III", "III-F", {"SSgr11", "Cas10", "Cas5a", "Csm3"}},
        
        {"Class 1", "Type IV", "IV-A1", {"Csf1", "Csf2", "Csf3", "Csf4", "Cas6"}},
        {"Class 1", "Type IV", "IV-A2", {"Csf2", "Csf3", "Csf4", "Cas6"}},
        {"Class 1", "Type IV", "IV-A3", {"Csf1", "Csf2", "Csf3", "Csf4", "Cas6"}},
        {"Class 1", "Type IV", "IV-B", {"Cas11b", "Csf1", "Csf2", "Csf3"}},
        {"Class 1", "Type IV", "IV-C", {"Csf2", "Csf3", "Cas10", "Cas11b", "Cas6"}},
        {"Class 1", "Type IV", "IV-D", {"Csf1", "Csf2", "Csf3", "RecD", "Cas6"}},
        {"Class 1", "Type IV", "IV-E", {"Csf2", "Csf3", "Csf4", "Cas6"}},
        
        {"Class 2", "Type II", "II-A", {"Cas9", "Cas1", "Cas2", "Csn2"}},
        {"Class 2", "Type II", "II-B", {"Cas9", "Cas1", "Cas2", "Cas4"}},
        {"Class 2", "Type II", "II-C", {"Cas9", "Cas1", "Cas2"}},
        {"Class 2", "Type II", "II-C2", {"Cas9", "Cas1", "Cas2"}},
        {"Class 2", "Type II", "II-D", {"Cas9", "Cas1", "Cas2"}},
        
        {"Class 2", "Type V", "V-A", {"Cas12a", "Cas1", "Cas2", "Cas4"}},
        {"Class 2", "Type V", "V-B1", {"Cas12b1", "Cas1", "Cas2", "Cas4"}},
        {"Class 2", "Type V", "V-B2", {"Cas12b2", "Cas1", "Cas2", "Cas4"}},
        {"Class 2", "Type V", "V-C", {"Cas12c", "Cas1"}},
        {"Class 2", "Type V", "V-D", {"Cas12d"}},
        {"Class 2", "Type V", "V-E", {"Cas12e", "Cas1", "Cas2", "Cas4"}},
        {"Class 2", "Type V", "V-F1", {"Cas12f1", "Cas1", "Cas2", "Cas4"}},
        {"Class 2", "Type V", "V-F2", {"Cas12f2", "Cas1", "Cas2", "Cas4"}},
        {"Class 2", "Type V", "V-F3", {"Cas12f3", "Cas1", "Cas2", "Cas4"}},
        {"Class 2", "Type V", "V-G", {"Cas12g"}},
        {"Class 2", "Type V", "V-H", {"Cas12h"}},
        {"Class 2", "Type V", "V-I", {"Cas12i"}},
        {"Class 2", "Type V", "V-J", {"Cas12j"}},
        {"Class 2", "Type V", "V-K", {"Cas12k", "TniQ"}},
        {"Class 2", "Type V", "V-L", {"Cas12l", "Cas1", "Cas2", "Cas4"}},
        {"Class 2", "Type V", "V-M", {"Cas12m"}},
        
        {"Class 2", "Type VI", "VI-A", {"Cas13a"}},
        {"Class 2", "Type VI", "VI-B1", {"Cas13b"}},
        {"Class 2", "Type VI", "VI-B2", {"Cas13b"}},
        {"Class 2", "Type VI", "VI-C", {"Cas13c"}},
        {"Class 2", "Type VI", "VI-D", {"Cas13d"}},
    };
    
    return subtypes;
}

std::string findProfileFile(const std::string& gene_name, const std::string& subtype, const std::string& profiles_dir) {
    std::string subtype_suffix = subtype;
    std::replace(subtype_suffix.begin(), subtype_suffix.end(), '-', '_');
    
    auto patterns_it = gene_patterns.find(gene_name);
    if (patterns_it == gene_patterns.end()) {
        return "";
    }
    
    for (const auto& pattern : patterns_it->second) {
        for (const auto& entry : fs::directory_iterator(profiles_dir)) {
            if (!entry.is_regular_file()) continue;
            
            std::string filename = entry.path().filename().string();
            if (filename.find(pattern) == 0 && filename.find(subtype) != std::string::npos && 
                filename.size() >= 4 && filename.substr(filename.size() - 4) == ".hmm") {
                return entry.path().string();
            }
        }
        
        // Try without subtype specificity
        for (const auto& entry : fs::directory_iterator(profiles_dir)) {
            if (!entry.is_regular_file()) continue;
            
            std::string filename = entry.path().filename().string();
            if (filename.find(pattern) == 0 && 
                filename.size() >= 4 && filename.substr(filename.size() - 4) == ".hmm") {
                return entry.path().string();
            }
        }
    }
    
    return "";
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <profiles_directory>" << std::endl;
        return 1;
    }
    
    std::string profiles_dir = argv[1];
    
    auto subtypes = parseSubtypes();
    
    std::cout << "CRISPR-Cas Subtype Total HMM Lengths\n";
    std::cout << "=====================================\n\n";
    
    std::string current_class = "";
    std::string current_type = "";
    
    for (const auto& subtype_info : subtypes) {
        if (subtype_info.class_name != current_class) {
            current_class = subtype_info.class_name;
            std::cout << "\n" << current_class << ":\n";
            std::cout << std::string(current_class.length() + 1, '-') << "\n";
        }
        
        if (subtype_info.type != current_type) {
            current_type = subtype_info.type;
            std::cout << "\n  " << current_type << ":\n";
        }
        
        int total_length = 0;
        std::vector<std::string> found_genes;
        std::vector<std::string> missing_genes;
        
        for (const auto& gene : subtype_info.genes) {
            std::string profile_file = findProfileFile(gene, subtype_info.subtype, profiles_dir);
            
            if (!profile_file.empty()) {
                Profile profile;
                {
                    StdoutSuppressor suppress;
                    if (profile.LoadFromFile(profile_file)) {
                        total_length += profile.GetLength();
                        found_genes.push_back(gene);
                    } else {
                        missing_genes.push_back(gene);
                    }
                }
            } else {
                missing_genes.push_back(gene);
            }
        }
        
        std::cout << "    " << subtype_info.subtype << ": ";
        
        for (size_t i = 0; i < found_genes.size(); i++) {
            std::cout << found_genes[i];
            if (i < found_genes.size() - 1) std::cout << ", ";
        }
        
        std::cout << " = " << total_length << " AA";
        
        if (!missing_genes.empty()) {
            std::cout << " (missing: ";
            for (size_t i = 0; i < missing_genes.size(); i++) {
                std::cout << missing_genes[i];
                if (i < missing_genes.size() - 1) std::cout << ", ";
            }
            std::cout << ")";
        }
        
        std::cout << "\n";
    }
    
    std::cout << "\n";
    
    return 0;
}
