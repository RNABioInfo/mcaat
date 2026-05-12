struct BeamPathInfo {
    std::vector<uint64_t> path;
    double total_mult;
};

class PhageCurator {
    // ...existing...
    std::vector<BeamPathInfo> BeamSearchPaths(uint64_t start, int length, int beam_width);
    // ...existing...
};