#ifndef SDBG_BUILD_H
#define SDBG_BUILD_H
#include <linux/limits.h>
#include <filesystem>
#include <sdbg/sdbg.h>
#include "settings.h"
#include <stdio.h>
//omp include
#include <omp.h>
#include <iostream>
#include <stdexcept>
#include <string>

#include "definitions.h"
#include "sorting/kmer_counter.h"
#include "sorting/read_to_sdbg.h"
#include "sorting/seq_to_sdbg.h"
#include "utils/options_description.h"
#include "utils/utils.h"

class SDBGBuild {
public:
    SDBGBuild(Settings settings);
    ~SDBGBuild();
private:
    std::ofstream lib_file;
    std::string WriteLibFile();
    void BuildLib();
    void BuildSDBG();
    bool BuildSDBGUtil(int argc, char **argv);
    bool BuildLibUtil(int argc, char **argv);
    Settings settings;
};

#endif