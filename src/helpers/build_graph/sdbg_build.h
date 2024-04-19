#include "../libs/megahit/src/main_sdbg_build.cpp"
#include "../libs/megahit/src/main_buildlib.cpp"
#include <iostream>
#include <fstream>
#include <string>
#include <../src/sdbg/sdbg.h>
#include "settings.h";

class SDBGBuild {
public:
    SDBGBuild();
    ~SDBGBuild();
private:
    std::ofstream lib_file;
    std::string _write_lib_file();
    void _build_lib(std::string lib_file_name);
    void _build_sdbg();
    settings settings;
};