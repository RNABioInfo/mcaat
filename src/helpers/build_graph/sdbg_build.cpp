#include "sdbg_builder.h"

SDBGBuild::SDBGBuild()
{
}

SDBGBuild::~SDBGBuild()
{
    this->_build_lib();
    this->_build_sdbg();
}

std::string SDBGBuild::_write_lib_file()
{
    //create a file under with the name settings.output_graph
    std::ofstream lib_file(settings.output_graph+".lib");
    lib_file << "#lib file for the SDBG from \n";
    //choose after last / in settings.input_folder
    int last_slash = settings.output_graph.find_last_of("/");
    //substring from last_slash to end
    lib_file << settings.input_folder.substr(last_slash+1);
    lib_file.close();
    return settings.output_graph+".lib";
}

void SDBGBuild::_build_lib(std::string lib_file_name)
{
    std::string path_name = _write_lib_file();
    char *argv[3];
    argv[0] = "buildlib";
    int last_slash = settings.output_graph.find_last_of("/");
    //convert string to char*
    argv[1] = const_cast<char*>(path_name.c_str());
    argv[2] = const_cast<char*>(settings.output_graph.c_str());
    
    main_build_lib(3, argv);

    return;

}

void SDBGBuild::_build_sdbg()
{
    // these are the arguments passed - make them into a char array
    // megahit_core read2sdbg --host_mem 128000000000 --read_lib_file outfile_prefix -m 1 -k 21 --num_cpu_threads 48

    char *argv[11];
    argv[0] = "read2sdbg";
    argv[1] = "--host_mem";
    argv[2] = const_cast<char*>(settings.host_memory.c_str());
    argv[3] = "--read_lib_file";
    argv[4] = const_cast<char*>(settings.output_graph.c_str());
    argv[5] = "-m";
    argv[6] = const_cast<char*>(settings.mem.c_str());
    argv[7] = "-k";
    argv[8] = const_cast<char*>(settings.kmer_size.c_str());
    argv[9] = "--num_cpu_threads";
    argv[10] = const_cast<char*>(settings.threads.c_str());
    main_read2sdbg(11, argv);
    return;
 }