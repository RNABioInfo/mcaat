#include "sdbg_build.h"

SDBGBuild::SDBGBuild(Settings settings) : settings(settings)
{
    //create a lib file
    //write to the lib file
    //build the lib file
    //build the sdbg

    this->BuildLib();
    this->BuildSDBG();
    
}

SDBGBuild::~SDBGBuild()
{
    //close the lib file
}

// @brief Write the lib file
// @details This function writes the lib file to the output folder
// @param None
// @return The path to the lib file

std::string SDBGBuild::WriteLibFile()
{   
    
    std::cout << "\n-----------------------------------------\n" << std::endl;
    std::cout << "2. Building the SDBG: " << std::endl;
    std::string data_lib_file_name = settings.graph_folder;  
    
    // Debugging: Log the intended directory
    std::cout << "Directory: " << data_lib_file_name << std::endl;

    namespace fs = std::filesystem;

    fs::path resolved_path_fs(data_lib_file_name);

    // Create all necessary parent directories
    try {
        if (!fs::exists(resolved_path_fs)) {
            fs::create_directories(resolved_path_fs);
            std::cout << "Directories created: " << resolved_path_fs << std::endl;
        } else {
            std::cout << "Directory already exists: " << resolved_path_fs << std::endl;
        }
    } catch (const fs::filesystem_error& e) {
        std::cerr << "Error creating directories: " << e.what() << std::endl;
        return "";
    }

    // Use std::filesystem to append to the path
    fs::path lib_file_path = resolved_path_fs / "data.lib";
    std::ofstream lib_file(lib_file_path, std::ios::out | std::ios::trunc);
    if (!lib_file) {
        std::cerr << "Error creating file: " << lib_file_path << std::endl;
        return "";
    }
    string switch_str = "se";
    if (settings.input_files.find(" ") != std::string::npos) {
        switch_str = "pe";
    }
    lib_file << "#lib file for the SDBG from " + settings.input_files + "\n";
    lib_file << switch_str+" " + settings.input_files;
    lib_file.close();
    std::cout << "Saved lib to file " << lib_file_path << std::endl;

    return lib_file_path.string();
}

// @brief Build the lib file from the input files
// @details This function builds the lib file from the input files and the settings
// provided by the user. It uses the buildlib tool from the Megahit library.
// @param None

void SDBGBuild::BuildLib()
{
    std::string path_name = WriteLibFile();
    // following are the arguments to call <tool_name> buildlib data.lib outfile_prefix
    char *argv[3];
    argv[0] = "buildlib";
    std::string full_out = settings.graph_folder + "/outfile_prefix";
    std::string outfile_prefix = settings.graph_folder + "/graph";
    //convert string to char*
    argv[1] = const_cast<char*>(path_name.c_str());
    argv[2] = const_cast<char*>(full_out.c_str());
    std::cout << "Building the library: " ;
    //print all the arguments
    for(int i=0; i<3; i++){
        std::cout << argv[i] << " ";
    }
    std::cout << std::endl;
    BuildLibUtil(3, argv);

    return;

}


bool SDBGBuild::BuildLibUtil(int argc, char **argv) {
  AutoMaxRssRecorder recorder;

  if (argc < 3) {
    return false;
  }
  SequenceLibCollection::Build(argv[1], argv[2]);

  return true;
}

// @brief Build the SDBGUtil
// @details This is a helper function to build the SDBG from the library file and the settings.
// This is a copy of the same function - main_build_sdbg defined in libs/megahit/src/main_sdbg_build.cpp
// @param argc The number of arguments 
// @param argv The arguments passed to the function
bool SDBGBuild::BuildSDBGUtil(int argc, char **argv) {
  AutoMaxRssRecorder recorder;

  // parse option the same as kmer_count
  OptionsDescription desc;
  Read2SdbgOption opt;

  desc.AddOption("kmer_k", "k", opt.k, "kmer size");
  desc.AddOption("min_kmer_frequency", "m", opt.solid_threshold,
                 "min frequency to output an edge");
  desc.AddOption(
      "host_mem", "", opt.host_mem,
      "Max memory to be used. 90% of the free memory is recommended.");
  desc.AddOption("num_cpu_threads", "", opt.n_threads,
                 "number of CPU threads. At least 2.");
  desc.AddOption("read_lib_file", "", opt.read_lib_file,
                 "input fast[aq] file, can be gzip'ed. \"-\" for stdin.");
  desc.AddOption("output_prefix", "", opt.output_prefix, "output prefix");
  desc.AddOption("mem_flag", "", opt.mem_flag,
                 "memory options. 0: minimize memory usage; 1: automatically "
                 "use moderate memory; "
                 "other: use all "
                 "available mem specified by '--host_mem'");
  desc.AddOption("need_mercy", "", opt.need_mercy, "to add mercy edges.");

  try {
    desc.Parse(argc, argv);

    if (opt.read_lib_file.empty()) {
      throw std::logic_error("No input file!");
    }

    if (opt.n_threads == 0) {
      opt.n_threads = omp_get_max_threads();
    }

    if (opt.host_mem == 0) {
      throw std::logic_error("Please specify the host memory!");
    }
  } catch (std::exception &e) {
    std::cerr << e.what() << std::endl;
    std::cerr
        << "Usage: sdbg_builder read2sdbg --read_lib_file fastx_file -o out"
        << std::endl;
    std::cerr << "Options:" << std::endl;
    std::cerr << desc << std::endl;
    exit(1);
  }

  SeqPkgWithSolidMarker pkg;

  {
    // stage 1
    Read2SdbgS1 runner(opt, &pkg);
    if (opt.solid_threshold > 1) {
      runner.Run();
    } else {
      runner.Initialize();
    }
  }

  {
    // stage 2
    Read2SdbgS2 runner(opt, &pkg);
    runner.Run();
  }

  return 0;
}

// @brief Build the SDBG
// @details This function builds the SDBG from the library file and the settings
// provided by the user. It uses the read2sdbg tool from the Megahit library.
// @param None
void SDBGBuild::BuildSDBG()
{
    // these are the arguments passed - make them into a char array
    // megahit_core read2sdbg --host_mem 128000000000 --read_lib_file outfile_prefix -m 1 -k 21 --num_cpu_threads 48
    // convert double to string 

    std::string host_mem = std::to_string(settings.ram);
    std::string num_cpu_threads = std::to_string(settings.threads);
    std::string read2sdbg = "read2sdbg";
    std::string host_mem_flag = "--host_mem";
    std::string read_lib_file_flag = "--read_lib_file";
    std::string full_out = settings.graph_folder + "/outfile_prefix";
    // convert string to char*  
    // Create a vector to hold the arguments
    std::vector<std::string> args = {
        read2sdbg,
        host_mem_flag,
        host_mem,
        read_lib_file_flag,
        full_out,
        "-m", "1",
        "-k", "23",
        "--num_cpu_threads", num_cpu_threads,
        "--o", settings.graph_folder + "/graph"
    };

    // Create a vector of char* to pass to the function
    std::vector<char*> argv;
    for (auto& arg : args) {
        argv.push_back(&arg[0]);
    }

    // Now argv can be used as a char* array
    this->BuildSDBGUtil(argv.size(), argv.data());
    std::cout << "\n-----------------------------------------\n" << std::endl;
    return;
}