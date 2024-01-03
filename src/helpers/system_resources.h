#include <iostream>
#include <fstream>
#include <algorithm>
#include <sys/sysinfo.h>
#include <pthread.h>
#include <thread>

using namespace std;

class SystemResources {
    
    private:
        long total_ram;
        int num_threads;
        int memory_percentage;
        int thread_count;
        long memory_to_use;
        int threads_to_use;
        
        long _GetTotalRam() {
            struct sysinfo sys_info;
            if(sysinfo(&sys_info) != 0) {
                cerr << "sysinfo: error reading system statistics";
                return -1;
            }
            return sys_info.totalram / 1024 / 1024;  // Convert bytes to MB
        }

        int _GetNumThreads() {
            return thread::hardware_concurrency();
        }

    public:
        SystemResources() {
            total_ram = this->_GetTotalRam();
            num_threads = this->_GetNumThreads();
        }

        

        void AskUser() {
            
            cout << "Enter the percentage of memory you want to use (1-100): ";
            cin >> memory_percentage;
            
            cout << "Enter the number of threads you want to use (1-" << num_threads << "): ";
            cin >> thread_count;
            this->CalculateUsage();
            cout << "You have chosen to use " << memory_to_use << " of memory and " 
                << thread_count << " threads." << endl;
            
        }

        void CalculateUsage() {
            memory_to_use = total_ram * memory_percentage;
            threads_to_use = min(thread_count, num_threads);
        }

        void WriteToFile() {
            ofstream file("data/system_resources.txt");
            file << "You have chosen to use " << memory_to_use << " MB of memory and " 
                << threads_to_use << " threads." << endl;
            file.close();
        }
    };
