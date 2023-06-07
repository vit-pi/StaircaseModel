///
/// INCLUSIONS
///

// Standard inclusions
#include <iostream>
#include <string>
#include <filesystem>
#include <vector>
#include <fstream>

using namespace std;
namespace fs = std::filesystem;

///
/// MAIN CODE
///

int main(void) {
	// File specifiers
    //vector<string> file_specifiers{"Mut", "DensMotSwitch", "CompBelowStair"};
    //vector<int> jobs{2339, 6759, 2859};
    vector<string> file_specifiers{"Standard", "Chemotax"};
    vector<int> jobs{259, 1169};
    // For each file_specifier create three version of folders with bash files
    // Folder: Snap_filespecifier_V1, ...
    // Bash: BashI.sh
    string dir_name;
    for (int version=1;version<6;version++){
        for(int file=0;file < file_specifiers.size(); file++) {
            // Create directory
            dir_name = "Adap"+file_specifiers[file]+"V"+to_string(version);
            fs::create_directory(dir_name);
            // Create bash files
            int bash_num = 1;
            while (bash_num <= 1+jobs[file]/5000){
                string file_name = dir_name+"/Bash"+to_string(bash_num)+".sh";
                ofstream BashFile(file_name);
                // Write the bash file
                BashFile << "#!/bin/bash\n";
                BashFile << "#SBATCH --time=0-18:00:00\n";
                BashFile << "#SBATCH -n 1\n";
                BashFile << "#SBATCH --partition=medium\n";
                if (bash_num == 1+jobs[file]/5000){
                    BashFile << "#SBATCH --array=0-"+to_string(jobs[file]%5000)+"\n";
                }
                else{
                    BashFile << "#SBATCH --array=0-4999\n";
                }
                BashFile << "/data/math-staircase-model/magd5825/StaircaseBuild/Code/"+file_specifiers[file]+".out $(($SLURM_ARRAY_TASK_ID+"+to_string((bash_num-1)*5000)+"))";
                bash_num++;
            }
        }
    }
	return 0;
}