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
    // vector<string> file_specifiers{"M1E2_M1E4_S1E3_D1E1", "M1E2_M1E4_S5_D3E1", "M1E2_M1E4_S5_D1E1",
    //"M3_M1E4_S1E3_D3E1", "M3_M1E4_S1E3_D1E1", "M3_M1E4_S5_D3E1", "M3_M1E4_S5_D1E1",
    //"M3_M1_S1E3_D3E1", "M3_M1_S1E3_D1E1", "M3_M1_S5_D3E1", "M3_M1_S5_D1E1",
    //"M1E2_C7E1", "M1_C7E1", "M1E2_C3E1", "M1_C3E1", "M1E2_SD7E1", "M1_SD7E1",
    //"SwM1E2_SwD1Ep3_D1E1", "SwM3_SwD1Ep3_D1E1", "SwM1E2_SwD1Ep4_D1E1",
    //"SwM3_SwD1Ep4_D1E1", "SwM1E2_SwD1Ep3_D3E1", "SwM3_SwD1Ep3_D3E1",
    //"SwM1E2_SwD1Ep4_D3E1", "SwM3_SwD1Ep4_D3E1"};
    // vector<string> file_specifiers{"SwM3_SwD5Ep4_D1E1", "SwM3_SwD5Ep4_D3E1", "SwM1E2_SwD5Ep4_D1E1", "SwM1E2_SwD5Ep4_D3E1"};
    // vector<string> file_specifiers{"M10_D4E1_R3", "M10_D4E1_R1", "M1E2_RC2E1", "M1_RC2E1"};
    // vector<string> file_specifiers{"M10_D4E1_R1", "M10_D4E1_R2", "M10_D4E1_R3"};
    // vector<string> file_specifiers{"M20_D3E1_R1", "M20_D3E1_R2", "M20_D3E1_R3"};
    // vector<string> file_specifiers{"SwM3_SwD5Ep4_D1E1_SwB1E3", "SwM3_SwD5Ep4_D3E1_SwB1E3", "SwM1E2_SwD5Ep4_D1E1_SwB1E3", "SwM1E2_SwD5Ep4_D3E1_SwB1E3"};
    //vector<string> file_specifiers{"M1_C6E1", "M1_C2E1", "M2_M1E2_S1E3", "M2_M1E2_S5"};
    vector<string> file_specifiers{"SwM3_SwD5Ep4_D1E1_SwB1E3", "SwM3_SwD5Ep4_D3E1_SwB1E3", "SwM1E2_SwD5Ep4_D1E1_SwB1E3", "SwM1E2_SwD5Ep4_D3E1_SwB1E3"};
    // For each file_specifier create three version of folders with bash files
    // Folder: Snap_filespecifier_V1, ...
    // Bash: BashI.sh
    string dir_name;
    for (int version=1;version<2;version++){
        for(int file=0;file < file_specifiers.size(); file++) {
            // Create directory
            dir_name = "Snap_"+file_specifiers[file]+"_V"+to_string(version);
            fs::create_directory(dir_name);
            // Create bash file
            string file_name = dir_name+"/BashI.sh";
	        ofstream BashFile(file_name);
            // Write the bash file
            BashFile << "#!/bin/bash\n";
            BashFile << "#SBATCH --time=0-8:00:00\n";
            BashFile << "#SBATCH -n 1\n\n";
            BashFile << "/store/DAMTP/vp381/SnapshotBuild/SnapCode/Snap_"+file_specifiers[file]+".out";
        }
    }
	return 0;
}