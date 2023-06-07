///
/// INCLUSIONS
///

// Standard inclusions
#include <iostream>
#include <string>
#include <vector>


using namespace std;

// Local inclusions
#include "Properties.h"
#include "Lattice.h"
#include "Functions.h"

///
/// DEFINE PROPERTIES
///
LatProp LProp;
vector<PopProp> PProp(1);

///
/// MAIN CODE
///

int main(int argc, char **argv) {
    int task_id = atoi(argv[1]);
	// Choose parameters
    int Repeat = 10;
	long int MaxIterate = 1E11;
    int AdapNum = 7;
    vector<double> move_params{1E-2, 1};
    vector<double> sdeath_params{0.9, 1.16, 1.48, 1.9, 2.41, 3.06, 3.88, 4.91, 6.2, 7.84, 9.9};
    vector<string> move_names{"1E2", "1"};
    vector<string> sdeath_names{"9E1", "1d16", "1d48", "1d9", "2d41", "3d06", "3d88", "4", "6", "7", "9"};
    // Compute the indices
    int sdeath_num = sdeath_params.size();
    int move_index = task_id/(Repeat*sdeath_num);
    int sdeath_index = task_id % sdeath_num;

    // Set parameters
	PProp[0].StressBirthRate = 1;
	PProp[0].Move = move_params[move_index];
    PProp[0].StressDeathRate = sdeath_params[sdeath_index];
    string FileSpecifier = to_string(task_id)+"_M"+move_names[move_index]+"_SB1_SD"+sdeath_names[sdeath_index];
	// Run AdapImage
    AdapImage(LProp, PProp, AdapNum, FileSpecifier, MaxIterate);
	return 0;
}
