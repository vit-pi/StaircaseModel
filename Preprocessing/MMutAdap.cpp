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
    vector<double> death_params{1E-1, 3E-1};
    vector<double> mut_params{1E-7, 3.16E-7, 1E-6, 3.16E-6, 1E-5, 3.16E-5, 1E-4, 3.16E-4, 1E-3};
    vector<double> move_params{1E-5, 3.16E-5, 1E-4, 3.16E-4, 1E-3, 3.16E-3, 1E-2, 3.16E-2, 1E-1, 3.16E-1, 1, 3.16, 10};
    vector<string> death_names{"1E1", "3E1"};
    vector<string> mut_names{"1E7", "3E7", "1E6", "3E6", "1E5", "3E5", "1E4", "3E4", "1E3"};
    vector<string> move_names{"1E5", "3E5", "1E4", "3E4", "1E3", "3E3", "1E2", "3E2", "1E1", "3E1", "1", "3", "10"};
    // Compute the indices
    int mov_num = move_params.size();
    int mut_num = mut_params.size();
    int death_index = task_id/(Repeat*mov_num*mut_num);
    int mut_index = (task_id % (mov_num*mut_num)) / mov_num;
    int move_index = (task_id % (mov_num*mut_num)) % mov_num;

    // Set parameters
	PProp[0].Move = move_params[move_index];
    PProp[0].DeathRate = death_params[death_index];
    PProp[0].MuteUp = mut_params[mut_index];
    string FileSpecifier = to_string(task_id)+"_M"+move_names[move_index]+"_MT"+mut_names[mut_index]+"_D"+death_names[death_index];
	// Run AdapImage
    AdapImage(LProp, PProp, AdapNum, FileSpecifier, MaxIterate);
	return 0;
}