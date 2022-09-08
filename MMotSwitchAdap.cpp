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
vector<PopProp> PProp(2);

///
/// MAIN CODE
///

int main(int argc, char **argv) {
    int task_id = atoi(argv[1]);
	// Choose parameters
    int Repeat = 10;
	long int MaxIterate = 1E11;
    int AdapNum = 7;
    vector<double> switch_params{0, 1E-3, 1E-2, 1E-1, 1};
    vector<double> move_params{1E-5, 3.16E-5, 1E-4, 3.16E-4, 1E-3, 3.16E-3, 1E-2, 3.16E-2, 1E-1, 3.16E-1, 1, 3.16, 10};
    vector<string> switch_names{"0", "1E3", "1E2", "1E1", "1"};
    vector<string> move_names{"1E5", "3E5", "1E4", "3E4", "1E3", "3E3", "1E2", "3E2", "1E1", "3E1", "1", "3", "10"};
    // Compute the indices
    int mov_num = move_params.size();
    int mov_options = mov_num*(mov_num+1)/2;
    int switch_index = task_id/(Repeat*mov_options);
    int x = task_id % mov_options;
    int move0_index = 0;
    while(x>=move0_index*(move0_index+1)/2){
        move0_index++;
    }
    move0_index--;
    int move1_index = x - move0_index*(move0_index+1)/2;

    // Set parameters
	LProp.PopulationNumber = 2;
    for (int population=0; population<LProp.PopulationNumber; population++) {
        PProp[population].SwitchUp = switch_params[switch_index];
        PProp[population].SwitchDown = switch_params[switch_index];
    }
	PProp[0].Move = move_params[move0_index];
	PProp[1].Move = move_params[move1_index];
    string FileSpecifier = to_string(task_id)+"_M"+move_names[move0_index]+"_M"+move_names[move1_index]+"_S"+switch_names[switch_index];
	// Run AdapImage
    AdapImage(LProp, PProp, AdapNum, FileSpecifier, MaxIterate);
	return 0;
}