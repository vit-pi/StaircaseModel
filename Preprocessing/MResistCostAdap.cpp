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
    vector<double> cost_params{1E-4, 3.16E-4, 1E-3, 3.16E-3, 1E-2, 3.16E-2, 1E-1, 1.78E-1, 3.16E-1, 5.62E-1, 6.3E-1, 7.94E-1};
    vector<string> move_names{"1E2", "1"};
    vector<string> cost_names{"1E4", "3E4", "1E3", "3E3", "1E2", "3E2", "1E1", "1d78E1", "3E1", "5E1", "6E1", "7E1"};
    // Compute the indices
    int cost_num = cost_params.size();
    int move_index = task_id/(Repeat*cost_num);
    int cost_index = task_id % cost_num;

    // Set parameters
	PProp[0].Move = move_params[move_index];
    PProp[0].ResistCost = cost_params[cost_index];
    string FileSpecifier = to_string(task_id)+"_M"+move_names[move_index]+"_RC"+cost_names[cost_index];
	// Run AdapImage
    AdapImage(LProp, PProp, AdapNum, FileSpecifier, MaxIterate);
	return 0;
}