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
    vector<double> chemotax_params{1E-1, 2E-1, 3E-1, 4E-1, 5E-1, 6E-1, 7E-1, 8E-1, 9E-1};
    vector<double> move_params{1E-5, 3.16E-5, 1E-4, 3.16E-4, 1E-3, 3.16E-3, 1E-2, 3.16E-2, 1E-1, 3.16E-1, 1, 3.16, 10};
    vector<string> chemotax_names{"1E1", "2E1", "3E1", "4E1", "5E1", "6E1", "7E1", "8E1", "9E1"};
    vector<string> move_names{"1E5", "3E5", "1E4", "3E4", "1E3", "3E3", "1E2", "3E2", "1E1", "3E1", "1", "3", "10"};
    // Compute the indices
    int mov_num = move_params.size();
    int chemotax_index = task_id/(Repeat*mov_num);
    int move_index = task_id % mov_num;

    // Set parameters
	PProp[0].Move = move_params[move_index];
    PProp[0].Chemotax = chemotax_params[chemotax_index];
    string FileSpecifier = to_string(task_id)+"_M"+move_names[move_index]+"_C"+chemotax_names[chemotax_index];
	// Run AdapImage
    AdapImage(LProp, PProp, AdapNum, FileSpecifier, MaxIterate);
	return 0;
}