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
    vector<double> sbirth_params{1E-3, 1E-2, 3.16E-2, 5.62E-2, 1E-1, 1.26E-1, 1.58E-1, 2E-1, 2.51E-1, 3.16E-1, 4E-1, 5.01E-1, 6.31E-1, 7.94E-1};
    vector<string> move_names{"1E2", "1"};
    vector<string> sbirth_names{"1E3", "1E2", "3E2", "5E2", "1E1", "1d26E1", "1d58E1", "2E1", "2d51E1", "3E1", "4E1", "5E1", "6E1", "7E1"};
    // Compute the indices
    int sbirth_num = sbirth_params.size();
    int move_index = task_id/(Repeat*sbirth_num);
    int sbirth_index = task_id % sbirth_num;

    // Set parameters
	PProp[0].Move = move_params[move_index];
    PProp[0].StressBirthRate = sbirth_params[sbirth_index];
    string FileSpecifier = to_string(task_id)+"_M"+move_names[move_index]+"_SB"+sbirth_names[sbirth_index];
	// Run AdapImage
    AdapImage(LProp, PProp, AdapNum, FileSpecifier, MaxIterate);
	return 0;
}
