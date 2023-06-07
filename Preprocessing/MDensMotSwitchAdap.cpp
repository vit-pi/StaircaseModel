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
    int switch_bias_index = 1;
    vector<double> switch_bias_params{0,1E-5/(1+1E-5),1E-4/(1+1E-4),1E-3/(1+1E-3),1E-2/(1+1E-2),1E-1/(1+1E-1)};
    vector<double> death_params{1E-1, 3E-1};
    vector<double> switch_density_params{1E2, 9.5E4};
    vector<double> move_params{1E-5, 3.16E-5, 1E-4, 3.16E-4, 1E-3, 3.16E-3, 1E-2, 3.16E-2, 1E-1, 3.16E-1, 1, 3.16, 10};
    vector<string> switch_bias_names{"0", "1E5", "1E4", "1E3", "1E2", "1E1"};
    vector<string> death_names{"1E1", "3E1"};
    vector<string> switch_density_names{"1E2", "9E4"};
    vector<string> move_names{"1E5", "3E5", "1E4", "3E4", "1E3", "3E3", "1E2", "3E2", "1E1", "3E1", "1", "3", "10"};
    // Compute the indices
    int dense_num = switch_density_params.size();
    int move_num = move_params.size();
    int death_index = task_id/(Repeat*move_num*move_num*dense_num);
    int dense_index = (task_id % (move_num*move_num*dense_num)) / (move_num*move_num);
    int move_low_dens_index = ((task_id % (move_num*move_num*dense_num)) % (move_num*move_num))/move_num;
    int move_high_dens_index = ((task_id % (move_num*move_num*dense_num)) % (move_num*move_num))%move_num;
    int AdapNum = 2;

    // Set parameters
    LProp.GridBound = 8;
	LProp.GenBound = 8;
	LProp.PopulationNumber = 2;
    PProp[0].InitCellsNum = 50;
    for(int population=0; population<LProp.PopulationNumber; population++){
        PProp[population].ConsiderDensSwitch = true;
        PProp[population].DeathRate = death_params[death_index];
        PProp[population].DensSwitchDens = switch_density_params[dense_index];
        PProp[population].DensSwitchTot = 10;
        PProp[population].DensSwitchBias = switch_bias_params[switch_bias_index];
    }
    PProp[0].Move = move_params[move_low_dens_index];
	PProp[1].Move = move_params[move_high_dens_index];
    string FileSpecifier = to_string(task_id)+"_M"+move_names[move_low_dens_index]+"_SwM"+move_names[move_high_dens_index]+"_SwD"+switch_density_names[dense_index]+"_D"+death_names[death_index]+"_SwB"+switch_bias_names[switch_bias_index];
	// Run AdapImage
    AdapImage(LProp, PProp, AdapNum, FileSpecifier, MaxIterate);
	return 0;
}
