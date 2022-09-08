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
	double StopTime = 4E3;
    vector<double> switch_params{1E-3, 1E-2, 1E-1, 1, 10};
    vector<double> move0_params{1E-3, 1E-2, 1E-2, 1E-1, 1E-1, 1E-1, 1, 1, 1, 1, 10, 10, 10, 10, 10};
    vector<double> move1_params{1E-3, 1E-3, 1E-2, 1E-3, 1E-2, 1E-1, 1E-3, 1E-2, 1E-1, 1, 1E-3, 1E-2, 1E-1, 1, 10};
    vector<string> switch_names{"1E3", "1E2", "1E1", "1", "10"};
    vector<string> move0_names{"1E3", "1E2", "1E2", "1E1", "1E1", "1E1", "1", "1", "1", "1", "10", "10", "10", "10", "10"};
    vector<string> move1_names{"1E3", "1E3", "1E2", "1E3", "1E2", "1E1", "1E3", "1E2", "1E1", "1", "1E3", "1E2", "1E1", "1", "10"};
    PProp[0].SwitchUp = switch_params[(task_id/15)];
    PProp[1].SwitchUp = switch_params[(task_id/15)];
    PProp[0].SwitchDown = switch_params[(task_id/15)];
    PProp[1].SwitchDown = switch_params[(task_id/15)];
	PProp[0].Move = move0_params[(task_id%15)];
	PProp[1].Move = move1_params[(task_id%15)];
	LProp.GridBound = 8;
	LProp.GenBound = 8;
	LProp.PopulationNumber = 2;
	PProp[0].StressDepChemotax = false;
    string name = "M"+move0_names[(task_id%15)]+"_M"+move1_names[(task_id%15)]+"_S"+switch_names[(task_id/15)];
	// Run AnimTot simualtion
	AnimTot(LProp, PProp, 200, StopTime, name, 0);
	// system("pause");
	return 0;
}