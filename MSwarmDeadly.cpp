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
    vector<double> death_params{1e-1, 3e-1};
    vector<double> swarm_density_params{1E2, 2.15E2, 4.64E2, 1E3, 2.15E3, 4.64E3, 1E4, 2.15E4, 4.64E4, 6E4, 7E4, 8E4, 9E4};
    vector<double> swarm_move_params{3.16E-3, 1E-2, 3.16E-2, 1E-1, 3.16E-1, 1, 3.16, 10};
    vector<string> death_names{"1E1", "3E1"};
    vector<string> swarm_density_names{"1E2", "2E2", "4E2", "1E3", "2E3", "4E3", "1E4", "2E4", "4E4", "6E4", "7E4", "8E4", "9E4"};
    vector<string> swarm_move_names{"3E3", "1E2", "3E2", "1E1", "3E1", "1", "3", "10"};
    // Compute the indices
    int dense_num = swarm_density_params.size();
    int move_num = swarm_move_params.size();
    int death_index = task_id/(Repeat*move_num*dense_num);
    int dense_index = (task_id % (move_num*dense_num)) / move_num;
    int move_index = (task_id % (move_num*dense_num)) % move_num;
    int AdapNum = 7;

    // Set parameters
    LProp.GridBound = 8;
	LProp.GenBound = 8;
	LProp.PopulationNumber = 1;
    PProp[0].AutomaticInitiation = false;
    PProp[0].Move = 1E-3;
    PProp[0].DeathRate = death_params[death_index];
	PProp[0].ConsiderSwarm = true;
    PProp[0].SwarmDens = swarm_density_params[dense_index];
	PProp[0].SwarmMove = swarm_move_params[move_index];
    string FileSpecifier = to_string(task_id)+"_SwM"+swarm_move_names[move_index]+"_SwD"+swarm_density_names[dense_index]+"_D"+death_names[death_index];
	// Run AdapImage
    AdapImage(LProp, PProp, AdapNum, FileSpecifier, MaxIterate);
	return 0;
}