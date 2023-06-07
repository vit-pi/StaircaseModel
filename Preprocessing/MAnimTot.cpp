///
/// INCLUSIONS
///

// Standard inclusions
#include <iostream>
#include <string>

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

int main(void) {
	// Choose parameters
	double StopTime = 1E2;
	PProp[0].Move = 1E-3;
	PProp[0].SwarmingMove = 1E-1;
	LProp.GridBound = 8;
	LProp.GenBound = 8;
	PProp[0].ConsiderSwarm = false;
	PProp[0].AutomaticInitiation = false;
	PProp[0].InitCellsNum = 9E4;
	// Run AnimTot simualtion
	AnimTot(LProp, PProp, 200, StopTime, "SwM1_SwD5E4_SwT3", 0);
	system("pause");
	return 0;
}