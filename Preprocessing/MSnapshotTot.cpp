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

int main(void) {
	// Choose parameters
	double StopTime = 1E4;
    double SimInterval = 0.1;
    string FileSpecifier = "M1E3_SwM3_SwD5Ep4_D3E1";
    LProp.GridBound = 8;
	LProp.GenBound = 8;
	LProp.PopulationNumber = 1;
    PProp[0].InitCellsNum = 6E4;
    PProp[0].DeathRate = 3E-1;
	PProp[0].ConsiderSwarm = true;
    PProp[0].SwarmDens = 5E4;
    PProp[0].Move = 1E-3;
	PProp[0].SwarmMove = 3.16;
	// Run SnapshotTot simualtion
    SnapshotTot(LProp, PProp, SimInterval, StopTime, FileSpecifier);
	return 0;
}

// List of used simulations (write ParameterValue_... for FileSpecifier, StopTime, SimInterval):
// Notice values: 1e-3=1E3, 1e3=1Ep3 (p=positive, most exponents negative -> blank)

// M1E3_SwM3_SwD5Ep4_S10_D3E1_SwB1E3, 1E4, 0.1 - new density switch, slow-fast
// M1E3_SwM3_SwD5Ep4_S10_D3E1_SwB1E3, 1E4, 0.1 - new density switch, fast-slow
// M1E3_SwM3_SwD5Ep4_D3E1, 1E4, 0.1 - old density switch, slow-fast
// M1E3_SwM3_SwD5Ep4_3E1, 1E4, 0.1 - old density switch, fast-slow

// M1_C6E1, 1E4, 1 - positive chemotaxis
// M1_C2E1, 1E4, 1 - negative chemotaxis

// M3E1_C7E1, 1E4, 1 - positive chemotaxis
// M10_C3E1, 1E4, 1 - negative chemotaxis

// M2_M1E2_S1E3, 1E4, 1 - low switching
// M2_M1E2_S5, 1E4, 1 - high switching

// M10_D4E1_R0, 1E4, 1E-1 - deadly motility, R=0
// M10_D4E1_R1, 1E4, 1E-1 - deadly motility, R=1
// M10_D4E1_R2, 1E4, 1E-1 - deadly motility, R=2
// M10_D4E1_R3, 1E4, 1E-1 - deadly motility, R=3

// M1E2_RC2E1, 1E4, 1 - high resistance cost, low motility
// M1_RC2E1, 1E4, 1 - high resistance cost, high motility

// M1E2_HGT1E7, 2E3, 1 - low HGT, low motility
// M1E2_HGT1E4, 2E3, 1 - high HGT, high motility
// M1_HGT1E7, 2E3, 1 - low HGT, high motility
// M1_HGT1E4, 2E3, 1 - high HGT, high motility
// M1E2, 1E4, 1 - low motility
// M1, 1E4, 1 - high motility

// M1E2_M1E4_S1E3_D3E1, 1E4, 1 - low motility, low switching, high death
// M1E2_M1E4_S1E3_D1E1, 1E4, 1 - low motility, low switching, low death
// M1E2_M1E4_S5_D3E1, 1E4, 1 - low motility, high switching, high death
// M1E2_M1E4_S5_D1E1, 1E4, 1 - low motility, high switching, low death
// M3_M1E4_S1E3_D3E1, 1E4, 1 - mixed motility, low switching, high death
// M3_M1E4_S1E3_D1E1, 1E4, 1 - mixed motility, low switching, low death
// M3_M1E4_S5_D3E1, 1E4, 1 - mixed motility, high switching, high death
// M3_M1E4_S5_D1E1, 1E4, 1 - mixed motility, high switching, low death
// M3_M1_S1E3_D3E1, 1E4, 1 - high motility, low switching, high death
// M3_M1_S1E3_D1E1, 1E4, 1 - high motility, low switching, low death
// M3_M1_S5_D3E1, 1E4, 1 - high motility, high switching, high death
// M3_M1_S5_D1E1, 1E4, 1 - high motility, high switching, low death

// M1E2_C7E1, 1E4, 1 - low motility, positive chemotaxis
// M1_C7E1, 1E4, 1 - high motility, positive chemotaxis
// M1E2_C3E1, 1E4, 1 - low motility, negative chemotaxis
// M1_C3E1, 1E4, 1 - high motility, negative chemotaxis

// M1E2_SD7E1, 1E4, 1 - low motility, bactericidal (sigma=8)
// M1_SD7E1, 1E4, 1 - high motility, bactericidal (sigma=8)

// SwM1E2_SwD1Ep3_D1E1, 1E4, 0.1 - low swarm motility, low swarm density, low death
// SwM3_SwD1Ep3_D1E1, 1E4, 0.1 - high swarm motility, low swarm density, low death
// SwM1E2_SwD1Ep4_D1E1, 1E4, 0.1 - low swarm motility, high swarm density, low death
// SwM3_SwD1Ep4_D1E1, 1E4, 0.1 - high swarm motility, high swarm density, low death
// SwM1E2_SwD1Ep3_D3E1, 1E4, 0.1 - low swarm motility, low swarm density, high death
// SwM3_SwD1Ep3_D3E1, 1E4, 0.1 - high swarm motility, low swarm density, high death
// SwM1E2_SwD1Ep4_D3E1, 1E4, 0.1 - low swarm motility, high swarm density, high death
// SwM3_SwD1Ep4_D3E1, 1E4, 0.1 - high swarm motility, high swarm density, high death
// SwM1E2_SwD5Ep4_D1E1, 1E4, 0.1 - low swarm moitility, extra high swarm density, low death
// SwM1E2_SwD5Ep4_D3E1, 1E4, 0.1 - low swarm moitility, extra high swarm density, high death
// SwM3_SwD5Ep4_D1E1, 1E4, 0.1 - high swarm moitility, extra high swarm density, low death
// SwM3_SwD5Ep4_D3E1, 1E4, 0.1 - high swarm moitility, extra high swarm density, high death