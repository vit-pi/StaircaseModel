# pragma once

using namespace std;

struct PopProp {
    // DECLARE STRAIN PROPERTIES (rates in 1/h units)
    int InitCellsNum = 1E2;
    int InitCellsPos = 0;
    int InitCellsGen = 0;
    bool AutomaticInitiation = false;
    double BirthRate = 1;
    double StressBirthRate = 0;
    double ResistCost = 0;
    double DeathRate = 1E-1;
    double StressDeathRate = 0;
    double MuteUp = 1E-7;
    double MuteDown = 1E-4;
    double StressMuteUp = 0;
    double StressMuteDown = 0;
    double Move = 1E-3;
    double Chemotax = 0.5; // in range [0, 1], probability to move right
    bool StressDepChemotax = false;
    double Chemokin = 0;
    double SwitchUp = 0;
    double SwitchDown = 0;
    bool ConsiderSwarm = false;
    double SwarmDens = 5E4;
    double SwarmMove = 1;
    bool ConsiderHGT = false;
    double HGTRate = 0;
};

struct LatProp {
    // LATTICE PROPERTIES
    int GridBound = 8;
    int GenBound = 8;
    int PopulationNumber = 1;
    double CarryingCapacity = 1E5;
    int D=1;
};
