#pragma once
#include "Properties.h"
#include "Lattice.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <ctime>
#include <random>
#include <limits>

using namespace std;

//
// DEFINE RANDOM ENGINE AND DISTRIBUTIONS
//

mt19937_64 RandEngine(time(0));
uniform_real_distribution<double> UnifDist(0, 1);

//
// INDEX PRIORITY QUEUE HANDLING FUNCTIONS
//

// Swap values at node1 and node2 and change the node array appropriately
void IndexPriorityQueue::Swap(int node1, int node2) {
    // Swap the nodes in the index array
	NodeArray[ActionQueue[node1]] = node2;
    NodeArray[ActionQueue[node2]] = node1;
    // Swap the values in the priority queue
	swap(TimeQueue[node1], TimeQueue[node2]);
    swap(ActionQueue[node1], ActionQueue[node2]);
}

// Heapify the tree rooted in node "node"
void IndexPriorityQueue::Heapify(int node) {
    int smallest = node; // Initialize smallest as root 
    int left_child = 2 * node + 1; // left = 2*i + 1 
    int right_child = 2 * node + 2; // right = 2*i + 2 

    // If left child is smaller than root 
    if (left_child < QueueSize && TimeQueue[left_child] < TimeQueue[smallest]) {
        smallest = left_child;
    }

    // If right child is smaller than smallest so far
    if (right_child < QueueSize && TimeQueue[right_child] < TimeQueue[smallest]) {
        smallest = right_child;
    }

    // If smallest is not root 
    if (smallest != node) {
        Swap(node, smallest);
        // Recursively heapify the affected sub-tree 
        Heapify(smallest);
    }
}

// Insert initial putative times for reactions in their order
void IndexPriorityQueue::Build(vector<double> delta_time) {
	TimeQueue = delta_time;
	QueueSize = delta_time.size();
	for (int i = 0; i < QueueSize; i++) {
		NodeArray.push_back(i);
        ActionQueue.push_back(i);
	}
    int start_node = (QueueSize/2)-1;
    for (int node = start_node; node >= 0; node--) {
        Heapify(node);
    }
}

// Update a given node
void IndexPriorityQueue::Update(int node, double new_delta_time) {
    // Update the putative time of a given node
    TimeQueue[node] = new_delta_time;
    // Reorganize the priority queue
    UpdateAux(node);
}

// Reorganize the priority queue
void IndexPriorityQueue::UpdateAux(int node) {
    int smallest_child = node;
    if (node < (QueueSize-1) / 2) {
        smallest_child = (TimeQueue[2 * node + 1]<TimeQueue[2 * node + 2]) ? 2 * node + 1 : 2 * node + 2;
    }
    else if (2 * node + 1 == QueueSize-1) {
        smallest_child = 2 * node + 1;
    }
    if ((node != 0) && (TimeQueue[node] < TimeQueue[(node - 1) / 2])) {
        Swap(node, (node - 1) / 2);
        UpdateAux((node - 1) / 2); 
    }
    else if ((node < QueueSize / 2) && (TimeQueue[node] > TimeQueue[smallest_child])) {
        Swap(node, smallest_child);
        UpdateAux(smallest_child);
    }
}
//
// NEWTON-RAPHSON - old version, do not use for non-standard staircase models
//

// Initialize the solver
NewtonSolve::NewtonSolve(double Iprecision, int IK, double Ir, double Is_r, double Id, double Is_d, double Inu_r, double Is_nu_r, double Inu_l, double Is_nu_l, int ID, int IGridBound){
    // Save variables
    precision=Iprecision;
    K=(double) IK;
    r=Ir;
    s_r=Is_r;
    d=Id;
    s_d=Is_d;
	nu_r=Inu_r;
	s_nu_r=Is_nu_r;
	nu_l=Inu_l;
	s_nu_l=Is_nu_l;
	D=ID;
	GridBound=IGridBound;
    // Prepare initial N and step
    step = vector<double>(GridBound, 2*precision);
    N = vector<double>(GridBound, 0);
    J = vector<double>(GridBound, 0);
    MinusF = vector<double>(GridBound, 0);
    nu_left = vector<double>(GridBound-1, s_nu_l);
    nu_right = vector<double>(GridBound-1, s_nu_r);
    for(int position=0; position<D; position++){
        N[position] = K;
        nu_left[position] = nu_l;
        nu_right[position] = nu_r;
    }
    nu_right[D-1]=s_nu_r;
}

// Find a magnitude of a vector in l1 metric
double NewtonSolve::Magnitude(vector<double> vec){
    double magnitude = 0;
    for (int i = 0; i<vec.size(); i++){
        magnitude += fabs(vec[i]);
    }
    return magnitude;
}

// Solve for the stable wild-type population
vector<double> NewtonSolve::Solve(){
    double w;
    // Approximate the root until the magnitude of the step is below some bound
    while(Magnitude(step)>=precision){
        // Compute N and J
        MinusF[0]=(nu_r+d-r*(1-N[0]/K))*N[0]-nu_l*N[1];
        J[0]=r-nu_r-d-2*r*N[0]/K;
        for (int position=1; position<D; position++){
            MinusF[position] = (nu_r+nu_l+d-r*(1-N[position]/K))*N[position]-nu_l*N[position+1]-nu_r*N[position-1];
            J[position] = r-nu_r-nu_l-d-2*r*N[position]/K;
        }
        for (int position=D; position<GridBound-1; position++){
            MinusF[position] = (s_nu_r+s_nu_l+s_d-s_r*(1-N[position]/K))*N[position]-s_nu_l*N[position+1]-s_nu_r*N[position-1];
            J[position] = s_r-s_nu_r-s_nu_l-s_d-2*s_r*N[position]/K;
        }
        MinusF[GridBound-1] = (s_nu_l+s_d-s_r*(1-N[GridBound-1]/K))*N[GridBound-1]-s_nu_r*N[GridBound-2];
        J[GridBound-1] = s_r-s_nu_l-s_d-2*s_r*N[GridBound-1]/K;
        // Compute step via tridiagonal matrix algorithm
        for(int position=1; position<GridBound; position++){
            w=nu_right[position-1]/J[position-1];
            J[position] -= w*nu_left[position-1];
            MinusF[position] -= w*MinusF[position-1];
        }
        step[GridBound-1]=MinusF[GridBound-1]/J[GridBound-1];
        for (int position = GridBound-2; position >= 0; position--){
            step[position] = (MinusF[position]-nu_left[position]*step[position+1])/J[position];
        }
        // Add the step to current N and repeat
        for (int position=0; position<GridBound; position++){
            N[position] += step[position];
        }
    }
    return N;
}

//
// LATTICE
//

// Compute stable wild-type population via bisection method - only for standard staircase
void Lattice::StablePopulation(double precision, int population){
    // Prepare variables
    double r = PProp[population].BirthRate;
	double s_r = PProp[population].StressBirthRate;
	double d = PProp[population].DeathRate+PProp[population].MuteUp;
	double s_d = PProp[population].DeathRate+PProp[population].MuteUp+PProp[population].StressDeathRate+PProp[population].StressMuteUp;
	double nu_r = PProp[population].Move;
	double s_nu_r = (PProp[population].Move+PProp[population].Chemokin)*2*PProp[population].Chemotax;
	double nu_l = PProp[population].Move;
	double s_nu_l = (PProp[population].Move+PProp[population].Chemokin)*2*(1-PProp[population].Chemotax);
    if (!PProp[population].StressDepChemotax){
        nu_r *= 2*PProp[population].Chemotax;
        nu_l *= 2*(1-PProp[population].Chemotax);
    }
    // Initialize solver and solve
    NewtonSolve solver(precision, LProp.CarryingCapacity, r, s_r, d, s_d, nu_r, s_nu_r, nu_l, s_nu_l, LProp.D, LProp.GridBound);
    vector<double> stable_pop = solver.Solve();
    // Write the result into StablePop
    for(int position=0;position<LProp.GridBound;position++){
        StablePop[population][position]=(int) stable_pop[position];
    }
}

// Find action rate
double Lattice::FindActionRate(int action) {
    double action_rate = 0;
    switch (Actions[action][3]) {
    case 0: // Divide
        if (AbsSpaceTot[Actions[action][2]]+(PProp[Actions[action][0]].CompetBelowStairs-1)*AbsSpaceTotUnderStair[Actions[action][2]] <= LProp.CarryingCapacity) {
            if (Actions[action][1] + LProp.D > Actions[action][2]) {
                action_rate = PProp[Actions[action][0]].BirthRate * pow(1 - PProp[Actions[action][0]].ResistCost, Actions[action][1] - Actions[action][2]+LProp.D-1);
            }
            else {
                action_rate = PProp[Actions[action][0]].StressBirthRate;
            }
            action_rate *= (1 - (AbsSpaceTot[Actions[action][2]]+(PProp[Actions[action][0]].CompetBelowStairs-1)*AbsSpaceTotUnderStair[Actions[action][2]]) / LProp.CarryingCapacity);
        }
        action_rate *= GenSpaceTot[Actions[action][0]][Actions[action][1]][Actions[action][2]];
        break;
    case 1: // Die
        action_rate = PProp[Actions[action][0]].DeathRate;
        if (Actions[action][1] + LProp.D - 1 < Actions[action][2]) {
            action_rate += PProp[Actions[action][0]].StressDeathRate;
        }
        action_rate *= GenSpaceTot[Actions[action][0]][Actions[action][1]][Actions[action][2]];
        break;
    case 2: // MuteUp
        action_rate = PProp[Actions[action][0]].MuteUp;
        if (Actions[action][1] + LProp.D - 1 < Actions[action][2]) {
            action_rate += PProp[Actions[action][0]].StressMuteUp;
        }
        action_rate *= GenSpaceTot[Actions[action][0]][Actions[action][1]][Actions[action][2]];
        break;
    case 3: // MuteDown
        action_rate = PProp[Actions[action][0]].MuteDown;
        if (Actions[action][1] + LProp.D - 1 < Actions[action][2]) {
            action_rate += PProp[Actions[action][0]].StressMuteDown;
        }
        action_rate *= GenSpaceTot[Actions[action][0]][Actions[action][1]][Actions[action][2]];
        break;
    case 4: // MoveRight
        // Compute usual or swarming motility
        if ((PProp[Actions[action][0]].ConsiderSwarm) && (SpaceTot[Actions[action][0]][Actions[action][2]] >= PProp[Actions[action][0]].SwarmDens)) {
            action_rate = PProp[Actions[action][0]].SwarmMove;
        }
        else {
            action_rate = PProp[Actions[action][0]].Move;
        }
        // Compute biased motility due to chemokinesis
        if (Actions[action][1] + LProp.D - 1 < Actions[action][2]) {
            action_rate += PProp[Actions[action][0]].Chemokin;
        }
        // Compute the biased motility due to chemotaxis
        if ((!PProp[Actions[action][0]].StressDepChemotax) || Actions[action][1] + LProp.D - 1 < Actions[action][2]) {
            action_rate *= 2*PProp[Actions[action][0]].Chemotax;
        }
        action_rate *= GenSpaceTot[Actions[action][0]][Actions[action][1]][Actions[action][2]];
        break;
    case 5: // MoveLeft
        // Compute usual or swarming motility
        if ((PProp[Actions[action][0]].ConsiderSwarm) && (SpaceTot[Actions[action][0]][Actions[action][2]] >= PProp[Actions[action][0]].SwarmDens)) {
            action_rate = PProp[Actions[action][0]].SwarmMove;
        }
        else {
            action_rate = PProp[Actions[action][0]].Move;
        }
        // Compute biased motility due to chemokinesis
        if (Actions[action][1] + LProp.D - 1 < Actions[action][2]) {
            action_rate += PProp[Actions[action][0]].Chemokin;
        }
        // Compute the biased motility due to chemotaxis (commented full version for later, uncomented positive=right, negative=left)
        if ((!PProp[Actions[action][0]].StressDepChemotax) || Actions[action][1] + LProp.D - 1 < Actions[action][2]) {
            action_rate *= 2*(1-PProp[Actions[action][0]].Chemotax);
        }
        action_rate *= GenSpaceTot[Actions[action][0]][Actions[action][1]][Actions[action][2]];
        break;
    case 6: // SwitchUp
        action_rate = PProp[Actions[action][0]].SwitchUp;
        if (PProp[Actions[action][0]].ConsiderDensSwitch){
            if (AbsSpaceTot[Actions[action][2]]<PProp[Actions[action][0]].DensSwitchDens){
                action_rate = 2*PProp[Actions[action][0]].DensSwitchBias*PProp[Actions[action][0]].DensSwitchTot;
            }
            else{
                action_rate = 2*(1-PProp[Actions[action][0]].DensSwitchBias)*PProp[Actions[action][0]].DensSwitchTot;
            }
        }
        action_rate *= GenSpaceTot[Actions[action][0]][Actions[action][1]][Actions[action][2]];
        break;
    case 7: // SwitchDown
        action_rate = PProp[Actions[action][0]].SwitchDown;
        if (PProp[Actions[action][0]].ConsiderDensSwitch){
            if (AbsSpaceTot[Actions[action][2]]<PProp[Actions[action][0]].DensSwitchDens){
                action_rate = 2*(1-PProp[Actions[action][0]].DensSwitchBias)*PProp[Actions[action][0]].DensSwitchTot;
            }
            else{
                action_rate = 2*PProp[Actions[action][0]].DensSwitchBias*PProp[Actions[action][0]].DensSwitchTot;
            }
        }
        action_rate *= GenSpaceTot[Actions[action][0]][Actions[action][1]][Actions[action][2]];
        break;
    case 8: // HGT
        action_rate = PProp[Actions[action][0]].HGTRate;
        action_rate *= GenSpaceTot[Actions[action][0]][Actions[action][4]][Actions[action][2]]; // donor
        action_rate *= GenSpaceTot[Actions[action][0]][Actions[action][1]][Actions[action][2]];
        break;
    }
    return action_rate;
}

// Execute action
void Lattice::ExecuteAction(int action) {
    switch (Actions[action][3]) {
    case 0: // Divide
        GenSpaceTot[Actions[action][0]][Actions[action][1]][Actions[action][2]] += 1;
        SpaceTot[Actions[action][0]][Actions[action][2]] += 1;
        GenTot[Actions[action][0]][Actions[action][1]] += 1;
        AbsSpaceTot[Actions[action][2]] += 1;
        if (Actions[action][1] + LProp.D - 1 < Actions[action][2]){
            AbsSpaceTotUnderStair[Actions[action][2]] += 1;
        }
        AbsGenTot[Actions[action][1]] += 1;
        AbsTot += 1;
        break;
    case 1: // Die
        GenSpaceTot[Actions[action][0]][Actions[action][1]][Actions[action][2]] -= 1;
        SpaceTot[Actions[action][0]][Actions[action][2]] -= 1;
        GenTot[Actions[action][0]][Actions[action][1]] -= 1;
        AbsSpaceTot[Actions[action][2]] -= 1;
        if (Actions[action][1] + LProp.D - 1 < Actions[action][2]){
            AbsSpaceTotUnderStair[Actions[action][2]] -= 1;
        }
        AbsGenTot[Actions[action][1]] -= 1;
        AbsTot -= 1;
        break;
    case 2: // MuteUp
        GenSpaceTot[Actions[action][0]][Actions[action][1]][Actions[action][2]] -= 1;
        GenTot[Actions[action][0]][Actions[action][1]] -= 1;
        AbsGenTot[Actions[action][1]] -= 1;
        if (Actions[action][1] + LProp.D == Actions[action][2]){
            AbsSpaceTotUnderStair[Actions[action][2]] -= 1;
        }

        GenSpaceTot[Actions[action][0]][Actions[action][1] + 1][Actions[action][2]] += 1;
        GenTot[Actions[action][0]][Actions[action][1] + 1] += 1;
        AbsGenTot[Actions[action][1] + 1] += 1;
        break;
    case 3: // MuteDown
        GenSpaceTot[Actions[action][0]][Actions[action][1]][Actions[action][2]] -= 1;
        GenTot[Actions[action][0]][Actions[action][1]] -= 1;
        AbsGenTot[Actions[action][1]] -= 1;
        if (Actions[action][1] + LProp.D - 1 == Actions[action][2]){
            AbsSpaceTotUnderStair[Actions[action][2]] += 1;
        }

        GenSpaceTot[Actions[action][0]][Actions[action][1] - 1][Actions[action][2]] += 1;
        GenTot[Actions[action][0]][Actions[action][1] - 1] += 1;
        AbsGenTot[Actions[action][1] - 1] += 1;
        break;
    case 4: // MoveRight
        GenSpaceTot[Actions[action][0]][Actions[action][1]][Actions[action][2]] -= 1;
        SpaceTot[Actions[action][0]][Actions[action][2]] -= 1;
        AbsSpaceTot[Actions[action][2]] -= 1;
        if (Actions[action][1] + LProp.D - 1 < Actions[action][2]){
            AbsSpaceTotUnderStair[Actions[action][2]] -= 1;
        }

        GenSpaceTot[Actions[action][0]][Actions[action][1]][Actions[action][2] + 1] += 1;
        SpaceTot[Actions[action][0]][Actions[action][2] + 1] += 1;
        AbsSpaceTot[Actions[action][2] + 1] += 1;
        if (Actions[action][1] + LProp.D - 1 < Actions[action][2] + 1){
            AbsSpaceTotUnderStair[Actions[action][2]+1] += 1;
        }
        break;
    case 5: // MoveLeft
        GenSpaceTot[Actions[action][0]][Actions[action][1]][Actions[action][2]] -= 1;
        SpaceTot[Actions[action][0]][Actions[action][2]] -= 1;
        AbsSpaceTot[Actions[action][2]] -= 1;
        if (Actions[action][1] + LProp.D -1 < Actions[action][2]){
            AbsSpaceTotUnderStair[Actions[action][2]] -= 1;
        }

        GenSpaceTot[Actions[action][0]][Actions[action][1]][Actions[action][2] - 1] += 1;
        SpaceTot[Actions[action][0]][Actions[action][2] - 1] += 1;
        AbsSpaceTot[Actions[action][2] - 1] += 1;
        if (Actions[action][1] + LProp.D - 1 < Actions[action][2] - 1){
            AbsSpaceTotUnderStair[Actions[action][2]-1] += 1;
        }
        break;
    case 6: // SwitchUp
        GenSpaceTot[Actions[action][0]][Actions[action][1]][Actions[action][2]] -= 1;
        SpaceTot[Actions[action][0]][Actions[action][2]] -= 1;
        GenTot[Actions[action][0]][Actions[action][1]] -= 1;

        GenSpaceTot[Actions[action][0] + 1][Actions[action][1]][Actions[action][2]] += 1;
        SpaceTot[Actions[action][0] + 1][Actions[action][2]] += 1;
        GenTot[Actions[action][0] + 1][Actions[action][1] + 1] += 1;
        break;
    case 7: // SwitchDown
        GenSpaceTot[Actions[action][0]][Actions[action][1]][Actions[action][2]] -= 1;
        SpaceTot[Actions[action][0]][Actions[action][2]] -= 1;
        GenTot[Actions[action][0]][Actions[action][1]] -= 1;

        GenSpaceTot[Actions[action][0] - 1][Actions[action][1]][Actions[action][2]] += 1;
        SpaceTot[Actions[action][0] - 1][Actions[action][2]] += 1;
        GenTot[Actions[action][0] - 1][Actions[action][1] + 1] += 1;
        break;
    case 8: // HGT
        GenSpaceTot[Actions[action][0]][Actions[action][1]][Actions[action][2]] -= 1;
        GenTot[Actions[action][0]][Actions[action][1]] -= 1;
        AbsGenTot[Actions[action][1]] -= 1;
        if ((Actions[action][4] + LProp.D > Actions[action][2]) && (Actions[action][1] + LProp.D <= Actions[action][2])){
            AbsSpaceTotUnderStair[Actions[action][2]] -= 1;
        }

        GenSpaceTot[Actions[action][0]][Actions[action][4]][Actions[action][2]] += 1;
        GenTot[Actions[action][0]][Actions[action][4]] += 1;
        AbsGenTot[Actions[action][4]] += 1;
        break;
    }
}

// Initiate the lattice object
Lattice::Lattice(LatProp l_prop, vector<PopProp> p_prop) {
    // Save properties and rescale appropriately whatever is needed
    LProp = l_prop;
    PProp = p_prop;

    // Initiate time
    Time = 0;

    // Initiate Cell Totals
    // Empty totals with 0s
    StablePop = vector<vector<int>>(LProp.PopulationNumber, vector<int>(LProp.GridBound,0));
    AbsTot = 0;
    AbsSpaceTot = vector<int>(LProp.GridBound, 0);
    AbsSpaceTotUnderStair = vector<int>(LProp.GridBound, 0);
    AbsGenTot = vector<int>(LProp.GenBound, 0);
    for (int population = 0; population < LProp.PopulationNumber; population++) {
        GenSpaceTot.push_back(vector<vector<int>>(LProp.GenBound, vector<int>(LProp.GridBound, 0)));
        GenTot.push_back(vector<int>(LProp.GenBound, 0));
        SpaceTot.push_back(vector<int>(LProp.GridBound, 0));
    }

    // Fill totals
    for (int population = 0; population < LProp.PopulationNumber; population++) {
        if (PProp[population].AutomaticInitiation){
            StablePopulation(1e-2, population);
            GenSpaceTot[population][0]=StablePop[population];
            SpaceTot[population]=StablePop[population];
            int stable_tot=0;
            for (int position = 0; position < LProp.GridBound; position ++){
                stable_tot += StablePop[population][position];
                AbsSpaceTot[position] += StablePop[population][position];
                if (position >=  LProp.D){
                    AbsSpaceTotUnderStair[position] += StablePop[population][position];
                }
            }
            GenTot[population][0] = stable_tot;
            AbsGenTot[0] += stable_tot;
            AbsTot += stable_tot;
        }
        else{
            GenSpaceTot[population][PProp[population].InitCellsGen][PProp[population].InitCellsPos] = PProp[population].InitCellsNum;
            SpaceTot[population][PProp[population].InitCellsPos] = PProp[population].InitCellsNum;
            GenTot[population][PProp[population].InitCellsGen] = PProp[population].InitCellsNum;
            AbsSpaceTot[PProp[population].InitCellsPos] += PProp[population].InitCellsNum;
            AbsGenTot[PProp[population].InitCellsGen] += PProp[population].InitCellsNum;
            AbsTot += PProp[population].InitCellsNum;
            if (PProp[population].InitCellsGen + LProp.D <= PProp[population].InitCellsPos){
                AbsSpaceTotUnderStair[PProp[population].InitCellsPos] += PProp[population].InitCellsNum;
            }
        } 
    }

    // Initiate Actions
    for (int position = 0; position < LProp.GridBound; position++) {
        for (int population = 0; population < LProp.PopulationNumber; population++) {
            for (int genotype = 0; genotype < LProp.GenBound; genotype++) {
                // 0 = Divide
                    Actions.push_back(vector<int> {population, genotype, position, 0});
                // 1 = Die
                    Actions.push_back(vector<int> {population, genotype, position, 1});
                // 2 = MuteUp
                    if (genotype < LProp.GenBound-1) {
                        Actions.push_back(vector<int> {population, genotype, position, 2});
                    }
                // 3 = MuteDown
                    if (genotype > 0) {
                        Actions.push_back(vector<int> {population, genotype, position, 3});
                    }
                // 4 = MoveRight
                    if (position < LProp.GridBound - 1) {
                        Actions.push_back(vector<int> {population, genotype, position, 4});
                    }
                // 5 = MoveLeft
                    if (position > 0) {
                        Actions.push_back(vector<int> {population, genotype, position, 5});
                    }
                // 6 = SwitchUp
                    if (population < LProp.PopulationNumber-1) {
                        Actions.push_back(vector<int> {population, genotype, position, 6});
                    }
                // 7 = SwitchDown
                    if (population > 0) {
                        Actions.push_back(vector<int> {population, genotype, position, 7});
                    }
                // 8 = HGT
                    if (PProp[population].ConsiderHGT) {
                        if (genotype < LProp.GenBound-1) {
                            for (int genotype_donor = genotype + 1;  genotype_donor < LProp.GenBound; genotype_donor++) {
                                Actions.push_back(vector<int> {population, genotype, position, 8, genotype_donor});
                            }
                        }
                    }
            }
        }
    }

    // Build Dependency Graph
    vector<int> container;
    for (int action = 0; action < Actions.size(); action++) {
        container.clear();
        switch (Actions[action][3])
        {
        case 0: // Divide
            for (int i = action; (i > -1) && (Actions[i][2] > Actions[action][2] - 1); i--) {
                if (Actions[i][3] == 0) {
                    container.push_back(i);
                }
                else if ((Actions[i][1] == Actions[action][1]) && (Actions[i][0] == Actions[action][0]) && (Actions[i][2] == Actions[action][2])) {
                    container.push_back(i);
                }
                else if (PProp[Actions[i][0]].ConsiderSwarm && (Actions[i][0] == Actions[action][0]) && ((Actions[i][3] == 4)||(Actions[i][3] == 5))) {
                    container.push_back(i);
                }
                else if (PProp[Actions[i][0]].ConsiderHGT && (Actions[i][0] == Actions[action][0]) && (Actions[i][2] == Actions[action][2]) && (Actions[i][3] == 8) && (Actions[i][4] == Actions[action][1])) {
                    container.push_back(i);
                }
                else if (PProp[Actions[i][0]].ConsiderDensSwitch && (Actions[i][2] == Actions[action][2]) && ((Actions[i][3] == 6)||(Actions[i][3] == 7))){
                    container.push_back(i);
                }
            }
            for (int i = action + 1; (i < Actions.size()) && (Actions[i][2] < Actions[action][2] + 1); i++) {
                if (Actions[i][3] == 0) {
                    container.push_back(i);
                }
                else if ((Actions[i][1] == Actions[action][1]) && (Actions[i][0] == Actions[action][0]) && (Actions[i][2] == Actions[action][2])) {
                    container.push_back(i);
                }
                else if (PProp[Actions[i][0]].ConsiderSwarm && (Actions[i][0] == Actions[action][0]) && ((Actions[i][3] == 4)||(Actions[i][3] == 5))) {
                    container.push_back(i);
                }
                else if (PProp[Actions[i][0]].ConsiderHGT && (Actions[i][0] == Actions[action][0]) && (Actions[i][2] == Actions[action][2]) && (Actions[i][3] == 8) && (Actions[i][4] == Actions[action][1])) {
                    container.push_back(i);
                }
                else if ((PProp[Actions[i][0]].ConsiderDensSwitch) && (Actions[i][2] == Actions[action][2]) && ((Actions[i][3] == 6)||(Actions[i][3] == 7))){
                    container.push_back(i);
                }
            }
            break;
        case 1: // Die
            for (int i = action; (i > -1) && (Actions[i][2] > Actions[action][2] - 1); i--) {
                if (Actions[i][3] == 0) {
                    container.push_back(i);
                }
                else if ((Actions[i][1] == Actions[action][1]) && (Actions[i][0] == Actions[action][0]) && (Actions[i][2] == Actions[action][2])) {
                    container.push_back(i);
                }
                else if (PProp[Actions[i][0]].ConsiderSwarm && (Actions[i][0] == Actions[action][0]) && ((Actions[i][3] == 4)||(Actions[i][3] == 5))) {
                    container.push_back(i);
                }
                else if (PProp[Actions[i][0]].ConsiderHGT && (Actions[i][0] == Actions[action][0]) && (Actions[i][2] == Actions[action][2]) && (Actions[i][3] == 8) && (Actions[i][4] == Actions[action][1])) {
                    container.push_back(i);
                }
                else if (PProp[Actions[i][0]].ConsiderDensSwitch && (Actions[i][2] == Actions[action][2]) && ((Actions[i][3] == 6)||(Actions[i][3] == 7))){
                    container.push_back(i);
                }
            }
            for (int i = action + 1; (i < Actions.size()) && (Actions[i][2] < Actions[action][2] + 1); i++) {
                if (Actions[i][3] == 0) {
                    container.push_back(i);
                }
                else if ((Actions[i][1] == Actions[action][1]) && (Actions[i][0] == Actions[action][0]) && (Actions[i][2] == Actions[action][2])) {
                    container.push_back(i);
                }
                else if (PProp[Actions[i][0]].ConsiderSwarm && (Actions[i][0] == Actions[action][0]) && ((Actions[i][3] == 4)||(Actions[i][3] == 5))) {
                    container.push_back(i);
                }
                else if (PProp[Actions[i][0]].ConsiderHGT && (Actions[i][0] == Actions[action][0]) && (Actions[i][2] == Actions[action][2]) && (Actions[i][3] == 8) && (Actions[i][4] == Actions[action][1])) {
                    container.push_back(i);
                }
                else if (PProp[Actions[i][0]].ConsiderDensSwitch && (Actions[i][2] == Actions[action][2]) && ((Actions[i][3] == 6)||(Actions[i][3] == 7))){
                    container.push_back(i);
                }
            }
            break;
        case 2: // MuteUp
            for (int i = action; (i > -1) && (Actions[i][2] > Actions[action][2] - 1); i--) {
                if ((Actions[i][1] == Actions[action][1]) && (Actions[i][0] == Actions[action][0]) && (Actions[i][2] == Actions[action][2])) {
                    container.push_back(i);
                }
                else if (PProp[Actions[i][0]].ConsiderHGT && (Actions[i][0] == Actions[action][0]) && (Actions[i][2] == Actions[action][2]) && (Actions[i][3] == 8) && ((Actions[i][4] == Actions[action][1])||(Actions[i][4] == Actions[action][1]+1))) {
                    container.push_back(i);
                }
                else if ((Actions[action][1] + LProp.D == Actions[action][2]) && (Actions[i][2] == Actions[action][2]) && (Actions[i][3] == 0)){
                    container.push_back(i);
                }
            }
            for (int i = action + 1; (i < Actions.size()) && (Actions[i][2] < Actions[action][2] + 1); i++) {
                if ((Actions[i][1] == Actions[action][1]) && (Actions[i][0] == Actions[action][0]) && (Actions[i][2] == Actions[action][2])) {
                    container.push_back(i);
                }
                else if ((Actions[i][1] == Actions[action][1] + 1) && (Actions[i][0] == Actions[action][0]) && (Actions[i][2] == Actions[action][2])) {
                    container.push_back(i);
                }
                else if (PProp[Actions[i][0]].ConsiderHGT && (Actions[i][0] == Actions[action][0]) && (Actions[i][2] == Actions[action][2]) && (Actions[i][3] == 8) && ((Actions[i][4] == Actions[action][1])||(Actions[i][4] == Actions[action][1]+1))) {
                    container.push_back(i);
                }
                else if ((Actions[action][1] + LProp.D == Actions[action][2]) && (Actions[i][2] == Actions[action][2]) && (Actions[i][3] == 0)){
                    container.push_back(i);
                }
            }
            break;
        case 3: // MuteDown
            for (int i = action; (i > -1) && (Actions[i][2] > Actions[action][2] - 1); i--) {
                if ((Actions[i][1] == Actions[action][1]) && (Actions[i][0] == Actions[action][0]) && (Actions[i][2] == Actions[action][2])) {
                    container.push_back(i);
                }
                else if ((Actions[i][1] == Actions[action][1] - 1) && (Actions[i][0] == Actions[action][0]) && (Actions[i][2] == Actions[action][2])) {
                    container.push_back(i);
                }
                else if (PProp[Actions[i][0]].ConsiderHGT && (Actions[i][0] == Actions[action][0]) && (Actions[i][2] == Actions[action][2]) && (Actions[i][3] == 8) && ((Actions[i][4] == Actions[action][1])||(Actions[i][4] == Actions[action][1]-1))) {
                    container.push_back(i);
                }
                else if ((Actions[action][1] + LProp.D == Actions[action][2] + 1) && (Actions[i][2] == Actions[action][2]) && (Actions[i][3] == 0)){
                    container.push_back(i);
                }
            }
            for (int i = action + 1; (i < Actions.size()) && (Actions[i][2] < Actions[action][2] + 1); i++) {
                if ((Actions[i][1] == Actions[action][1]) && (Actions[i][0] == Actions[action][0]) && (Actions[i][2] == Actions[action][2])) {
                    container.push_back(i);
                }
                else if (PProp[Actions[i][0]].ConsiderHGT && (Actions[i][0] == Actions[action][0]) && (Actions[i][2] == Actions[action][2]) && (Actions[i][3] == 8) && ((Actions[i][4] == Actions[action][1])||(Actions[i][4] == Actions[action][1]-1))) {
                    container.push_back(i);
                }
                else if ((Actions[action][1] + LProp.D == Actions[action][2] + 1) && (Actions[i][2] == Actions[action][2]) && (Actions[i][3] == 0)){
                    container.push_back(i);
                }
            }
            break;
        case 4: // MoveRight
            for (int i = action; (i > -1) && (Actions[i][2] > Actions[action][2] - 1); i--) {
                if (Actions[i][3] == 0) {
                    container.push_back(i);
                }
                else if ((Actions[i][1] == Actions[action][1]) && (Actions[i][0] == Actions[action][0]) && (Actions[i][2] == Actions[action][2])) {
                    container.push_back(i);
                }
                else if (PProp[Actions[i][0]].ConsiderSwarm && (Actions[i][0] == Actions[action][0]) && ((Actions[i][3] == 4)||(Actions[i][3] == 5))) {
                    container.push_back(i);
                }
                else if (PProp[Actions[i][0]].ConsiderHGT && (Actions[i][0] == Actions[action][0]) && ((Actions[i][2] == Actions[action][2])||(Actions[i][2] == Actions[action][2]+1)) && (Actions[i][3] == 8) && (Actions[i][4] == Actions[action][1])) {
                    container.push_back(i);
                }
                else if (PProp[Actions[i][0]].ConsiderDensSwitch && ((Actions[i][2] == Actions[action][2])||(Actions[i][2] == Actions[action][2]+1)) && ((Actions[i][3] == 6)||(Actions[i][3] == 7))){
                    container.push_back(i);
                }
            }
            for (int i = action + 1; (i < Actions.size()) && (Actions[i][2] < Actions[action][2] + 2); i++) {
                if (Actions[i][3] == 0) {
                    container.push_back(i);
                }
                else if ((Actions[i][1] == Actions[action][1]) && (Actions[i][0] == Actions[action][0]) && (Actions[i][2] == Actions[action][2])) {
                    container.push_back(i);
                }
                else if ((Actions[i][1] == Actions[action][1]) && (Actions[i][0] == Actions[action][0]) && (Actions[i][2] == Actions[action][2] + 1)) {
                    container.push_back(i);
                }
                else if (PProp[Actions[i][0]].ConsiderSwarm && (Actions[i][0] == Actions[action][0]) && ((Actions[i][3] == 4)||(Actions[i][3] == 5))) {
                    container.push_back(i);
                }
                else if (PProp[Actions[i][0]].ConsiderHGT && (Actions[i][0] == Actions[action][0]) && ((Actions[i][2] == Actions[action][2])||(Actions[i][2] == Actions[action][2]+1)) && (Actions[i][3] == 8) && (Actions[i][4] == Actions[action][1])) {
                    container.push_back(i);
                }
                else if (PProp[Actions[i][0]].ConsiderDensSwitch && ((Actions[i][2] == Actions[action][2])||(Actions[i][2] == Actions[action][2]+1)) && ((Actions[i][3] == 6)||(Actions[i][3] == 7))){
                    container.push_back(i);
                }
            }
            break;
        case 5: // MoveLeft
            for (int i = action; (i > -1) && (Actions[i][2] > Actions[action][2] - 2); i--) {
                if (Actions[i][3] == 0) {
                    container.push_back(i);
                }
                else if ((Actions[i][1] == Actions[action][1]) && (Actions[i][0] == Actions[action][0]) && (Actions[i][2] == Actions[action][2])) {
                    container.push_back(i);
                }
                else if ((Actions[i][1] == Actions[action][1]) && (Actions[i][0] == Actions[action][0]) && (Actions[i][2] == Actions[action][2] - 1)) {
                    container.push_back(i);
                }
                else if (PProp[Actions[i][0]].ConsiderSwarm && (Actions[i][0] == Actions[action][0]) && ((Actions[i][3] == 4)||(Actions[i][3] == 5))) {
                    container.push_back(i);
                }
                else if (PProp[Actions[i][0]].ConsiderHGT && (Actions[i][0] == Actions[action][0]) && ((Actions[i][2] == Actions[action][2])||(Actions[i][2] == Actions[action][2]-1)) && (Actions[i][3] == 8) && (Actions[i][4] == Actions[action][1])) {
                    container.push_back(i);
                }
                else if (PProp[Actions[i][0]].ConsiderDensSwitch && ((Actions[i][2] == Actions[action][2])||(Actions[i][2] == Actions[action][2]-1)) && ((Actions[i][3] == 6)||(Actions[i][3] == 7))){
                    container.push_back(i);
                }
            }
            for (int i = action + 1; (i < Actions.size()) && (Actions[i][2] < Actions[action][2] + 1); i++) {
                if (Actions[i][3] == 0) {
                    container.push_back(i);
                }
                else if ((Actions[i][1] == Actions[action][1]) && (Actions[i][0] == Actions[action][0]) && (Actions[i][2] == Actions[action][2])) {
                    container.push_back(i);
                }
                else if (PProp[Actions[i][0]].ConsiderSwarm && (Actions[i][0] == Actions[action][0]) && ((Actions[i][3] == 4)||(Actions[i][3] == 5))) {
                    container.push_back(i);
                }
                else if (PProp[Actions[i][0]].ConsiderHGT && (Actions[i][0] == Actions[action][0]) && ((Actions[i][2] == Actions[action][2])||(Actions[i][2] == Actions[action][2]-1)) && (Actions[i][3] == 8) && (Actions[i][4] == Actions[action][1])) {
                    container.push_back(i);
                }
                else if (PProp[Actions[i][0]].ConsiderDensSwitch && ((Actions[i][2] == Actions[action][2])||(Actions[i][2] == Actions[action][2]-1)) && ((Actions[i][3] == 6)||(Actions[i][3] == 7))){
                    container.push_back(i);
                }
            }
            break;
        case 6: // SwitchUp
            for (int i = action; (i > -1) && (Actions[i][2] > Actions[action][2] - 1); i--) {
                if ((Actions[i][1] == Actions[action][1]) && (Actions[i][0] == Actions[action][0]) && (Actions[i][2] == Actions[action][2])) {
                    container.push_back(i);
                }
                else if (PProp[Actions[i][0]].ConsiderSwarm && (Actions[i][0] == Actions[action][0]) && ((Actions[i][3] == 4)||(Actions[i][3] == 5))) {
                    container.push_back(i);
                }
                else if (PProp[Actions[i][0]].ConsiderHGT && ((Actions[i][0] == Actions[action][0])||(Actions[i][0] == Actions[action][0]+1)) && (Actions[i][2] == Actions[action][2])  && (Actions[i][3] == 8) && (Actions[i][4] == Actions[action][1])) {
                    container.push_back(i);
                }
            }
            for (int i = action + 1; (i < Actions.size()) && (Actions[i][2] < Actions[action][2] + 1); i++) {
                if ((Actions[i][1] == Actions[action][1]) && (Actions[i][0] == Actions[action][0]) && (Actions[i][2] == Actions[action][2])) {
                    container.push_back(i);
                }
                else if ((Actions[i][1] == Actions[action][1]) && (Actions[i][0] == Actions[action][0] + 1) && (Actions[i][2] == Actions[action][2])) {
                    container.push_back(i);
                }
                else if (PProp[Actions[i][0]].ConsiderSwarm && ((Actions[i][0] == Actions[action][0])||(Actions[i][0] == Actions[action][0] + 1)) && ((Actions[i][3] == 4)||(Actions[i][3] == 5))) {
                    container.push_back(i);
                }
                else if (PProp[Actions[i][0]].ConsiderHGT && ((Actions[i][0] == Actions[action][0])||(Actions[i][0] == Actions[action][0]+1)) && (Actions[i][2] == Actions[action][2])  && (Actions[i][3] == 8) && (Actions[i][4] == Actions[action][1])) {
                    container.push_back(i);
                }
            }
            break;
        case 7: // SwitchDown
            for (int i = action; (i > -1) && (Actions[i][2] > Actions[action][2] - 1); i--) {
                if ((Actions[i][1] == Actions[action][1]) && (Actions[i][0] == Actions[action][0]) && (Actions[i][2] == Actions[action][2])) {
                    container.push_back(i);
                }
                else if ((Actions[i][1] == Actions[action][1]) && (Actions[i][0] == Actions[action][0] - 1) && (Actions[i][2] == Actions[action][2])) {
                    container.push_back(i);
                }
                else if (PProp[Actions[i][0]].ConsiderSwarm && ((Actions[i][0] == Actions[action][0])||(Actions[i][0] == Actions[action][0] - 1)) && ((Actions[i][3] == 4)||(Actions[i][3] == 5))) {
                    container.push_back(i);
                }
                else if (PProp[Actions[i][0]].ConsiderHGT && ((Actions[i][0] == Actions[action][0])||(Actions[i][0] == Actions[action][0]-1)) && (Actions[i][2] == Actions[action][2])  && (Actions[i][3] == 8) && (Actions[i][4] == Actions[action][1])) {
                    container.push_back(i);
                }
            }
            for (int i = action + 1; (i < Actions.size()) && (Actions[i][2] < Actions[action][2] + 1); i++) {
                if ((Actions[i][1] == Actions[action][1]) && (Actions[i][0] == Actions[action][0]) && (Actions[i][2] == Actions[action][2])) {
                    container.push_back(i);
                }
                else if (PProp[Actions[i][0]].ConsiderSwarm && (Actions[i][0] == Actions[action][0]) && ((Actions[i][3] == 4)||(Actions[i][3] == 5))) {
                    container.push_back(i);
                }
                else if (PProp[Actions[i][0]].ConsiderHGT && ((Actions[i][0] == Actions[action][0])||(Actions[i][0] == Actions[action][0]-1)) && (Actions[i][2] == Actions[action][2])  && (Actions[i][3] == 8) && (Actions[i][4] == Actions[action][1])) {
                    container.push_back(i);
                }
            }
            break;
        case 8: // HGT
            for (int i = action; (i > -1) && (Actions[i][2] > Actions[action][2] - 1); i--) {
                if ((Actions[i][1] == Actions[action][1]) && (Actions[i][0] == Actions[action][0]) && (Actions[i][2] == Actions[action][2])) {
                    container.push_back(i);
                }
                else if (PProp[Actions[i][0]].ConsiderHGT && (Actions[i][0] == Actions[action][0]) && (Actions[i][2] == Actions[action][2]) && (Actions[i][3] == 8) && ((Actions[i][4] == Actions[action][1])||(Actions[i][4] == Actions[action][1]+1))) {
                    container.push_back(i);
                }
                else if ((Actions[action][1] + LProp.D <= Actions[action][2])&&(Actions[action][4] + LProp.D > Actions[action][2])&& (Actions[i][2] == Actions[action][2]) && (Actions[i][3] == 0)){
                    container.push_back(i);
                }
            }
            for (int i = action + 1; (i < Actions.size()) && (Actions[i][2] < Actions[action][2] + 1); i++) {
                if ((Actions[i][1] == Actions[action][1]) && (Actions[i][0] == Actions[action][0]) && (Actions[i][2] == Actions[action][2])) {
                    container.push_back(i);
                }
                else if ((Actions[i][1] == Actions[action][4]) && (Actions[i][0] == Actions[action][0]) && (Actions[i][2] == Actions[action][2])) {
                    container.push_back(i);
                }
                else if (PProp[Actions[i][0]].ConsiderHGT && (Actions[i][0] == Actions[action][0]) && (Actions[i][2] == Actions[action][2]) && (Actions[i][3] == 8) && ((Actions[i][4] == Actions[action][1])||(Actions[i][4] == Actions[action][1]+1))) {
                    container.push_back(i);
                }
                else if ((Actions[action][1] + LProp.D <= Actions[action][2])&&(Actions[action][4] + LProp.D > Actions[action][2])&& (Actions[i][2] == Actions[action][2]) && (Actions[i][3] == 0)){
                    container.push_back(i);
                }
            }
            break;
        }
        DependencyGraph.push_back(container);
    }

    // Calculate Action Rates
    for (int action = 0; action < Actions.size(); action++) {
        ActionRates.push_back(FindActionRate(action));
    }

    // Calculate Putative Times and Organize Them into PutativeTimesQueue
    vector<double> putative_times = vector<double>(Actions.size(), 0);
    PreInfiniteTime = vector<double>(Actions.size(), 0);
    for (int action = 0; action < Actions.size(); action++) {
        if (ActionRates[action]==0){
            putative_times[action]=numeric_limits<double>::infinity();
        }
        else{
            putative_times[action]=-log(UnifDist(RandEngine))/ActionRates[action];
        }
    }
    PutativeTimesQueue.Build(putative_times);
}

// Update lattice object
void Lattice::Update() {
    Time = PutativeTimesQueue.TimeQueue[0];
    int executed_action = PutativeTimesQueue.ActionQueue[0];
    ExecuteAction(executed_action);
    double new_action_rate;
    double new_putative_time;
    for (auto action : DependencyGraph[executed_action]) {
        new_action_rate = FindActionRate(action);
        if (action == executed_action) {
            if (new_action_rate == 0) {
                new_putative_time = numeric_limits<double>::infinity();
            }
            else {
                new_putative_time = -log(UnifDist(RandEngine))/new_action_rate+Time;
            }
            ActionRates[action] = new_action_rate;
        }
        else {
            if (new_action_rate == 0) {
                if (PreInfiniteTime[action] == 0) {
                    PreInfiniteTime[action] = PutativeTimesQueue.TimeQueue[PutativeTimesQueue.NodeArray[action]]-Time;
                }
                new_putative_time = numeric_limits<double>::infinity();
            }
            else {
                if (ActionRates[action] == 0) {
                    new_putative_time = -log(UnifDist(RandEngine))/new_action_rate+Time;
                    PreInfiniteTime[action] = 0;
                }
                else if (PreInfiniteTime[action] != 0){
                    new_putative_time = ActionRates[action]/new_action_rate*PreInfiniteTime[action]+Time;
                    PreInfiniteTime[action] = 0;
                }
                else{
                    new_putative_time = ActionRates[action]/new_action_rate*(PutativeTimesQueue.TimeQueue[PutativeTimesQueue.NodeArray[action]]-Time)+Time;
                }
                ActionRates[action] = new_action_rate;
            }
        }
        PutativeTimesQueue.Update(PutativeTimesQueue.NodeArray[action], new_putative_time);
    }
}

// Check whether a given population is extinct
bool Lattice::IsExtinct(int population) {
    bool extinct = true;
    for (int genotype = 0; genotype < LProp.GenBound - 1; genotype++) {
        if (GenTot[population][genotype] > 0) {
            extinct = false;
            break;
        }
    }
    return extinct;
}