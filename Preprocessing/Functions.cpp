#pragma once
#include "Functions.h"
#include "Lattice.h"
#include "Properties.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <string>

using namespace std;

// Function to check whether adaptation numbers were overcome in at least one population.
// NotAdapted (for T, for D, for P, for at least one of them)
vector<bool> NotAdapted(int pop_num, vector<vector<int>> adap_gen, vector<int> adap_num){
	vector<bool> not_adapted_array(4,false);
	bool exists_not_adapted;
	bool not_adapted;
	for (int time_type=0; time_type<3; time_type++) {
		not_adapted = true;
		for (int population=0; population<pop_num; population++) {
			if (adap_gen[time_type][population]>=adap_num[time_type]){
				not_adapted = false;
			}
		}
		not_adapted_array[time_type] = not_adapted;
	}
	if (not_adapted_array[0] || not_adapted_array[1]) {
		not_adapted_array[3] = true;
	}
	return not_adapted_array;
}

// Save GenSpaceTot to a .csv file
void SaveGenSpaceTot(Lattice &latt, string file_name){
	ofstream OutputFile(file_name);
	OutputFile << "GenSpaceTot\n";
	OutputFile << "Time," << latt.Time << "\n";
	for (int population=0; population < latt.LProp.PopulationNumber; population++) {
		OutputFile << "Population," << population << "\n";
		for (int genotype = 0; genotype < latt.LProp.GenBound; genotype++) {
			for (int position = 0; position < latt.LProp.GridBound; position++) {
				OutputFile << latt.GenSpaceTot[population][genotype][position] << ", ";
			}
			OutputFile << "\n";
		}
	}
	OutputFile.close();
}

// Write preamble to a .csv file
void WritePreamble(Lattice &latt, ofstream &output_file){
	// Summary of lattice properties
	output_file << "LatticeProperties\n";
	output_file << "PopulationNumber," << latt.LProp.PopulationNumber << "\n";
	output_file << "GridBound," << latt.LProp.GridBound << "\n";
	output_file << "GenBound," << latt.LProp.GenBound << "\n";
	output_file << "CarryingCapacity," << latt.LProp.CarryingCapacity << "\n";
	output_file << "D," << latt.LProp.D << "\n";
	output_file << "\n";
	// Summary of initial conditions
	for (int population=0; population<latt.LProp.PopulationNumber; population++) {
		output_file << "Population," << to_string(population) << "\n";
		if(latt.PProp[population].AutomaticInitiation){
		output_file << "StablePop,";
		for (int position = 0; position<latt.LProp.GridBound; position++){
			output_file << latt.StablePop[population][position] << ",";
		}
		output_file << "\n";
		}
		else{
			output_file << "InitCellsNum," << latt.PProp[population].InitCellsNum << "\n";
			output_file << "InitCellsGen," << latt.PProp[population].InitCellsGen << "\n";
			output_file << "InitCellsPos," << latt.PProp[population].InitCellsPos << "\n";
		}
		output_file << "BirthRate," << latt.PProp[population].BirthRate << "\n";
		output_file << "StressBirthRate," << latt.PProp[population].StressBirthRate << "\n";
		output_file << "ResistCost,"  << latt.PProp[population].ResistCost << "\n";
		output_file << "DeathRate," << latt.PProp[population].DeathRate << "\n";
		output_file << "StressDeathRate," << latt.PProp[population].StressDeathRate << "\n";
		output_file << "MuteUp," << latt.PProp[population].MuteUp << "\n";
		output_file << "MuteDown," << latt.PProp[population].MuteDown << "\n";
		output_file << "StressMuteUp," << latt.PProp[population].StressMuteUp << "\n";
		output_file << "StressMuteDown," << latt.PProp[population].StressMuteDown << "\n";
		output_file << "Move," << latt.PProp[population].Move << "\n";
		output_file << "Chemotax," << latt.PProp[population].Chemotax << "\n";
		output_file << "StressDepChemotax," << latt.PProp[population].StressDepChemotax << "\n";
		output_file << "Chemokin," << latt.PProp[population].Chemokin << "\n";
		output_file << "SwitchUp," << latt.PProp[population].SwitchUp << "\n";
		output_file << "SwitchDown," << latt.PProp[population].SwitchDown << "\n";
		output_file << "ConsiderSwarm," << latt.PProp[population].ConsiderSwarm << "\n";
		output_file << "SwarmDens," << latt.PProp[population].SwarmDens << "\n";
		output_file << "SwarmMove," << latt.PProp[population].SwarmMove << "\n";
		output_file << "ConsiderHGT," << latt.PProp[population].ConsiderHGT << "\n";
		output_file << "HGTRate," << latt.PProp[population].HGTRate << "\n"; // added 5 new parameters
		output_file << "CompetBelowStairs," << latt.PProp[population].CompetBelowStairs << "\n";
		output_file << "ConsiderDensSwitch," << latt.PProp[population].ConsiderDensSwitch << "\n";
		output_file << "DensSwitchBias," << latt.PProp[population].DensSwitchBias << "\n";
		output_file << "DensSwitchTot," << latt.PProp[population].DensSwitchTot << "\n";
		output_file << "DensSwitchDens," << latt.PProp[population].DensSwitchDens << "\n";
		output_file << "\n";
	}
}

// Find the adaption rate and snapshots of a given population with given l_prop, and p_prop.
// Adap_num adaptations are computed, unless the lattice time goes over max_iterate
// Returns times: 0 = mutant appears in O, 1 = mutants outcompete wildtype in O!!!
void AdapImage(LatProp l_prop, vector<PopProp> p_prop, int adap_number,
	string file_specifier, long int max_iterate) {
	// Prepare variables
	vector<vector<vector<double>>> WaitingTime(l_prop.PopulationNumber, vector<vector<double>>(l_prop.GridBound, vector<double>(3,0)));
	long int iterate = 0;
	// Prepare adap_gen for T=0, D=1, P=2
	vector<int> adap_num(3,adap_number);
	adap_num[2] = l_prop.GridBound;
	vector<vector<int>> adap_gen(3, vector<int>(l_prop.PopulationNumber, 0));
	string file_name;
	// Do the simulation
	// Create the lattice and save the initial conditions into a separate file
	Lattice latt(l_prop, p_prop);
	file_name = "OutAdap_" + file_specifier + ".csv";
	ofstream OutputFileA(file_name);
	WritePreamble(latt, OutputFileA);
	OutputFileA.close();

	// Do the lattice simulation
	for (int population=0; population<l_prop.PopulationNumber; population++) {
		adap_gen[0][population] = latt.PProp[population].InitCellsGen + 1; // controlling the first mutant in the sink
		adap_gen[1][population] = latt.PProp[population].InitCellsGen + 1; // controlling an established pop. in the sink
		adap_gen[2][population] = latt.PProp[population].InitCellsPos; // controlling a phenotypic resistance (presence of swarm in adap_gen compartment) when ConsiderSwarm
	}
	adap_num[0] += 1;
	adap_num[1] += 1;
	vector<bool> not_adapted = NotAdapted(l_prop.PopulationNumber,adap_gen, adap_num);
	while (not_adapted[3]) {
		for (int population=0; population<l_prop.PopulationNumber; population++) {
			if ((not_adapted[0]) && (latt.GenSpaceTot[population][adap_gen[0][population]][adap_gen[0][population]+latt.LProp.D-1] > 0)) {
				WaitingTime[population][adap_gen[0][population] - 1 - latt.PProp[population].InitCellsGen][0] = latt.Time;
				file_name= "OutAdap_" + file_specifier + "_" + to_string(population) + "_"+to_string(adap_gen[0][population]) + "T.csv";
				SaveGenSpaceTot(latt, file_name);
				adap_gen[0][population]++;
				//if (adap_gen[0][population]-adap_gen[1][population]>1){
				//	adap_gen[1][population]++;
				//}	
			}
			if ((not_adapted[1]) &&(latt.GenSpaceTot[population][adap_gen[1][population]][adap_gen[1][population]+latt.LProp.D-1] > latt.GenSpaceTot[population][adap_gen[1][population]-1][adap_gen[1][population]+latt.LProp.D-1])) {
				WaitingTime[population][adap_gen[1][population] - 1 - latt.PProp[population].InitCellsGen][1] = latt.Time;
				file_name= "OutAdap_" + file_specifier + "_" + to_string(population) + "_"+to_string(adap_gen[1][population]) + "D.csv";
				SaveGenSpaceTot(latt, file_name);
				adap_gen[1][population]++;
			}
			if (latt.PProp[population].ConsiderSwarm){
				if ((not_adapted[2]) && (latt.SpaceTot[population][adap_gen[2][population]]>latt.PProp[population].SwarmDens)){
					WaitingTime[population][adap_gen[2][population] - latt.PProp[population].InitCellsPos][2] = latt.Time;
					file_name= "OutAdap_" + file_specifier + "_" + to_string(population) + "_"+to_string(adap_gen[2][population]) + "PA.csv";
					SaveGenSpaceTot(latt, file_name);
					adap_gen[2][population]++;
				}
			}
			if (latt.PProp[population].ConsiderDensSwitch){
				if ((not_adapted[2]) && (latt.AbsSpaceTot[adap_gen[2][population]]>latt.PProp[population].DensSwitchDens)){
					WaitingTime[population][adap_gen[2][population] - latt.PProp[population].InitCellsPos][2] = latt.Time;
					file_name= "OutAdap_" + file_specifier + "_" + to_string(population) + "_"+to_string(adap_gen[2][population]) + "PA.csv";
					SaveGenSpaceTot(latt, file_name);
					adap_gen[2][population]++;
				}
			}
		}
		latt.Update();
		not_adapted = NotAdapted(l_prop.PopulationNumber,adap_gen, adap_num);
		iterate++;
		// If one simulation takes more iterates than max_iterate, terminate
		if ((iterate > max_iterate) || (latt.AbsTot == 0)) {
			not_adapted[3]=false;
		}
	}

	// Write final report
	file_name = "OutAdap_" + file_specifier + "_Times.csv";
	ofstream OutputFileB(file_name);
	//  Write the adaptation times into the file
	for (int population=0; population<l_prop.PopulationNumber; population++) {
		OutputFileB << "Population," << to_string(population)<<"\n";
		// If not ConsiderSwarm
		if (!p_prop[population].ConsiderSwarm) {
			for (int type=0; type<2; type++){
				OutputFileB << "TimeType(0=T,1=D)," << to_string(type) << "\n";
				for (int adap_gen=0; adap_gen<adap_number; adap_gen++){
					OutputFileB << WaitingTime[population][adap_gen][type] << ",";
				}
				OutputFileB << "\n";
			}
		}
		// If ConsiderSwarm
		else {
			for (int type=0; type<3; type++){
				OutputFileB << "TimeType(0=T,1=D,2=PA)," << to_string(type) << "\n";
				for (int adap_gen=0; adap_gen<l_prop.GridBound; adap_gen++){
					OutputFileB << WaitingTime[population][adap_gen][type] << ",";
				}
				OutputFileB << "\n";
			}
		}
	}
	OutputFileB.close();
}

// Produce plotable output for animation with sim_length frames, and of stop_time lattice time duration.
void SnapshotTot(LatProp l_prop, vector<PopProp> p_prop, double sim_interval, double stop_time, string file_specifier) {
	// Generate lattice object and initialize variables
	Lattice latt(l_prop, p_prop);
	int sim_index = 0;
	
	// Create output file and initialize it
	string file_name = "OutSnapshot_" + file_specifier + ".csv";
	ofstream OutputFile(file_name);
	OutputFile << "Snapshot properties\n";
	OutputFile << "SimInterval," << sim_interval << "\n";
	OutputFile << "StopTime," << stop_time << "\n";
	OutputFile << "\n";
	WritePreamble(latt, OutputFile);
	OutputFile.close();

	
	// Run the lattice model and write the GenSpaceTot to the file
	while (latt.Time < stop_time)
	{
		// Write current GenSpaceTot array to file
		if (latt.Time >= (double)sim_index * sim_interval) {
			// Write time on console
			cout << latt.Time << "\n";
			// Write to file
			SaveGenSpaceTot(latt, "OutSnapshot_" + file_specifier + "_"+to_string(sim_index)+".csv");	
			sim_index++;
		}
		// Update lattice
		latt.Update();
	}
}

// Compute the stable wild-type population with a given absolute precision.
// Output to .csv file
void FindStablePop(LatProp l_prop, vector<PopProp> p_prop, int population, double precision, string file_specifier){
	Lattice latt(l_prop, p_prop);
	latt.StablePopulation(precision, population);
	string file_name = "OutStablePop_" + file_specifier + ".csv";
	ofstream OutputFile(file_name);
	for (int position = 0; position<latt.LProp.GridBound; position++){
		OutputFile << latt.StablePop[population][position] << ",";
	}
	OutputFile << "\n";
	WritePreamble(latt, OutputFile);
}

// Produce plotable output for animation with sim_length frames, and of stop_time lattice time duration.
void AnimTot(LatProp l_prop, vector<PopProp> p_prop, int sim_length, double stop_time, string file_specifier, int population) {
	// Generate lattice object
	Lattice latt(l_prop, p_prop);

	// Decide on plotting parameters
	double sim_interval = stop_time / sim_length;
	int sim_index = 0;

	// Create output file and initialize it
	string FileName = "OutAnimTot_" + file_specifier + ".py";
	ofstream OutputFile(FileName);
	OutputFile << "import numpy as np\n\n";
	OutputFile << "SimLength = " << sim_length << "\n";
	OutputFile << "ConsiderFirstArrival = False\n";
	OutputFile << "StopTime = " << stop_time << "\n";
	OutputFile << "Times = np.zeros(" << sim_length << ", dtype=float)\n";
	OutputFile << "GenSpaceTot = np.zeros((" << sim_length << ", " << l_prop.PopulationNumber << ", " << l_prop.GenBound << ", " << l_prop.GridBound << "), dtype=float)\n";
	// Run the lattice model and write the GenSpaceTot to the file
	while (latt.Time < stop_time)
	{
		// Write current GenSpaceTot array to file
		if (latt.Time >= (double)sim_index * sim_interval) {
			// Write time on console
			cout << latt.Time << "\n";
			// Write to file
			OutputFile << "Times[" << sim_index << "] = " << latt.Time << "\n";
			for (population=0; population < l_prop.PopulationNumber; population++) {
				OutputFile << "GenSpaceTot[" << sim_index << "][" << population << "] = np.asarray([";
				for (int genotype = 0; genotype < l_prop.GenBound; genotype++) {
					OutputFile << "[";
					for (int position = 0; position < l_prop.GridBound; position++) {
						OutputFile << latt.GenSpaceTot[population][genotype][position] << ", ";
					}
					OutputFile << "],\n";
				}
				OutputFile << "])\n";
			}	
			sim_index++;
		}
		// Update lattice
		latt.Update();
	}
}
