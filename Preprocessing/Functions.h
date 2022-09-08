# pragma once
# include "Properties.h"
# include "Lattice.h"
# include <iostream>
# include <vector>
# include <vector>

using namespace std;

// Function to check whether adaptation numbers were overcome in at least one population.
vector<bool> IsAdapted(int pop_num, vector<vector<int>> adap_gen, vector<int> adap_num);

// Save GenSpaceTot to a .csv file
void SaveGenSpaceTot(Lattice &latt, string file_name);

// Write preamble to a .csv file
void WritePreamble(Lattice &latt, ofstream &output_file);

// Find the adaption rate and snapshots of a given population with given l_prop, and p_prop.
// Adap_num adaptations are computed for (T, D, P - same for each population), unless the lattice time goes over max_iterate.
void AdapImage(LatProp l_prop, vector<PopProp> p_prop, int adap_number, string file_specifier, long int max_iterate);

// Produces snapshots of GenSpaceTot at multiples of sim_interval until stop_time, save under file_specifier.
void SnapshotTot(LatProp l_prop, vector<PopProp> p_prop, double sim_interval, double stop_time, string file_specifier);

// Compute the stable wild-type population with a given absolute precision.
// Output to .csv file.
void FindStablePop(LatProp l_prop, vector<PopProp> p_prop, int population, double precision, string file_specifier);

// Produce plotable output for animation with sim_length frames, and of stop_time lattice time duration.
void AnimTot(LatProp l_prop, vector<PopProp> p_prop, int sim_length, double stop_time, string file_specifier, int population);