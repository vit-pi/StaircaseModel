#pragma once
#include "Properties.h"
#include <vector>

using namespace std;

// INDEX PRIORITY QUEUE HANDLING FUNCTIONS
class IndexPriorityQueue
{
public:
	vector<double> TimeQueue; //saves times
	vector<int> ActionQueue; //saves actions
	vector<int> NodeArray; //for each action saves the node index
	void Build(vector<double> delta_time);
	void Update(int node, double new_delta_time);
private:
	int QueueSize = 0;
	void Heapify(int node);
	void Swap(int node1, int node2);
	void UpdateAux(int node);
};

// NEWTON-RAPHSON
class NewtonSolve // Solve for stable population via Newton-Raphson and Tridiagonal matrix algorithm
{
public:
	NewtonSolve(double Iprecision, int IK, double Ir, double Is_r, double Id, double Is_d, double Inu_r, double Is_nu_r, double Inu_l, double Is_nu_l, int ID, int IGridBound);
	vector<double> Solve();
private:
	double precision;
	double r;
	double s_r;
	double d;
	double s_d;
	double nu_r;
	double s_nu_r;
	double nu_l;
	double s_nu_l;
	double K;
	int D;
	int GridBound;
	vector<double> step;  // Newton-Raphson step
	vector<double> J;  // Diagonal elements of the Jacobian matrix
	vector<double> MinusF; // Minus the functions whose root we seek
	vector<double> N; // Current estimate on N
	vector<double> nu_right;  // Vector of lower diagonal Jacobian entries
	vector<double> nu_left; // Vector of upper diagonal Jacobian entries
	double Magnitude(vector<double> vec);
};

// LATTICE CLASS
class Lattice
{
public:
	// Properties and Constant Variables
	LatProp LProp;
	vector<PopProp> PProp;

	// Dynamical Variables
	vector<vector<vector<int>>> GenSpaceTot;
	vector<vector<int>> SpaceTot;
	vector<vector<int>> GenTot;
	vector<int> AbsSpaceTot;
	vector<int> AbsSpaceTotUnderStair;
	vector<int> AbsGenTot;
	vector<vector<int>> StablePop;
	int AbsTot;
	double Time;

	// Methods
	Lattice(LatProp l_prop, vector<PopProp> p_prop);
	void Update();
	bool IsExtinct(int population);
	void StablePopulation(double precision, int population);
	
private:
	vector<vector<int>> Actions;
	vector<vector<int>> DependencyGraph;
	vector<double> ActionRates;
	vector<double> PreInfiniteTime;
	vector<double> PreInfinitePutativeTime;
	IndexPriorityQueue PutativeTimesQueue;
	double FindActionRate(int action);
	void ExecuteAction(int action);
};