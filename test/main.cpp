#include <ilcplex/ilocplex.h>

#include "type.h"
#include "smp.h"
#include "graph.h"
#include "readers.h"

#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>

ILOSTLBEGIN

int main()
{
	int time_limit = 1200;
	int max_cuts = 1200;
	double epsilon = 0;
	SmpForm formulation = NS;
	string filename = "part_graph_binary.txt";

	// Read graph into G:
	Reader myReader;
	std::shared_ptr<Graph>G = std::make_shared<Graph>();
	myReader.read_graph(filename, G);

	// Build Model
	map<NODE, double> cost;
	cost = G->nodes_value();
	IloEnv env;
	cout << "Begin to execute SmpSolver() ..." << endl;
	SmpSolver smp_solver = SmpSolver(env, G, formulation, epsilon, time_limit, max_cuts);
	smp_solver.update_problem(cost);

	// Solve in cplex
	smp_solver.solve();

	env.end();


	cout << endl << endl;
	return 0;
}