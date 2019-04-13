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

int main(int argc, char** argv)
{

	SmpForm formulation = NS;
	string filename = "";

	// Check the command line arguments
	// [filename] [formulation] [callbackOption] [relax] [time_limit] [max_cuts] [epsilon]
	// parsing formulation arg
	filename = argv[1];
	int formulationOption = stoi(argv[2]);
	switch (formulationOption)
	{
	case 1:
		cout << "formulation : SCF" << endl;
		formulation = SCF;
		break;
	case 2:
		cout << "formulation : MCF" << endl;
		formulation = MCF;
		break;
	case 3:
		cout << "formulation : STEINER" << endl;
		formulation = STEINER;
		break;
	case 4:
		cout << "formulation : NS" << endl;
		formulation = NS;
		break;
	default:
		break;
	}

	// parsing callback option arg
	int callbackOption = atoi(argv[3]);
	cout << "formulation option " << formulationOption << endl;
	cout << "callback option " << callbackOption << endl;

	bool relax = atoi(argv[4]);
	if (relax)
		cout << "relax applied" << endl;
	else
		cout << "no relax" << endl;

	bool ns_sep_opt = atoi(argv[5]);
	if (ns_sep_opt)
		cout << " use both seperation " << endl;
	else
		cout << " use their own " << endl;

	//set timelimit and cuts number and constraint add tolerance index
	int time_limit = atoi(argv[6]);
	int max_cuts = atoi(argv[7]);
	double epsilon = atof(argv[8]);

	// Read graph into G:
	Reader myReader;
	std::shared_ptr<Graph>G = std::make_shared<Graph>();
	myReader.read_graph(filename, G);

	// Build Model
	map<NODE, double> cost;
	cost = G->nodes_value();
	IloEnv env;
	cout << "Begin to execute SmpSolver() ..." << endl;
	SmpSolver smp_solver = SmpSolver(env, G, formulation, epsilon, time_limit, max_cuts, callbackOption, relax, ns_sep_opt, filename);
	smp_solver.update_problem(cost);

	// Solve in cplex
	smp_solver.solve();

	env.end();

	return 0;
}