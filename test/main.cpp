#include <ilcplex/ilocplex.h>

#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>

#include "graph.h"
#include "readers.h"
#include "smp.h"
#include "type.h"

ILOSTLBEGIN

int main(int argc, char** argv) {
    // data/part_graph_1.txt 4 1 0 0 1200 1156 0.36 1156 0.36
    // data/random_graph/plan_random/group_1/dataset1_1_1_1/animal_1.txt 4 1 0 0
    // 3 4 10 20 10 2 0 1200 1156 0.36 1156 0.36
    SmpForm formulation = NS;
    string filename = "";

    // Check the command line arguments
    // [filename] [formulation] [callbackOption] [relax] [time_limit] [max_cuts]
    // [epsilon] parsing formulation arg
    filename = argv[1];
    cout << filename << endl;

    int formulationOption = stoi(argv[2]);
    switch (formulationOption) {
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
        cout << " use both seperation " << endl;  // 1
    else
        cout << " use their own " << endl;  // 0 own

    // Local Branch parameter
    int LB_MaxRestart = atoi(argv[6]);
    int LB_MaxIter = atoi(argv[7]);
    int Rmin = atoi(argv[8]);
    int Rmax = atoi(argv[9]);
    int BCSolNum = atoi(argv[10]);
    int BCTime = atoi(argv[11]);
    int MIPDisplayLevel = atoi(argv[12]);

    // set timelimit and cuts number and constraint add tolerance index
    int time_limit = atoi(argv[13]);
    int max_cuts_lazy = atoi(argv[14]);
    double epsilon_lazy = atof(argv[15]);
    int max_cuts_user = atoi(argv[16]);
    double epsilon_user = atof(argv[17]);
    int UseLocalBranch = atoi(argv[18]);
    int LB_CP_Option = atoi(argv[19]);
    int lazy_sep_opt = atoi(argv[20]);

	cout << "max_cuts_lazy: " << max_cuts_lazy << endl;
	cout << "epsilon_lazy: " << epsilon_lazy << endl;
	cout << "max_cuts_user: " << max_cuts_user << endl;
	cout << "epsilon_user: " << epsilon_user << endl;

	if (UseLocalBranch) {
		cout << "LB_MaxRestart: " << LB_MaxRestart << endl;
		cout << "LB_MaxIter: " << LB_MaxIter << endl;
		cout << "Rmin: " << Rmin << endl;
		cout << "Rmax: " << Rmax << endl;
		cout << "BCSolNum: " << BCSolNum << endl;
		cout << "BCTime: " << BCTime << endl;
		cout << "MIP Display Level: " << MIPDisplayLevel << endl;
		cout << "Use Local Branch: " << UseLocalBranch << endl;
		cout << "Use Cutpool as: ";
		if (!LB_CP_Option)cout << " constraint " << endl;
		else cout << " cut " << endl;
	}

	cout << "seperate lazy constraint as: ";
	if (!lazy_sep_opt)cout << " one to multi " << endl;
	else cout << " multi to multi " << endl;

    // Read graph into G:
    Reader myReader;
    std::shared_ptr<Graph> G = std::make_shared<Graph>();
    myReader.read_graph(filename, G);

    // Build Model
    map<NODE, double> cost;
    cost = G->nodes_value();
    IloEnv SMPenv, LBenv;
    if (UseLocalBranch) {
        cout << "Begin to execute LBSolver() ..." << endl;
        LBSolver lb_solver =
            LBSolver(LBenv, G, formulation, callbackOption, relax, ns_sep_opt,
                     LB_MaxRestart, LB_MaxIter, Rmin, Rmax, BCSolNum, BCTime,
                     epsilon_lazy, epsilon_user, max_cuts_lazy, max_cuts_user,
                     filename, MIPDisplayLevel, LB_CP_Option, lazy_sep_opt);
        lb_solver.update_LB_problem();

        // Solve in cplex
        lb_solver.LocalBranchSearch();

        LBenv.end();
    } else {
        cout << "Begin to execute SmpSolver() ..." << endl;
        SmpSolver smp_solver =
            SmpSolver(SMPenv, G, formulation, epsilon_lazy, epsilon_user,
                      time_limit, max_cuts_lazy, max_cuts_user, callbackOption,
                      relax, ns_sep_opt, filename, LB_CP_Option, lazy_sep_opt);
        smp_solver.update_problem(cost, formulation);

        // Solve in cplex
        smp_solver.solve();

        SMPenv.end();
    }

    return 0;
}