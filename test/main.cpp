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
    // 1200 1156 0.36 1156 0.36
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

    // set timelimit and cuts number and constraint add tolerance index
    int time_limit = atoi(argv[6]);
    int max_cuts_lazy = atoi(argv[7]);
    double epsilon_lazy = atof(argv[8]);
    int max_cuts_user = atoi(argv[9]);
    double epsilon_user = atof(argv[10]);
    int LB_MaxRestart = 5;
    int LB_MaxIter = 5;
    int Rmin = 10;
    int Rmax = 30;
    int BCSolNum = 10;
    int BCTime = 20;

    // Read graph into G:
    Reader myReader;
    std::shared_ptr<Graph> G = std::make_shared<Graph>();
    myReader.read_graph(filename, G);

    // Build Model
    map<NODE, double> cost;
    cost = G->nodes_value();
    IloEnv SMPenv, LBenv;
    cout << "Begin to execute LBSolver() ..." << endl;
    LBSolver lb_solver =
        LBSolver(LBenv, G, formulation, callbackOption, relax, ns_sep_opt,
                 LB_MaxRestart, LB_MaxIter, Rmin, Rmax, BCSolNum, BCTime,
                 epsilon_lazy, epsilon_user, max_cuts_lazy, max_cuts_user);
    lb_solver.update_LB_problem();

    // Solve in cplex
    lb_solver.LocalBranchSearch();

    cout << "Begin to execute SmpSolver() ..." << endl;
    SmpSolver smp_solver =
        SmpSolver(SMPenv, G, formulation, epsilon_lazy, epsilon_user,
                  time_limit, max_cuts_lazy, max_cuts_user, callbackOption,
                  relax, ns_sep_opt, filename);
    smp_solver.update_problem(cost, formulation);

    // Solve in cplex
    smp_solver.solve();

    SMPenv.end();
    LBenv.end();

    return 0;
}