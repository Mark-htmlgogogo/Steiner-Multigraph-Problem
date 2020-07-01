#pragma once

#include <ilcplex/ilocplex.h>

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "graph.h"
#include "separation.h"
#include "type.h"

using namespace std;

class LBSolver {
   public:
    LBSolver(IloEnv env, std::shared_ptr<Graph> g_ptr, SmpForm formulation,
             int callbackOption, bool relax, bool ns_sep_out, int LB_CallTime,
             int LB_InterationTime, int R, int BCSolNum, int BCTime,
             double epsilon_lazy, double epsilon_user, int max_cuts_lazy,
             int max_cuts_user, IloNumVarArray x_vararray,
             IloNumVarArray x_vararray_primal,
             map<NODE, IloNumVar> primal_node_vars,
             map<pair<NODE, INDEX>, IloNumVar> partition_node_vars,
             map<INDEX, NODE> ns_root,
             map<pair<NODE, INDEX>, int> x_varindex_ns,
             map<NODE, int> x_varindex_ns_primal);

    // Operation function
    void Floyd(map<NODE, vector<int>>& subGnodesIdx,
               map<int, NODE>& rev_subGnodesIdx, vector<vector<int>>& distance,
               vector<vector<int>>& path, int& idx, int k);
    void GenerateInitialSolution(int k, map<NODE, bool>& xPartSol);
    void LocalBranchSearch();
    void LocalBranch();

    // Varaible interface
    const vector<IloExpr> cutPoolLhs() const { return _cutPoolLhs; }
    const vector<IloExpr> cutPoolRhs() const { return _cutPoolRhs; }
    const vector<double> violation() const { return _violation; }

   private:
    const int MAXN = 1000;
    const int INF = 0x3f3f3f3f;

    // CPLEX varaible
    IloModel model;
    IloCplex cplex;
    IloObjective objective;

    // LB varaible
    int LB_CallTime;        // K-th time call Local Branch
    int LB_InterationTime;  // Local Branch search time
    int R;                  // Maximum replaceable neighbor
    int BCSolNum;           // Cplex solving B&C number
    int BCTime;             // Cplex solving B&C time

    vector<IloExpr> _cutPoolLhs, _cutPoolRhs;  // Cut pool expression.
    vector<double> _violation;                 // Cut pool violation

    // General Varaible
    std::shared_ptr<Graph> G;
    SmpForm formulation;
    int callbackOption;
    bool relax;
    bool ns_sep_opt;
    int max_cuts_lazy;
    int max_cuts_user;
    double tol_lazy;
    double tol_user;

    // NS varaible
    IloNumVarArray x_vararray;
    IloNumVarArray x_vararray_primal;
    map<NODE, IloNumVar> primal_node_vars;                  // x_i
    map<pair<NODE, INDEX>, IloNumVar> partition_node_vars;  // x_i^k
    map<INDEX, NODE> ns_root;
    map<pair<NODE, INDEX>, int> x_varindex_ns;
    map<NODE, int> x_varindex_ns_primal;
};