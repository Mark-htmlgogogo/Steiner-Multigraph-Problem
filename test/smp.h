#pragma once
#include <ilcplex/ilocplex.h>

#include <map>
#include <memory>
#include <string>

#include "graph.h"
#include "separation.h"
#include "type.h"

using namespace std;

typedef map<NODE, IloNumVar> NODE_VARS;
typedef map<NODE, vector<NODE>> ADJ_LIST;

class SmpSolver {
   public:
    SmpSolver(IloEnv env, std::shared_ptr<Graph> g_ptr, SmpForm formulation,
              double epsilon_lazy, double epsilon_user, int time_limit,
              int max_cuts_lazy, int max_cuts_user, int callbackOption,
              bool relax, bool ns_sep_opt, string filename, int LB_CP_Option,
              int lazy_sep_opt, int MIPDisplayLevel);

    void update_problem(const map<NODE, double>& obj_coeff,
                        SmpForm formulation);
    void print_to_file();

    double elpased_time() { return elapsed_time; }
    double elpased_ticks() { return elapsed_ticks; }

    void solve();

    void clear();

   private:
    void build_problem_scf();
    void build_problem_mcf();
    void build_problem_mcf_terminal();
    void build_problem_steiner();
    void build_problem_ns();

    void printInfo(IloNumVar var);

    /* All used variable */
    map<NODE, IloNumVar> primal_node_vars;                  // x_i
    map<pair<NODE, INDEX>, IloNumVar> partition_node_vars;  // x_i^k
    IloNumVarArray x_vararray;
    IloNumVarArray x_vararray_primal;

    /* SCF */
    map<INDEX, IloNumVar> source_node_vars;                      // x_s_k
    map<pair<NODE_PAIR, INDEX>, IloNumVar> partition_flow_vars;  // y_ij^k

    /* MCF */
    map<pair<NODE_PAIR, NODE_PAIR>, IloNumVar> multi_flow_vars;  // y_ij_km

    /*	MCF teimanial */
    map<pair<NODE, NODE_PAIR>, IloNumVar> path_flow_vars;

    /* STRINER */
    map<pair<NODE_PAIR, INDEX>, IloNumVar> edge_vars;  // y_ij_k
    map<pair<NODE_PAIR, INDEX>, int> x_varindex_steiner;
    map<INDEX, NODE> Steiner_root;

    /* NS */
    map<INDEX, NODE> ns_root;
    map<pair<NODE, INDEX>, int> x_varindex_ns;
    map<NODE, int> x_varindex_ns_primal;

    /* Local Branch */
    int Obj;                     // Obj value corresponding to final solution
    map<NODE, bool> xPrimalSol;  // Final primal solution
    map<INDEX, map<NODE, bool>> xPartSol;  // Final partiton solution

    IloModel model;
    IloCplex cplex;
    IloObjective objective;

    int time_limit;
    int ncuts;
    int max_cuts_lazy;
    int max_cuts_user;
    double tol_lazy;
    double tol_user;
    SmpForm formulation;
    int callbackOption;
    bool relax;
    bool ns_sep_opt;
    std::shared_ptr<Graph> G;
    double elapsed_time;
    double elapsed_ticks;
    string filename;
    int LB_CP_Option;
    int fianlsolveflag;
    int lazy_sep_opt;
	int MIPDisplayLevel;
};

class LBSolver {
   public:
    LBSolver(IloEnv env, std::shared_ptr<Graph> g_ptr, SmpForm formulation,
             int callbackOption, bool relax, bool ns_sep_out,
             int LB_MaxRestarts, int LB_MaxIter, int Rmin, int Rmax,
             int BCSolNum, int BCTime, double epsilon_lazy, double epsilon_user,
             int max_cuts_lazy, int max_cuts_user, string filename,
             int MIPDisplayLevel, int LB_CP_Option, int lazy_sep_opt);
    void Floyd(map<NODE, int>& subGnodesIdx, map<int, NODE>& rev_subGnodesIdx,
               vector<vector<int>>& distance, vector<vector<int>>& path,
               int& idx, int k);
    void GenerateInitialSolution(int k);
    void update_LB_problem();
    void update_LB_problem_final();
    void build_problem_ns_simplifer();
    void build_problem_ns_final();
    void LocalBranchSearch();
    void LocalBranch(int& ObjValue);
    void FinalSolve();
    void CheckSolution();
    void print_to_file();

    // interference of varaible
    const map<NODE, map<NODE, bool>>& FxPartSol() const { return xPartSol; }
    const map<NODE, bool>& FxPrimalSol() const { return xPrimalSol; }
    const int& Fobj() const { return Final_Obj; }

   private:
    const int MAXN = 1000;
    const int INF = 0x3f3f3f3f;

    // CPLEX varaible
    IloModel LBmodel;
    IloCplex LBcplex;
    IloObjective LBobjective;

    IloModel FLBmodel;
    IloCplex FLBcplex;
    IloObjective FLBobjective;

    // LB varaible
    int LB_MaxRestarts;  // K-th time call Local Branch
    int LB_MaxIter;      // Local Branch search time
    int Rmin;            // Minimum replaceable neighbor
    int Rmax;            // Maximum replaceable neighbor
    int BCSolNum;        // Cplex solving B&C number
    int BCTime;          // Cplex solving B&C time

    double TOT_LB_TIME;
    double TOT_TIME;
    double FINAL_SOLVE_TIME;
    int LocalBranchTime;

    int Final_Obj;  // Obj value corresponding to final solution
    map<NODE, bool> Final_xPrimalSol;            // Final primal solution
    map<INDEX, map<NODE, bool>> Final_xPartSol;  // Final partiton solution
    map<NODE, bool> xPrimalSol;
    map<INDEX, map<NODE, bool>> xPartSol;
	NODE_SET HeuristicPool;

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
    string filename;
    int MIPDisplayLevel;
    int LB_CP_Option;
    int fianlsolveflag;
    int lazy_sep_opt;

    // NS varaible
    IloNumVarArray x_vararray;
    IloNumVarArray x_vararray_primal;
    map<NODE, IloNumVar> primal_node_vars;                  // x_i
    map<pair<NODE, INDEX>, IloNumVar> partition_node_vars;  // x_i^k
    map<INDEX, NODE> ns_root;
    map<pair<NODE, INDEX>, int> x_varindex_ns;
    map<NODE, int> x_varindex_ns_primal;
};