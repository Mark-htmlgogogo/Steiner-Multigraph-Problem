#pragma once
#include "type.h"
#include "graph.h"
#include "separation.h"
#include <ilcplex/ilocplex.h>
#include <map>
#include <memory>
#include <map>
#include <string>

using namespace std;

typedef map<NODE, IloNumVar> NODE_VARS;
typedef map<NODE, vector<NODE>> ADJ_LIST;

class SmpSolver
{
public:
	SmpSolver(IloEnv env,
	          std::shared_ptr<Graph>g_ptr,
	          SmpForm formulation,
	          double epsilon_lazy,
			  double epsilon_user,
	          int time_limit,
	          int max_cuts_lazy,
			  int max_cuts_user,
	          int callbackOption,
	          bool relax,
	          bool ns_sep_opt,
	          string filename);

	void update_problem(const map<NODE, double> &obj_coeff);
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
	map<NODE, IloNumVar> primal_node_vars;					  // x_i
	map<pair<NODE, INDEX>, IloNumVar> partition_node_vars;	  // x_i^k
	IloNumVarArray x_vararray;

	/* SCF */
	map<INDEX, IloNumVar> source_node_vars;					  // x_s_k
	map<pair<NODE_PAIR, INDEX>, IloNumVar> partition_flow_vars; // y_ij^k

	/* MCF */
	map<pair<NODE_PAIR, NODE_PAIR>, IloNumVar> multi_flow_vars; //y_ij_km

	/*	MCF teimanial */
	map<pair<NODE, NODE_PAIR>, IloNumVar>path_flow_vars;

	/* STRINER */
	map<pair<NODE_PAIR, INDEX>, IloNumVar> edge_vars;			  //y_ij_k
	map<pair<NODE_PAIR, INDEX>, int> x_varindex_steiner;
	map<INDEX, NODE> Steiner_root;

	/* NS */
	map<INDEX, NODE> ns_root;
	map<pair<NODE, INDEX>, int> x_varindex_ns;

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
	std::shared_ptr<Graph>G;
	double elapsed_time;
	double elapsed_ticks;
	string filename;
};
