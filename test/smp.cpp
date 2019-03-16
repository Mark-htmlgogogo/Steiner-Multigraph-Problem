#include "smp.h"
#include "callback.h"
#include "type.h"
#include "separation.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <numeric>

#define LOG    \
	if (false) \
	cerr
#define TOL 0.001

using namespace std;

/*******************************************************/
/************* Build the model *************************/
/*******************************************************/

/* Constructor */

SmpSolver::SmpSolver(
	IloEnv env,
	std::shared_ptr<Graph> g_ptr,
	SmpForm formulation_,
	double epsilon,
	int time_limit_,
	int max_cuts_)
{
	/* Initialize Cplex Sturctures */
	model = IloModel(env);
	objective = IloObjective();

	/* settings */
	formulation = formulation_;
	G = g_ptr;
	time_limit = time_limit_;
	max_cuts = max_cuts_;
	tol = epsilon;

	/* Add x_i variables: primal_node_vars */
	char var_name[255];
	set<NODE> T = G->t_total();
	int idx = 0;
	for (auto &node : G->nodes())
	{
		IloNumVar var;
		snprintf(var_name, 255, "x_%d", node);
		// set x_i = 1 if i belongs to T, others to {0,1}
		/*if (T.find(node) != T.end()) // Constraint(2)
			var = IloNumVar(env, 1, 1, IloNumVar::Bool, var_name);
		else
			var = IloNumVar(env, 0, 1, IloNumVar::Bool, var_name);*/
		var = IloNumVar(env, 0, 1, IloNumVar::Bool, var_name);
		primal_node_vars[node] = var;
		printInfo(var);
	}

	LOG << "Building problem...";
	cout << "Begin to build problem..." << endl;

	switch (formulation)
	{
	case SCF:
		build_problem_scf();
		break;
	case MCF:
		build_problem_mcf();
		break;
	case STEINER:
		build_problem_steiner();
		break;
	case NS:
		build_problem_ns();
		break;
	}

	LOG << "DONE." << endl;

	/* Cplex settings */
	cplex = IloCplex(model);				 // create a ILOG CPLEX algorithm and extract a model
	cplex.setParam(IloCplex::MIPDisplay, 3); // set display level
	if (formulation > 0)
		cplex.setParam(IloCplex::AdvInd, 1);	 // start value: 1
	cplex.setParam(IloCplex::EpGap, 1e-09);		 // set MIP gap tolerance
	cplex.setParam(IloCplex::Threads, 1);		 // set the number of parallel threads
	cplex.setParam(IloCplex::TreLim, 2048);		 // set the limit of tree memory in megabytes
	cplex.setParam(IloCplex::TiLim, time_limit); // set time limit in secs, default:1200
	cplex.setParam(IloCplex::EpInt, 1e-06);		 // set integrality tolerance

	/* Implement user cut or lazy constraint here*/
	switch (formulation)
	{
	case NONE:
	case SCF:
	case MCF:
		break;
	case STEINER:
	case NS:
		cplex.use(StrongComponentLazyCallback(env, G, edge_vars, x_vararray, x_varindex, tol, max_cuts, formulation, Steiner_root, primal_node_vars));
		cplex.use(SmpCutCallback(env,G,edge_vars,x_vararray,x_varindex,tol,max_cuts,formulation,Steiner_root,primal_node_vars));
	default:
		break;
	}
}

/* Update the objective function */
void SmpSolver::update_problem(const const map<NODE, double> &obj_coeff)
{
	/* Build objective function */
	IloEnv env = model.getEnv();
	IloExpr totalCost(env);
	double objcoeff = 0.0;
	IloNumVar var;

	for (auto &node : G->nodes())
	{
		var = primal_node_vars[node];
		objcoeff = obj_coeff.at(node);
		totalCost += var * objcoeff;
	}

	objective = IloObjective(env, totalCost, IloObjective::Minimize);
	model.add(objective);
}

void SmpSolver::solve()
{
	LOG << "--------------------- SOLVING SMP " << endl;

	//double start_time = cplex.getCplexTime();
	//double start_ticks = cplex.getDetTime();

	cout << "Number of constraints: " << cplex.getNrows() << endl;
	cout << "Number of variables(int):   " << cplex.getNintVars() << endl;
	cout << "Number of variables(binary):   " << cplex.getNbinVars() << endl;
	cout << "Number of variables(cols):   " << cplex.getNcols() << endl;

	cplex.solve();

	cout << "Solution status = " << cplex.getStatus() << endl;
	cout << "Objectvie value = " << cplex.getObjValue() << endl;
	for (auto var : primal_node_vars)
		cout << var.second.getName() << "\t" << cplex.getValue(var.second) << endl;
	for (auto var : source_node_vars)
		cout << var.second.getName() << "\t" << cplex.getValue(var.second) << endl;
	for (auto var : partition_node_vars)
		cout << var.second.getName() << "\t" << cplex.getValue(var.second) << endl;
	for (auto var : partition_flow_vars)
		cout << var.second.getName() << "\t" << cplex.getValue(var.second) << endl;
	for (auto var : multi_flow_vars)
		cout << var.second.getName() << "\t" << cplex.getValue(var.second) << endl;
	for (auto var : edge_vars)
		cout << var.second.getName() << "\t" << cplex.getValue(var.second) << endl;

	//elapsed_time  = cplex.getCplexTime() - start_time;
	//elapsed_ticks = cplex.getDetTime() - start_ticks;
}

void SmpSolver::solveLP_Steiner()
{
	LOG << "--------Solving LP for..." << endl;
	double start_time = cplex.getCplexTime();
	double start_ticks = cplex.getDetTime();

	bool separated = true;
	ncuts = 0;
	while (separated)
	{
		separated = false;
		LOG << "Time elapsed: " << cplex.getCplexTime() - start_time << endl;
		LOG << "Set new time limit of: " << time_limit - (cplex.getCplexTime() - start_time) << endl;
		cplex.setParam(IloCplex::TiLim, max(0.0, time_limit - (cplex.getCplexTime() - start_time))); //TILIM - elapsed time
		LOG << "Solve updated LP" << endl;
		cplex.solve();
		if (cplex.getStatus() != IloAlgorithm::Optimal)
			break;
		LOG << "Current obj value: " << cplex.getObjValue() << endl;

		IloEnv env = cplex.getEnv();
		IloNumArray val = IloNumArray(env, edge_vars.size());
		cplex.getValues(val, x_vararray);

		pair<NODE_PAIR, INDEX> pair_ij_k;
		SUB_Graph subG;
		map<pair<NODE_PAIR, INDEX>, double>xSol;
		for (auto k : G->p_set())
		{
			subG = G->get_subgraph()[k];
			pair_ij_k.second = k;
			for (auto arc : subG.arcs())
			{
				pair_ij_k.first.first = arc.first;
				pair_ij_k.first.second = arc.second;
				xSol[pair_ij_k] = val[x_varindex[pair_ij_k]];
			}
		}

		/* Build the cut */
		vector<IloExpr> cutLhs, cutRhs;
		vector<IloRange> cons;
		vector<double> violation;

		if (!separate_sc_Steiner(env, xSol, G, edge_vars, cutLhs, cutRhs, violation))
			seperate_min_cut_Steiner(env, xSol, G, edge_vars, cutLhs, cutRhs, violation, Steiner_root, primal_node_vars);

		// Only need to get the max_cuts maximally-violated inequalities
		vector<int> p(violation.size()); /* vector with indices */
		iota(p.begin(), p.end(), 0);     /* increasing */
		bool sorted = false;

		int attempts;
		if (max_cuts < 0)
			attempts = violation.size();
		else
		{
			attempts = min(max_cuts, int(violation.size()));
			partial_sort(p.begin(), p.begin() + attempts, p.end(),[&](int i, int j) 
			{ return violation[i] > violation[j]; });/* sort indices according to violation */
			sorted = true;
		}

		for (unsigned int i = 0; i < attempts; ++i)
		{
			if (violation[p[i]] >= TOL)
			{
				LOG << "Adding user cut for the " << i + 1 << "-th maximally violated constraint. Violation: "
					<< violation[p[i]] << endl;
				try
				{
					LOG << (cutLhs[p[i]] >= cutRhs[p[i]]) << endl;
					model.add(cutLhs[p[i]] >= cutRhs[p[i]]);
					separated = true;
					++ncuts;
				}
				catch (IloException e)
				{
					cerr << "Cannot add cut" << endl;
				}
			}
			else /* sorted, so no further violated ineq exist */
				if (sorted)
					break;
		}
		for (unsigned int i = 0; i < cutLhs.size(); ++i)
		{
			cutLhs[i].end();
			cutRhs[i].end();
		}
		val.end();
	}

	elapsed_time = cplex.getCplexTime() - start_time;
	elapsed_ticks = cplex.getDetTime() - start_ticks;
}

/* Single commodity flow formulation */
void SmpSolver::build_problem_scf()
{
	/*****************/
	/* Add variables */
	/*****************/

	cout << "Begin to build problem scf..." << endl;
	cout << "Begin to add variables..." << endl;

	IloEnv env = model.getEnv();

	char var_name[255];
	char con_name[255];
	int V_k_size;
	SUB_Graph subG;
	map<INDEX, NODE_SET> V_k_set = G->v_set();
	map<INDEX, NODE_SET> T_k_set = G->t_set();
	pair<NODE, INDEX> pair_i_k;
	pair<NODE_PAIR, INDEX> pair_ij_k;
	map<INDEX, NODE> root;
	NODE_PAIR pair_source_root;
	IloNumVar temp_var;

	// Add varibales x_i^k ,x_s_k, y_ij^k,
	for (auto k : G->p_set())
	{
		//cout <<  "k = " << k << endl;
		V_k_size = static_cast<int>(V_k_set[k].size());
		subG = G->get_subgraph()[k];

		// Add  x_i^k (binary), for each node in G[V_k]:
		pair_i_k.second = k;
		for (auto i : subG.nodes())
		{
			IloNumVar var;
			snprintf(var_name, 255, "x_%d^%d", i, k);
			if (T_k_set[k].find(i) != T_k_set[k].end())
			{
				var = IloNumVar(env, 1, 1, IloNumVar::Int, var_name);
			}

			else
			{
				var = IloNumVar(env, 0, 1, IloNumVar::Bool, var_name);
			}
			pair_i_k.first = i;
			partition_node_vars[pair_i_k] = var;
			printInfo(var);
		}

		// Add x_s_k (float [0, |V_k_set|]), for each k in P_index_set:
		snprintf(var_name, 255, "x_s%d", k);
		temp_var = IloNumVar(env, 0, V_k_size, IloNumVar::Float, var_name);
		source_node_vars[k] = temp_var;
		printInfo(temp_var);

		// for each T_k, choose a root r_k
		auto firstElement = T_k_set.at(k).begin();
		root[k] = *firstElement;
		cout << "r" << k << " = " << root[k] << endl;

		// Add y_ij^k (float), for arc (s_k, r_k) and each arc in G[V_k]:
		pair_ij_k.second = k;
		pair_source_root.first = 0;
		pair_source_root.second = root[k];
		pair_ij_k.first = pair_source_root;
		// Add y_(s_k,r_k)
		snprintf(var_name, 255, "y_s%dr%d", k, k);
		temp_var = IloNumVar(env, 0.0, IloInfinity, IloNumVar::Float, var_name);
		partition_flow_vars[pair_ij_k] = temp_var;
		printInfo(temp_var);

		// Add  y_ij^k
		for (auto &arc : subG.arcs())
		{
			IloNumVar var;
			snprintf(var_name, 255, "y_%d%d^%d", arc.first, arc.second, k);
			var = IloNumVar(env, 0.0, IloInfinity, IloNumVar::Float, var_name);
			pair_ij_k.first = arc;
			partition_flow_vars[pair_ij_k] = var;
			printInfo(var);
		}
	}

	/*******************/
	/* Add constraints */
	/*******************/

	cout << "Begin to add constraints..." << endl;

	// Constraint��5��: x_i >= x_i^k k \in P, i \in T_k
	for (auto i : G->v_total())
	{
		cout << "constraint(5)" << endl;
		for (auto k : G->nodes_of_v().at(i))
		{
			pair_i_k.first = i;
			pair_i_k.second = k;
			snprintf(con_name, 255, "%s >= %s (5)",
					 primal_node_vars[i].getName(),
					 partition_node_vars[pair_i_k].getName());
			model.add(primal_node_vars[i] >= partition_node_vars[pair_i_k]).setName(con_name);
			cout << con_name << endl;
		}
	}

	// Constraints (6) - (9)
	for (auto k : G->p_set())
	{
		// Set variables's k part
		V_k_size = static_cast<int>(V_k_set.at(k).size());
		subG = G->get_subgraph()[k];
		pair_i_k.second = k;
		pair_ij_k.second = k;
		pair_source_root.second = root[k];

		// Constraint (6): flow conservation at s_k
		//  y_(s_k, r_k) + x_s_k = |V_k|
		cout << "constraint(6)" << endl;
		pair_source_root.first = 0;
		pair_source_root.second = root[k];
		pair_ij_k.first = pair_source_root;

		snprintf(con_name, 255, "%s + %s >= %d (6)",
				 partition_flow_vars[pair_ij_k].getName(),
				 source_node_vars[k].getName(),
				 V_k_size);
		model.add(
				 partition_flow_vars[pair_ij_k] +
					 source_node_vars[k] >=
				 V_k_size)
			.setName(con_name);
		cout << con_name << endl;

		snprintf(con_name, 255, "%s + %s <= %d (6)",
				 partition_flow_vars[pair_ij_k].getName(),
				 source_node_vars[k].getName(),
				 V_k_size);
		model.add(
				 partition_flow_vars[pair_ij_k] +
					 source_node_vars[k] <=
				 V_k_size)
			.setName(con_name);
		cout << con_name << endl;

		// Constraint (7): y_ij^k <= |V_k| * x_j^k
		cout << "constraint(7)" << endl;
		for (auto &arc : subG.arcs())
		{
			pair_ij_k.first = arc;
			pair_i_k.first = arc.second;

			snprintf(con_name, 255, "%s <= %d * %s (7)",
					 partition_flow_vars[pair_ij_k].getName(),
					 V_k_size,
					 partition_node_vars[pair_i_k].getName());
			model.add(
					 partition_flow_vars[pair_ij_k] <=
					 V_k_size * partition_node_vars[pair_i_k])
				.setName(con_name);
			cout << con_name << endl;
		}

		// Constraint (8): flow conservation at x_j^k
		// (k-flow)->j = x_j_k + y->(k-flow)

		cout << "constraint(8)" << endl;
		for (auto j : subG.nodes())
		{
			IloExpr flow_in_j(env);
			IloExpr flow_out_j(env);
			pair_ij_k.first.second = j;
			char flow_in_j_name[255];
			memset(flow_in_j_name, 0, 255);
			for (auto i : subG.adj_nodes_list().at(j))
			{

				cout << "in" << j << ": " << i << endl;
				pair_ij_k.first.first = i;
				flow_in_j += partition_flow_vars[pair_ij_k];
				strcat(flow_in_j_name, partition_flow_vars[pair_ij_k].getName());
				strcat(flow_in_j_name, " + ");
			}
			// consider r_k and (s_k, r_k)
			if (j == root[k])
			{
				cout << "in" << j << ": "
					 << "s" << k << endl;
				pair_ij_k.first = pair_source_root;
				flow_in_j += partition_flow_vars[pair_ij_k];
				strcat(flow_in_j_name, partition_flow_vars[pair_ij_k].getName());
			}
			flow_in_j.setName(flow_in_j_name);

			pair_ij_k.first.first = j;
			char flow_out_j_name[255];
			memset(flow_out_j_name, 0, 255);
			for (auto i : subG.adj_nodes_list().at(j))
			{
				cout << "out" << j << ": " << i << endl;
				pair_ij_k.first.second = i;
				flow_out_j += partition_flow_vars[pair_ij_k];
				strcat(flow_out_j_name, partition_flow_vars[pair_ij_k].getName());
				strcat(flow_out_j_name, " + ");
			}
			flow_out_j.setName(flow_out_j_name);

			pair_i_k.first = j;

			snprintf(con_name, 255, "%s >= %s + %s (8)",
					 flow_in_j.getName(),
					 flow_out_j.getName(),
					 partition_node_vars[pair_i_k].getName());
			model.add(
					 flow_in_j >=
					 flow_out_j +
						 partition_node_vars[pair_i_k])
				.setName(con_name);
			cout << con_name << endl;

			snprintf(con_name, 255, "%s <= %s + %s (8)",
					 flow_in_j.getName(),
					 flow_out_j.getName(),
					 partition_node_vars[pair_i_k].getName());
			model.add(
					 flow_in_j <=
					 flow_out_j +
						 partition_node_vars[pair_i_k])
				.setName(con_name);
			cout << con_name << endl;
		}

		// Constraint (9): total system flow = flow of (0, r_k)
		// sum(x_j^k) = y_(s_k,r_k)
		cout << "constraint(9)" << endl;
		IloExpr selected_node_size(env);
		for (auto i : subG.nodes())
		{
			pair_i_k.first = i;
			selected_node_size += partition_node_vars[pair_i_k];
			cout << partition_node_vars[pair_i_k].getName() << " + ";
		}
		pair_ij_k.first = pair_source_root;

		cout << " = " << partition_flow_vars[pair_ij_k].getName() << endl;

		snprintf(con_name, 255, "selected_node_size_%d", k);
		selected_node_size.setName(con_name);

		snprintf(con_name, 255, "%s >= %s (9)",
				 selected_node_size.getName(),
				 partition_flow_vars[pair_ij_k].getName());
		model.add(
				 selected_node_size >=
				 partition_flow_vars[pair_ij_k])
			.setName(con_name);
		cout << con_name << endl;

		snprintf(con_name, 255, "%s <= %s (9)",
				 selected_node_size.getName(),
				 partition_flow_vars[pair_ij_k].getName());
		model.add(
				 selected_node_size <=
				 partition_flow_vars[pair_ij_k])
			.setName(con_name);
		cout << con_name << endl;
	}
}

/* Multi-commodity flow formulation */
void SmpSolver::build_problem_mcf()
{
	/*****************/
	/* Add variables */
	/*****************/

	cout << "Begin to build problem mcf..." << endl;
	cout << "Begin to add variables..." << endl;

	IloEnv env = model.getEnv();

	char var_name[255];
	char con_name[255];
	int V_k_size;
	SUB_Graph subG;
	map<INDEX, NODE_SET> V_k_set = G->v_set();
	map<INDEX, NODE_SET> T_k_set = G->t_set();
	pair<NODE, INDEX> pair_i_k;
	pair<NODE_PAIR, NODE_PAIR> pair_ij_km;
	map<INDEX, NODE> root;
	IloNumVar temp_var;

	// Add variables x_i^k,
	for (auto k : G->p_set())
	{
		cout << endl
			 << "Add " << k << " partition graph varaible..." << endl;
		V_k_size = static_cast<int>(V_k_set[k].size());
		subG = G->get_subgraph()[k];
		int subG_nodesNum = subG.nodes().size();

		// Add  x_i^k (binary), for each node in G[V_k]:
		pair_i_k.second = k;
		for (auto i : subG.nodes())
		{
			IloNumVar var;
			snprintf(var_name, 255, "x_%d^%d", i, k);
			if (T_k_set[k].find(i) != T_k_set[k].end())
				var = IloNumVar(env, 1, 1, IloNumVar::Int, var_name);
			else
				var = IloNumVar(env, 0, 1, IloNumVar::Int, var_name);
			pair_i_k.first = i;
			partition_node_vars[pair_i_k] = var;
			printInfo(var);
		}

		// For each T_k, choose a root r_k
		auto firstElement = T_k_set[k].begin();
		root[k] = *firstElement;
		cout << "r" << k << " = " << root[k] << endl;

		// Add y_ij_km (float)
		pair_ij_km.second.first = k;
		for (auto m : subG.nodes())
		{
			if (m == root[k])
				continue;
			pair_ij_km.second.second = m;
			for (auto arc : subG.arcs())
			{
				IloNumVar var;
				snprintf(var_name, 255, "x_%d,%d^%d,%d", arc.first, arc.second, k, m);
				var = IloNumVar(env, 0.0, IloInfinity, IloNumVar::Float, var_name);
				pair_ij_km.first = arc;
				multi_flow_vars[pair_ij_km] = var;
				printInfo(var);
			}
		}
	}

	/*******************/
	/* Add constraints */
	/*******************/

	cout << "Begin to add constraints..." << endl;

	// Add constraint (13-1)
	for (auto i : G->v_total())
	{
		cout << "constraint(13-1)" << endl;
		for (auto k : G->nodes_of_v().at(i))
		{
			pair_i_k.first = i;
			pair_i_k.second = k;
			snprintf(con_name, 255, "%s >= %s (5)", primal_node_vars[i].getName(),
					 partition_node_vars[pair_i_k].getName());
			model.add(primal_node_vars[i] >= partition_node_vars[pair_i_k]).setName(con_name);
			cout << con_name << endl;
		}
	}

	// Add constraint (14-21)
	cout << "Add constraint (14-16)" << endl;
	for (auto k : G->p_set())
	{
		cout << endl
			 << "For partition: " << k << endl;
		V_k_size = static_cast<int>(V_k_set[k].size());
		subG = G->get_subgraph()[k];
		pair_ij_km.second.first = k;

		// Add constraint (15-18)
		for (auto m : subG.nodes())
		{
			if (m == root[k])
				continue;
			pair_ij_km.second.second = m;

			// Partition for con_15： y_ir^k_km
			string con_15_name = "";
			IloRange con_15 = IloRange(env, 0, 0);
			for (auto i : subG.adj_nodes_list().at(root[k]))
			{
				pair_ij_km.first.first = i;
				pair_ij_km.first.second = root[k];
				con_15.setLinearCoef(multi_flow_vars[pair_ij_km], 1);
				con_15_name = con_15_name + " + " + multi_flow_vars[pair_ij_km].getName();
			}
			con_15.setName(con_15_name.c_str());
			model.add(con_15);
			cout << "con_15_name: " << con_15_name << endl;

			// Partition for con_16 ： y_im_km
			string con_16_name = "";
			pair_i_k.first = m;
			pair_i_k.second = k;
			IloNumVar coeff = partition_node_vars[pair_i_k];
			IloExpr flow_in_m(env);
			con_16_name = "";
			for (auto i : subG.adj_nodes_list().at(m))
			{
				pair_ij_km.first.first = i;
				pair_ij_km.first.second = m;
				flow_in_m += multi_flow_vars[pair_ij_km];
				con_16_name = con_16_name + " + " + multi_flow_vars[pair_ij_km].getName();
			}
			snprintf(con_name, 255, "con_16: %s <= %s", con_16_name.c_str(),
				partition_node_vars[pair_i_k].getName());
			cout << con_name << endl;
			model.add(flow_in_m <= coeff).setName(con_name);

			snprintf(con_name, 255, "con_16: %s >= %s", con_16_name.c_str(),
				partition_node_vars[pair_i_k].getName());
			cout << con_name << endl;
			model.add(flow_in_m >= coeff).setName(con_name);

			// Partition for con_17： y_mj_km
			string con_17_name = "";
			IloRange con_17 = IloRange(env, 0, 0);
			for (auto j : subG.adj_nodes_list().at(m))
			{
				pair_ij_km.first.first = m;
				pair_ij_km.first.second = j;
				con_17.setLinearCoef(multi_flow_vars[pair_ij_km], 1);
				con_17_name = con_17_name + " + " + multi_flow_vars[pair_ij_km].getName();
			}
			con_17.setName(con_17_name.c_str());
			model.add(con_17);
			cout << "con_17_name: " << con_17_name << endl;

			//Partition for con_18: sigma{y_iq_km} = sigma{y_kj_km}
			string con_18_name_left = "", con_18_name_right = "";
			IloExpr flow_in_q(env);
			IloExpr flow_out_q(env);
			for (auto q : subG.nodes())
			{
				if (q == m || q == root[k])
					continue;
				for (auto i : subG.adj_nodes_list().at(q))
				{
					pair_ij_km.first.first = i;
					pair_ij_km.first.second = q;
					flow_in_q += multi_flow_vars[pair_ij_km];
					con_18_name_left = con_18_name_left + "+" + multi_flow_vars[pair_ij_km].getName();
				}

				for (auto j : subG.adj_nodes_list().at(q))
				{
					pair_ij_km.first.first = q;
					pair_ij_km.first.second = j;
					flow_out_q += multi_flow_vars[pair_ij_km];
					con_18_name_right = con_18_name_right + "+" + multi_flow_vars[pair_ij_km].getName();
				}
			}
			snprintf(con_name, 255, "con18: %s <= %s",
				con_18_name_left.c_str(), con_18_name_right.c_str());
			cout << con_name << endl;
			model.add(flow_in_q <= flow_out_q).setName(con_name);

			snprintf(con_name, 255, "con18: %s >= %s",
				con_18_name_left.c_str(), con_18_name_right.c_str());
			cout << con_name << endl;
			model.add(flow_in_q >= flow_out_q).setName(con_name);
		}

		// Add constraint (19-20)
		pair_i_k.second = k;
		for (auto m : subG.nodes())
		{
			if (m == root[k])
				continue;
			pair_ij_km.second.second = m;
			for (auto arc : subG.arcs())
			{
				pair_ij_km.first.first = arc.first;
				pair_ij_km.first.second = arc.second;

				// constraint 18
				pair_i_k.first = arc.first;
				snprintf(con_name, 255, "con_19: %s <= %s", multi_flow_vars[pair_ij_km].getName(),
						 partition_node_vars[pair_i_k].getName());
				cout << con_name << endl;
				model.add(multi_flow_vars[pair_ij_km] <= partition_node_vars[pair_i_k]).setName(con_name);

				// constraint 19
				pair_i_k.first = arc.second;
				snprintf(con_name, 255, "con_20: %s <= %s", multi_flow_vars[pair_ij_km].getName(),
						 partition_node_vars[pair_i_k].getName());
				cout << con_name << endl;
				model.add(multi_flow_vars[pair_ij_km] <= partition_node_vars[pair_i_k]).setName(con_name);
			}
		}
	}

	cout << "end" << endl;
}

/* Steiner forest formulation */
void SmpSolver::build_problem_steiner()
{
	/*****************/
	/* Add variables */
	/*****************/

	cout << "Begin to build problem Steiner..." << endl;
	cout << "Begin to add variables..." << endl;

	IloEnv env = model.getEnv();

	char var_name[255];
	char con_name[255];
	int V_k_size;
	SUB_Graph subG;
	map<INDEX, NODE_SET> V_k_set = G->v_set();
	map<INDEX, NODE_SET> T_k_set = G->t_set();
	pair<NODE, INDEX> pair_i_k;
	pair<NODE_PAIR, INDEX> pair_ij_k;
	IloNumVar temp_var;

	//Add varaible y_ij_k
	int idx = 0;
	x_vararray = IloNumVarArray(env);
	for (auto k : G->p_set())
	{
		pair_ij_k.second = k;
		subG = G->get_subgraph()[k];
		for (auto &arc : subG.arcs())
		{
			IloNumVar var;
			snprintf(var_name, 255, "y_%d%d_%d", arc.first, arc.second, k);
			pair_ij_k.first = arc;
			var = IloNumVar(env, 0, 1, IloNumVar::Int, var_name);
			edge_vars[pair_ij_k] = var;
			x_vararray.add(var);
			x_varindex[pair_ij_k] = idx++;
			model.add(var);
			printInfo(var);
		}

		// For each T_k, choose a root r_k
		auto firstElement = T_k_set[k].begin();
		Steiner_root[k] = *firstElement;
		cout << "r" << k << " = " << Steiner_root[k] << endl;
	}

	/*******************/
	/* Add constraints */
	/*******************/
	
	cout << "Begin to Add the Constraint..." << endl;
	for (auto k : G->p_set())
	{
		cout << "For part " << k << "..." << endl;
		subG = G->get_subgraph()[k];

		// Partition for con_23: x_j >= sigma{y_ij_k}
		for (auto j : subG.nodes())
		{
			string con_23_name = "";
			IloExpr sigma_vars(env);
			pair_ij_k.second = k;
			for (auto i : subG.adj_nodes_list().at(j))
			{
				pair_ij_k.first.first = i;
				pair_ij_k.first.second = j;
				sigma_vars += edge_vars[pair_ij_k];
				con_23_name = con_23_name + " + " + edge_vars[pair_ij_k].getName();
			}
			snprintf(con_name, 255, "con_23: %s >= %s", primal_node_vars[j].getName(),
				con_23_name.c_str());
			cout << con_name << endl;
			model.add(primal_node_vars[j] >= sigma_vars);
		}

		// Partition for con_24: sigma{y_it_k = 1}
		for (auto t : T_k_set[k])
		{
			pair_ij_k.second = k;

			string con_24_name = "";
			IloRange con_24 = IloRange(env, 1, 1);
			for (auto i : subG.adj_nodes_list().at(t))
			{
				pair_ij_k.first.first = i;
				pair_ij_k.first.second = t;
				con_24.setLinearCoef(edge_vars[pair_ij_k], 1);
				con_24_name = con_24_name + " + " + edge_vars[pair_ij_k].getName();
			}
			con_24.setName(con_24_name.c_str());
			model.add(con_24);
			cout << "con_24: " << con_24_name << "= 1" << endl;
		}

		// Partition for con_25: sigma{y_iq_k <= 1}
		for (auto q : V_k_set[k])
		{
			if (std::find(T_k_set[k].begin(), T_k_set[k].end(), q) != T_k_set[k].end())
				continue;
			pair_ij_k.second = k;

			string con_25_name = "";
			IloRange con_25 = IloRange(env, 0, 1);
			for (auto i : subG.adj_nodes_list().at(q))
			{
				pair_ij_k.first.first = i;
				pair_ij_k.first.second = q;
				con_25.setLinearCoef(edge_vars[pair_ij_k], 1);
				con_25_name = con_25_name + " + " + edge_vars[pair_ij_k].getName();
			}
			con_25.setName(con_25_name.c_str());
			model.add(con_25);
			cout << "con_25: 0 <= " << con_25_name << " <= 1" << endl;
		}

		// Partition for con_26: y_ij_k + y_ji_k <= 1;
		for (auto arc : subG.arcs())
		{
			string con_26_name = "";
			pair_ij_k.second = k;
			NODE i = arc.first, j = arc.second;
			IloRange con_26 = IloRange(env, 0, 1);

			pair_ij_k.first = make_pair(i, j);
			con_26.setLinearCoef(edge_vars[pair_ij_k], 1);
			con_26_name = con_26_name + edge_vars[pair_ij_k].getName();

			pair_ij_k.first = make_pair(j, i);
			con_26.setLinearCoef(edge_vars[pair_ij_k], 1);
			con_26_name = con_26_name + " + " + edge_vars[pair_ij_k].getName();

			con_26.setName(con_26_name.c_str());
			cout << "con_26: " << con_26_name << " = 1" << endl;
			model.add(con_26);
		}
		cout << endl;
	}
	cout << "end" << endl;
}

/* node separator formulation */
void SmpSolver::build_problem_ns()
{
}

/* frequent used function */
void SmpSolver::printInfo(IloNumVar var)
{
	printf("add %s: (%f, %f), type: %d\n", var.getName(), var.getLb(), var.getUb(), var.getType());
}
