#pragma once

#include "smp.h"
#include "graph.h"
#include "type.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <lemon/smart_graph.h>
#include <lemon/concepts/maps.h>
#include <lemon/connectivity.h>
#include <lemon/preflow.h>
#include <lemon/time_measure.h>

using namespace std;
using namespace lemon;

typedef SmartDigraph::Node LemonNode;
typedef SmartDigraph::Arc LemonArc;

void build_support_graph(SmartDigraph& support_graph, map<NODE, LemonNode>& v_nodes, map<LemonNode, NODE>& rev_nodes,
	const map<NODE_PAIR, double>& xSol, std::shared_ptr<Graph>, INDEX k);

void build_cap_graph(SmartDigraph& cap_graph, SmartDigraph::ArcMap<double>& x_capacities, map<NODE, LemonNode>& v_nodes, map<LemonNode, NODE>& rev_nodes,
	const map<NODE_PAIR, double>& xSol, std::shared_ptr<Graph>, INDEX k);

bool separate_sc_Steiner(IloEnv masterEnv, const map<pair<NODE_PAIR,INDEX>, double>& xSol, std::shared_ptr<Graph>,
	const map<pair<NODE_PAIR, INDEX>, IloNumVar>& edge_vars, vector<IloExpr>& cutLhs, vector<IloExpr>& cutRhs, vector<double>& violation);

bool seperate_min_cut_Steiner(IloEnv masterEnv, const map<pair<NODE_PAIR, INDEX>, double>& xSol, std::shared_ptr<Graph>,
	const map<pair<NODE_PAIR, INDEX>, IloNumVar>& edge_vars, vector<IloExpr>& cutLhs, vector<IloExpr>& cutRhs, vector<double>& violation,
	const map<INDEX, NODE>& root, const map<NODE, IloNumVar>& primal_node_vars);
