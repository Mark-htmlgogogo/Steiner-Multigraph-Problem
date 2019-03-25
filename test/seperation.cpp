/*
file : separation.h
*/

#include "smp.h"
#include "graph.h"
#include "type.h"
#include "separation.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <lemon/smart_graph.h>
#include <lemon/concepts/maps.h>
#include <lemon/connectivity.h>
#include <lemon/preflow.h>
#include <lemon/time_measure.h>

#define LOG if(false) cerr
#define TOL 0.001

using namespace std;
using namespace lemon;

void build_support_graph_Steiner(SmartDigraph& support_graph, map<NODE, LemonNode>& v_nodes, map<LemonNode, NODE>& rev_nodes,
                                 const map<pair<NODE_PAIR, INDEX>, double>& xSol, std::shared_ptr<Graph> G, INDEX k)
{
	SUB_Graph subG = G->get_subgraph()[k];
	LemonNode a, b;
	for (NODE i : subG.nodes())
	{
		for (NODE j : subG.adj_nodes_list().at(i))
		{
			pair<NODE_PAIR, INDEX> pair_ij_k;
			pair_ij_k.first.first = i;
			pair_ij_k.first.second = j;
			pair_ij_k.second = k;

			if (xSol.at(pair_ij_k) > TOL)
			{
				if (v_nodes.count(i) == 0)
				{
					a = support_graph.addNode();
					v_nodes[i] = a;
					rev_nodes[a] = i;
				}
				if (v_nodes.count(j) == 0)
				{
					b = support_graph.addNode();
					v_nodes[j] = b;
					rev_nodes[b] = j;
				}
				support_graph.addArc(v_nodes[i], v_nodes[j]);
				support_graph.addArc(v_nodes[j], v_nodes[i]);
				LOG << "added arc: " << i << " " << j << endl;
			}
		}
	}
}

void build_cap_graph_Steiner(SmartDigraph& cap_graph, SmartDigraph::ArcMap<double>& x_capacities, map<NODE, LemonNode>& v_nodes, map<LemonNode, NODE>& rev_nodes,
                             const map<pair<NODE_PAIR, INDEX>, double>& xSol, std::shared_ptr<Graph>G, INDEX k)
{
	pair<NODE_PAIR, INDEX>pair_ij_k;
	pair_ij_k.second = k;
	SUB_Graph subG = G->get_subgraph()[k];
	LemonNode a, b;
	LemonArc arc, rev_arc;
	for (NODE i : subG.nodes())
	{
		for (NODE j : subG.adj_nodes_list().at(i))
		{
			if (v_nodes.count(i) == 0)
			{
				a = cap_graph.addNode();
				v_nodes[i] = a;
				rev_nodes[a] = i;
			}
			if (v_nodes.count(j) == 0)
			{
				b = cap_graph.addNode();
				v_nodes[j] = b;
				rev_nodes[b] = j;
			}
			arc = cap_graph.addArc(v_nodes[i], v_nodes[j]);
			rev_arc = cap_graph.addArc(v_nodes[j], v_nodes[i]);

			pair_ij_k.first.first = i;
			pair_ij_k.first.second = j;
			x_capacities[arc] = xSol.at(pair_ij_k);
			LOG << "added arc: " << i << " " << j;
			LOG << " with capacity: " << xSol.at(pair_ij_k) << endl;

			pair_ij_k.first.first = j;
			pair_ij_k.first.second = i;
			x_capacities[rev_arc] = xSol.at(pair_ij_k);

			LOG << "added arc: " << j << " " << i;
			LOG << " with capacity: " << xSol.at(pair_ij_k) << endl;
		}
	}
}

void build_cap_graph_ns(ListGraph& cap_graph, ListGraph::EdgeMap<double>& x_capacities, map<NODE, pair<ListNode, ListNode>>& v_nodes,
                        map<pair<ListNode, ListNode>, NODE>& rev_nodes, const map<pair<NODE, INDEX>, double>&xSol, std:: shared_ptr<Graph>G, INDEX k, map<INDEX, NODE>& ns_root)
{
	pair<NODE, INDEX>pair_i_k;
	map<INDEX, NODE_SET> T_k_set = G->t_set();
	SUB_Graph subG = G->get_subgraph()[k];

	ListNode a, b;
	ListEdge arc, rev_arc;
	ListNode_Pair list_node_pair;

	// Add NODE first
	for (NODE i : subG.nodes()) {
		if (i == ns_root[k])
			continue;

		// if NODE i is a terminal node
		if (std::find(T_k_set[k].begin(), T_k_set[k].end(), i) == T_k_set[k].end() && v_nodes.count(i) != 0) {
			a = cap_graph.addNode();
			list_node_pair = make_pair(a, a);
			v_nodes[i] = list_node_pair;
			rev_nodes[list_node_pair] = i;
		}

		// if NDOE i is not a terminal node
		else {
			a = cap_graph.addNode();
			b = cap_graph.addNode();
			list_node_pair = make_pair(a, b);
			v_nodes[i] = list_node_pair;
			rev_nodes[list_node_pair] = i;

			//Add edge for terminal node directly
			pair_i_k.first = i;
			pair_i_k.second = k;
			arc = cap_graph.addEdge(v_nodes[i].first, v_nodes[i].second);
			x_capacities[arc] = xSol.at(pair_i_k);
			LOG << "added arc: " << i << "' " << i << "' ";
			LOG << "with capacity: " << xSol.at(pair_i_k) << endl;
		}
	}

	// Add Arc for those edge initially exist in the sub_graph
	for (NODE_PAIR edge : subG.arcs()) {
		NODE i = edge.first, j = edge.second;
		arc = cap_graph.addEdge(v_nodes[i].first, v_nodes[i].second);
		x_capacities[arc] = INF;
		LOG << "added arc: " << i << " " << j << " ";
		LOG << "with capacity: " << INF << endl;
	}
}

/*  Strong Component separation for Steiner  */
bool separate_sc_Steiner(IloEnv masterEnv, const map<pair<NODE_PAIR, INDEX>, double>& xSol, std::shared_ptr<Graph>G,
                         const map<pair<NODE_PAIR, INDEX>, IloNumVar>& edge_vars, vector<IloExpr>& cutLhs, vector<IloExpr>& cutRhs, vector<double>& violation)
{
	bool ret = false;
	pair<NODE_PAIR, INDEX>pair_ij_k;

	for (auto k : G->p_set())
	{
		cout << "For part : " << k << endl << endl;
		pair_ij_k.second = k;

		/* Build Support subGraph */
		SmartDigraph support_graph;
		map<NODE, LemonNode>v_nodes;
		map<LemonNode, NODE>rev_nodes;
		SUB_Graph subG = G->get_subgraph()[k];
		build_support_graph_Steiner(support_graph, v_nodes, rev_nodes, xSol, G, k);

		/* Search for strong components */
		SmartDigraph::NodeMap<int> nodemap(support_graph);
		int components = stronglyConnectedComponents(support_graph, nodemap);

		cout << "Number of components is : " << components << endl;

		cutLhs = vector<IloExpr>(components);
		cutRhs = vector<IloExpr>(components);
		violation = vector<double>(components);
		vector<int> cardinality(components, 0);
		vector<double> out_degree(components, 0);
		vector<double> max_node_degree(components, 0);
		//vector<NODE> max_node(components, 0);
		SmartDigraph::NodeMap<double> node_out_degree(support_graph, 0);

		for (SmartDigraph::NodeIt i(support_graph); i != INVALID; ++i)
		{
			LOG << "--------------- node " << rev_nodes[i] << endl;
			int comp = nodemap[i];
			if (0 == cardinality[comp]++)
			{
				LOG << "Initialized lhs" << endl;
				cutLhs[comp] = IloExpr(masterEnv);
			}
			for (NODE j : subG.adj_nodes_list().at(rev_nodes[i]))
			{
				pair_ij_k.first.first = rev_nodes[i];
				if (v_nodes.count(j) == 0 || comp != nodemap[v_nodes[j]])
				{
					pair_ij_k.first.second = j;
					node_out_degree[i] += (xSol.at(pair_ij_k));
					cutLhs[comp] += (edge_vars.at(pair_ij_k));
				}
			}
			LOG << "degree: " << node_out_degree[i] << endl;
			if (node_out_degree[i] >= max_node_degree[comp] + TOL)
			{
				LOG << "add rhs " << endl;
				max_node_degree[comp] = node_out_degree[i];
				cutRhs[comp] = IloExpr(masterEnv);
				for (NODE j : subG.adj_nodes_list().at(rev_nodes[i]))
				{
					pair_ij_k.first.first = rev_nodes[i];
					pair_ij_k.first.second = j;
					cutRhs[comp] += (edge_vars.at(pair_ij_k));
				}
			}
		}

		for (int i = 0; i < components; ++i)
		{
			LOG << "component " << i << ": " << cardinality[i] << endl;
			LOG << "delta(S) " << out_degree[i] << endl;
			LOG << "max(delta(i)) " << max_node_degree[i] << endl;
			LOG << "lhs: " << cutLhs[i] << endl;
			LOG << "rhs: " << cutRhs[i] << endl;
			violation[i] = max_node_degree[i] - out_degree[i];
			if (violation[i] >= TOL)
				ret = true;
		}
	}
	return ret;
}

/*  Min cut separation for Steiner  */
bool seperate_min_cut_Steiner(IloEnv masterEnv, const map<pair<NODE_PAIR, INDEX>, double>& xSol, std::shared_ptr<Graph>G,
                              const map<pair<NODE_PAIR, INDEX>, IloNumVar>& edge_vars, vector<IloExpr>& cutLhs, vector<IloExpr>& cutRhs, vector<double>& violation,
                              const map<INDEX, NODE>& root, const map<NODE, IloNumVar>& primal_node_vars)
{
	bool ret = false;
	pair<NODE_PAIR, INDEX>pair_ij_k;
	cutLhs = vector<IloExpr>();
	cutRhs = vector<IloExpr>();
	violation = vector<double>();

	for (auto k : G->p_set())
	{
		LOG << "Ran min-cut..." << endl;
		pair_ij_k.second = k;

		// Build graph with x values as capacities
		SmartDigraph cap_graph;
		SmartDigraph::ArcMap<double> x_capacities(cap_graph);
		map<NODE, LemonNode> v_nodes;
		map<LemonNode, NODE> rev_nodes;
		SUB_Graph subG = G->get_subgraph()[k];
		build_cap_graph_Steiner(cap_graph, x_capacities, v_nodes, rev_nodes, xSol, G, k);

		LOG << "Built graph..." << endl;

		IloExpr newCutLhs;
		IloExpr newCutRhs;
		double newViolation;
		double min_cut_value;

		//rv_cut >= delta(v)
		for (NODE q : subG.nodes())
		{
			if (q == root.at(k))
				continue;
			Preflow<SmartDigraph, SmartDigraph::ArcMap<double>>min_cut(cap_graph, x_capacities, v_nodes[root.at(k)], v_nodes[q]);
			min_cut.runMinCut();
			min_cut_value = min_cut.flowValue();

			// Compute the out degree of v
			double node_out_degree = 0.0;
			for (auto j : subG.adj_nodes_list().at(q))
			{
				pair_ij_k.first.first = j;
				pair_ij_k.first.second = q;
				node_out_degree += xSol.at(pair_ij_k);
			}

			LOG << q << endl;
			LOG << "Min-cut " << min_cut_value << endl;
			LOG << "Node degree " << node_out_degree << endl;

			if (node_out_degree > min_cut_value)
			{
				newCutLhs = IloExpr(masterEnv);
				newCutRhs = IloExpr(masterEnv);
				newViolation = node_out_degree - min_cut_value;
				for (SmartDigraph::NodeIt i(cap_graph); i != INVALID; ++i)
				{
					if (min_cut.minCut(i))
						for (NODE j : subG.adj_nodes_list().at(rev_nodes[i]))
							if (v_nodes.count(j) == 0 || !min_cut.minCut(v_nodes[j]))
							{
								pair_ij_k.first.first = rev_nodes[i];
								pair_ij_k.first.second = j;
								newCutLhs += (edge_vars.at(pair_ij_k));
							}
				}

				for (auto j : subG.adj_nodes_list().at(q))
				{
					pair_ij_k.first.first = j;
					pair_ij_k.first.second = q;
					pair_ij_k.second = k;
					newCutRhs += (edge_vars.at(pair_ij_k));
				}

				cutLhs.push_back(newCutLhs);
				cutRhs.push_back(newCutRhs);
				violation.push_back(newViolation);

				LOG << "node " << q << endl;
				LOG << "cut " << cutLhs.size() << endl;
				LOG << "delta(S) " << min_cut_value << endl;
				LOG << "delta(v) " << node_out_degree << endl;
				LOG << "lhs: " << newCutLhs << endl;
				LOG << "rhs: " << newCutRhs << endl;
			}
		}
	}
	ret = cutLhs.size() > 0 ? 1 : 0;
	return ret;
}

/*  Min cut seperation for NS  */
bool seperate_min_cut_ns(IloEnv masterEnv, const map<pair<NODE, INDEX>, double>&xSol, std::shared_ptr<Graph>,
                         const map<pair<NODE, INDEX>, IloNumVar>& partition_node_vars, vector<IloExpr>& cutLhs, vector<IloExpr>& cutRhs, vector<double>& violation,
                         const map<INDEX, NODE>& ns_root)
{
	bool ret = false;
}
