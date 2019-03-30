/*
file : separation.h
*/

#include "smp.h"
#include "graph.h"
#include "type.h"
#include "separation.h"
#include "unionfind.h"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <lemon/smart_graph.h>
#include <lemon/list_graph.h>
#include <lemon/concepts/maps.h>
#include <lemon/connectivity.h>
#include <lemon/preflow.h>
#include <lemon/time_measure.h>

#define LOG if(false) cerr
//#define LOG cout
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

void build_support_graph_ns(ListDigraph& support_graph, map<NODE, ListNode>& v_nodes, map<ListNode, NODE>&rev_nodes,
                            const map<pair<NODE, INDEX>, double>&xSol, std::shared_ptr<Graph> G, INDEX k)
{
	SUB_Graph subG = G->get_subgraph()[k];
	ListNode a, b;
	pair<NODE, INDEX>pair_i_k;

	for (NODE i : subG.nodes())
	{
		pair_i_k.first = i;
		pair_i_k.second = k;
		if (xSol.at(pair_i_k) > TOL && v_nodes.count(i) == 0)
		{
			a = support_graph.addNode();
			v_nodes[i] = a;
			rev_nodes[a] = i;
		}
		LOG << "added NODE: " << i << endl;
	}

	for (auto& arc : subG.arcs())
	{
		NODE u = arc.first;
		NODE v = arc.second;
		pair<NODE, INDEX>pair_u_k = make_pair(u, k);
		pair<NODE, INDEX>pair_v_k = make_pair(v, k);
		if (xSol.at(pair_u_k) > TOL && xSol.at(pair_v_k) > TOL)
		{
			support_graph.addArc(v_nodes[u], v_nodes[v]);
			support_graph.addArc(v_nodes[v], v_nodes[u]);
			LOG << "added arc: " << u << " " << v << endl;
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

void build_cap_graph_ns(ListDigraph& cap_graph, ListDigraph::ArcMap<double>& x_capacities, map<NODE, pair<ListNode, ListNode>>& v_nodes,
                        map<ListNode, NODE>& rev_nodes, const map<pair<NODE, INDEX>, double>&xSol, std:: shared_ptr<Graph>G, INDEX k,
                        const map<INDEX, NODE>& ns_root)
{
	pair<NODE, INDEX>pair_i_k;
	map<INDEX, NODE_SET> T_k_set = G->t_set();
	SUB_Graph subG = G->get_subgraph()[k];
	const double INF = subG.nodes().size() + 1;

	ListNode a, b;
	ListArc arc, rev_arc;
	ListNode_Pair list_node_pair;

	//Add Node and Arc
	for (NODE i : subG.nodes())
	{
		if (v_nodes.count(i) == 0) {
			a = cap_graph.addNode();
			v_nodes[i] = make_pair(a, a);
			rev_nodes[a] = i;
		}
	}

	for (auto& Arc : subG.arcs())
	{
		NODE u = Arc.first;
		NODE v = Arc.second;
		arc = cap_graph.addArc(v_nodes[u].first, v_nodes[v].first);
		x_capacities[arc] = INF;
		LOG << "added arc: " << u << " " << v;
		LOG << " with capacity: " << INF << endl;
	}

	//split node
	for (NODE i : subG.nodes())
	{
		if (i == ns_root.at(k))
			continue;
		ListNode new_node = cap_graph.split(v_nodes[i].first, false);
		v_nodes[i].second = new_node;
		rev_nodes[new_node] = i;

		pair_i_k.first = i;
		pair_i_k.second = k;
		arc = cap_graph.addArc(v_nodes[i].first, v_nodes[i].second);
		x_capacities[arc] = xSol.at(pair_i_k);

		LOG << "added arc: " << i << "' " << i << "'' ";
		LOG << " with capacity: " << xSol.at(pair_i_k) << endl;
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


bool seperate_sc_ns(IloEnv masterEnv, const map<pair<NODE, INDEX>, double>& xSol, std::shared_ptr<Graph>G, const map<pair<NODE, INDEX>, IloNumVar>& partition_node_vars,
                    vector<IloExpr>& cutLhs, vector<IloExpr>& cutRhs, vector<double>& violation, const map<INDEX, NODE>& ns_root)
{
	bool ret = false;
	pair<NODE, INDEX>pair_i_k;

	for (auto k : G->p_set())
	{
		cout << "For part : " << k << endl << endl;
		pair_i_k.second = k;

		// Build Support subGraph
		ListDigraph support_graph;
		map<NODE, ListNode>v_nodes;
		map<ListNode, NODE>rev_nodes;
		SUB_Graph subG = G->get_subgraph()[k];
		build_support_graph_ns(support_graph, v_nodes, rev_nodes, xSol, G, k);

		// Search for strong components
		ListDigraph::NodeMap<int> nodemap(support_graph);
		int components = stronglyConnectedComponents(support_graph, nodemap);

		cout << "Number of components is : " << components << endl;

		vector<int> cardinality(components, 0);
		vector<double> value_comp(components, 0);
		vector<NODE_SET> comp_set(components);

		// Add the divide the nodes into different component
		int root_comp;
		for (ListDigraph::NodeIt i(support_graph); i != INVALID; ++i)
		{
			LOG << "--------------- node " << rev_nodes[i] << endl;
			int comp = nodemap[i];
			if (cardinality[comp] == 0)
				cardinality[comp]++;
			if (rev_nodes[i] == ns_root.at(k))
				root_comp = comp;
			comp_set[comp].insert(rev_nodes[i]);
		}

		// Add the new initial graph with unionfind set
		NODE_SET root_adj_nodes;
		UnionFind<NODE>forest(subG.nodes());
		map<NODE, bool>reached;
		for (auto& arc : subG.arcs())
		{
			NODE u = arc.first;
			NODE v = arc.second;
			bool u_selected = v_nodes.count(u);
			bool v_selected = v_nodes.count(v);
			reached[u] = true;
			reached[v] = true;

			LOG << "for pair:" << u << " and  " << v << ", " << u_selected << " " << v_selected << endl;

			if (u_selected || v_selected)
			{
				if (u_selected && v_selected && nodemap[v_nodes[u]] == root_comp)
					continue;
				// one node is adj to the source
				else if (u_selected && nodemap[v_nodes[u]] == root_comp)
					root_adj_nodes.insert(v);
				else if (v_selected && nodemap[v_nodes[v]] == root_comp)
					root_adj_nodes.insert(u);
				// selected but not adj to the source
				else
				{
					// add arc to unionfind
					if (forest.find_set(u) != forest.find_set(v))
					{
						forest.join(u, v);
						LOG << "join :" << u << " and " << v << endl;
					}
				}
			}

			else if (forest.find_set(u) != forest.find_set(v))
				forest.join(u, v);
		}

		// Perform the check procedure (whether s and t is connected)
		for (auto target_set : comp_set)
		{
			IloExpr newCutLhs(masterEnv);
			IloExpr newCutRhs(masterEnv);
			double newCutValue = 0;
			double newViolation = 0;
			double totvalue = 1;

			auto firstElement = target_set.begin();
			auto t = *firstElement;

			if (nodemap[v_nodes[t]] == root_comp)
				continue;

			for (auto s : root_adj_nodes)
			{
				if (reached[s] && reached[t] && forest.find_set(s) == forest.find_set(t))
				{
					cout << "union " << s << "and" << t << endl;
					pair_i_k.first = s;
					newCutLhs += (partition_node_vars.at(pair_i_k));
					cout << partition_node_vars.at(pair_i_k).getName() << endl;
					newCutValue += xSol.at(pair_i_k);
				}
				else
					continue;
			}

			IloNumVar temp_var = IloNumVar(masterEnv, 1, 1, IloNumVar::Float);
			newCutRhs += temp_var;

			newViolation = 1.0 - newCutValue;

			if (newCutValue < 1)
			{
				cutLhs.push_back(newCutLhs);
				cutRhs.push_back(newCutRhs);
				violation.push_back(newViolation);

				LOG << "lhs: " << newCutLhs << endl;
				LOG << "rhs: " << newCutRhs << endl;
				LOG << "violation: " << newViolation << endl;

				if (newViolation > TOL)
					ret = true;
			}
		}
	}
	return ret;
}

/*  Min cut seperation for NS  */
bool seperate_min_cut_ns(IloEnv masterEnv, const map<pair<NODE, INDEX>, double>&xSol, std::shared_ptr<Graph>G, const map<pair<NODE, INDEX>, IloNumVar>& partition_node_vars,
                         vector<IloExpr>& cutLhs, vector<IloExpr>& cutRhs, vector<double>& violation, const map<INDEX, NODE>& ns_root)
{
	bool ret = false;
	pair<NODE, INDEX>pair_i_k;
	map<INDEX, NODE_SET> V_k_set = G->v_set();
	map<INDEX, NODE_SET> T_k_set = G->t_set();

	cutLhs = vector<IloExpr>();
	cutRhs = vector<IloExpr>();
	violation = vector<double>();

	for (auto k : G->p_set())
	{
		LOG << "Ran min-cut..." << endl;
		pair_i_k.second = k;

		ListDigraph cap_graph;
		ListDigraph::ArcMap<double>x_capacities(cap_graph);
		map<NODE, pair<ListNode, ListNode>> v_nodes;
		map<ListNode, NODE> rev_nodes;
		SUB_Graph subG = G->get_subgraph()[k];

		build_cap_graph_ns(cap_graph, x_capacities, v_nodes, rev_nodes, xSol, G, k, ns_root);

		LOG << "Built graph..." << endl;

		IloExpr newCutLhs;
		IloExpr newCutRhs;
		double newViolation;
		double min_cut_value;

		for (auto q : subG.t_set())
		{
			if (q == ns_root.at(k))
				continue;

			Preflow<ListDigraph, ListDigraph::ArcMap<double>>min_cut(cap_graph, x_capacities, v_nodes[ns_root.at(k)].first, v_nodes[q].first);
			min_cut.runMinCut();
			min_cut_value = min_cut.flowValue();

			LOG << q << endl;
			LOG << "Min-cut " << min_cut_value << endl;

			if (min_cut_value < 1)
			{
				newCutLhs = IloExpr(masterEnv);
				newCutRhs = IloExpr(masterEnv);
				newViolation = 1 - min_cut_value;

				/*for (ListDigraph::NodeIt i(cap_graph); i != INVALID; ++i)
				{
					if (min_cut.minCut(i))
					{
						bool is_cut_arc = false;

						ListNode a = v_nodes[rev_nodes[i]].first;
						ListNode b = v_nodes[rev_nodes[i]].second;
						if (a == i)
							is_cut_arc = !min_cut.minCut(b);
						else if (b == i)
							is_cut_arc = !min_cut.minCut(a);
						if (is_cut_arc) {
							LOG << "find cur arc: " << rev_nodes[i] << endl;
							pair_i_k.first = rev_nodes[i];
							newCutLhs += (partition_node_vars.at(pair_i_k));
						}
					}
				}*/

				for (NODE i : subG.nodes())
				{
					if (i == ns_root.at(k))
						continue;
					ListNode a = v_nodes[i].first;
					ListNode b = v_nodes[i].second;

					if ( (min_cut.minCut(a) && !min_cut.minCut(b))
					        || (!min_cut.minCut(a) && min_cut.minCut(b)))
					{
						LOG << "find cur arc: " << i << endl;
						pair_i_k.first = i;
						newCutLhs += (partition_node_vars.at(pair_i_k));
					}

				}
				IloNumVar temp_var = IloNumVar(masterEnv, 1, 1, IloNumVar::Int);
				newCutRhs += temp_var;

				cutLhs.push_back(newCutLhs);
				cutRhs.push_back(newCutRhs);
				violation.push_back(newViolation);

				LOG << "node " << q << endl;
				LOG << "cut " << cutLhs.size() << endl;
				LOG << "flowValue " << min_cut_value << endl;
				LOG << "lhs: " << newCutLhs << endl;
				LOG << "rhs: " << newCutRhs << endl;
			}
		}
	}
	ret = cutLhs.size() > 0 ? 1 : 0;
	return ret;
}
