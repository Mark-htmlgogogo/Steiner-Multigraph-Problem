/*
file : separation.h
*/

#include <lemon/concepts/maps.h>
#include <lemon/connectivity.h>
#include <lemon/list_graph.h>
#include <lemon/preflow.h>
#include <lemon/smart_graph.h>
#include <lemon/time_measure.h>

#include <algorithm>
#include <fstream>
#include <iostream>

#include "separation.h"

#define LOG \
    if (false) cerr
//#define LOG cout
#define TOL 0.001
#define rep(i, a, b) for (int i = (a); i < (b); i++)

using namespace std;
using namespace lemon;

extern CutPool cutpool;

void build_support_graph_Steiner(
    SmartDigraph& support_graph, map<NODE, LemonNode>& v_nodes,
    map<LemonNode, NODE>& rev_nodes,
    const map<pair<NODE_PAIR, INDEX>, double>& xSol, std::shared_ptr<Graph> G,
    INDEX k) {
    SUB_Graph subG = G->get_subgraph()[k];
    LemonNode a, b;
    for (NODE i : subG.nodes()) {
        for (NODE j : subG.adj_nodes_list().at(i)) {
            pair<NODE_PAIR, INDEX> pair_ij_k;
            pair_ij_k.first.first = i;
            pair_ij_k.first.second = j;
            pair_ij_k.second = k;

            if (xSol.at(pair_ij_k) > TOL) {
                if (v_nodes.count(i) == 0) {
                    a = support_graph.addNode();
                    v_nodes[i] = a;
                    rev_nodes[a] = i;
                }
                if (v_nodes.count(j) == 0) {
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

void build_support_graph_ns(ListDigraph& support_graph,
                            map<NODE, ListNode>& v_nodes,
                            map<ListNode, NODE>& rev_nodes,
                            const map<pair<NODE, INDEX>, double>& xSol,
                            std::shared_ptr<Graph> G, INDEX k) {
    SUB_Graph subG = G->get_subgraph()[k];
    ListNode a, b;
    pair<NODE, INDEX> pair_i_k;

    for (NODE i : subG.nodes()) {
        pair_i_k.first = i;
        pair_i_k.second = k;
        if (xSol.at(pair_i_k) > TOL && v_nodes.count(i) == 0) {
            a = support_graph.addNode();
            v_nodes[i] = a;
            rev_nodes[a] = i;
        }
        LOG << "added NODE: " << i << endl;
    }
    for (auto& arc : subG.arcs()) {
        NODE u = arc.first;
        NODE v = arc.second;
        pair<NODE, INDEX> pair_u_k = make_pair(u, k);
        pair<NODE, INDEX> pair_v_k = make_pair(v, k);
        if (xSol.at(pair_u_k) > TOL && xSol.at(pair_v_k) > TOL) {
            support_graph.addArc(v_nodes[u], v_nodes[v]);
            support_graph.addArc(v_nodes[v], v_nodes[u]);
            LOG << "added arc: " << u << " " << v << endl;
        }
    }
}

void build_cap_graph_Steiner(SmartDigraph& cap_graph,
                             SmartDigraph::ArcMap<double>& x_capacities,
                             map<NODE, LemonNode>& v_nodes,
                             map<LemonNode, NODE>& rev_nodes,
                             const map<pair<NODE_PAIR, INDEX>, double>& xSol,
                             std::shared_ptr<Graph> G, INDEX k) {
    pair<NODE_PAIR, INDEX> pair_ij_k;
    pair_ij_k.second = k;
    SUB_Graph subG = G->get_subgraph()[k];
    LemonNode a, b;
    LemonArc arc, rev_arc;
    for (NODE i : subG.nodes()) {
        for (NODE j : subG.adj_nodes_list().at(i)) {
            if (v_nodes.count(i) == 0) {
                a = cap_graph.addNode();
                v_nodes[i] = a;
                rev_nodes[a] = i;
            }
            if (v_nodes.count(j) == 0) {
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

// Initial build method for NS
// void build_cap_graph_ns(ListDigraph& cap_graph,
//	ListDigraph::ArcMap<double>& x_capacities,
//	map<NODE, pair<ListNode, ListNode>>& v_nodes,
//	map<ListNode, NODE>& rev_nodes,
//	const map<pair<NODE, INDEX>, double>& xSol,
//	std::shared_ptr<Graph> G, INDEX k,
//	const map<INDEX, NODE>& ns_root) {
//
//	pair<NODE, INDEX> pair_i_k;
//	map<INDEX, NODE_SET> T_k_set = G->t_set();
//	SUB_Graph subG = G->get_subgraph()[k];
//	const double INF = 10000000000000.0;
//
//	ListNode a, b;
//	ListArc arc, rev_arc;
//	ListNode_Pair list_node_pair;
//
//	// Add Node and Arc
//	for (NODE i : subG.nodes()) {
//		if (v_nodes.count(i) == 0) {
//			a = cap_graph.addNode();
//			v_nodes[i] = make_pair(a, a);
//			rev_nodes[a] = i;
//		}
//	}
//
//	for (auto& Arc : subG.arcs()) {
//		NODE u = Arc.first;
//		NODE v = Arc.second;
//		arc = cap_graph.addArc(v_nodes[u].first, v_nodes[v].first);
//		x_capacities[arc] = INF;
//		LOG << "added arc: " << u << " " << v;
//		LOG << " with capacity: " << INF << endl;
//	}
//
//	// split node
//	for (NODE i : subG.nodes()) {
//		if (i == ns_root.at(k)) continue;
//		ListNode new_node = cap_graph.split(v_nodes[i].first, false);
//		v_nodes[i].second = new_node;
//		rev_nodes[new_node] = i;
//
//		pair_i_k.first = i;
//		pair_i_k.second = k;
//		arc = cap_graph.addArc(v_nodes[i].first, v_nodes[i].second);
//		x_capacities[arc] = xSol.at(pair_i_k);
//
//		LOG << "added arc: " << i << "' " << i << "'' ";
//		LOG << " with capacity: " << xSol.at(pair_i_k) << endl;
//	}
//
//	return;
//}

// NS mincut graph pre-construct varaible
map<INDEX, ListDigraph> ns_mincut_capgraph;
map<INDEX, map<NODE, pair<ListNode, ListNode>>> ns_mincut_v_nodes;
map<INDEX, map<ListNode, NODE>> ns_mincut_rev_nodes;
map<INDEX, map<NODE, ListArc>> ns_mincut_split_arc;
// NS mincut graph pre-construct procedure
void generate_ns_mincut_graph(std::shared_ptr<Graph> G,
                              const map<INDEX, NODE>& ns_root) {
    // clear the varaible
    ns_mincut_capgraph.clear();
    ns_mincut_v_nodes.clear();
    ns_mincut_rev_nodes.clear();
    ns_mincut_split_arc.clear();

    for (auto k : G->p_set()) {
        ListNode a, b;
        ListArc arc, rev_arc;
        SUB_Graph subG = G->get_subgraph()[k];

        // Add Node and Arc
        for (NODE i : subG.nodes()) {
            if (ns_mincut_v_nodes[k].count(i) == 0) {
                a = ns_mincut_capgraph[k].addNode();
                ListNode_Pair LP = make_pair(a, a);
                ns_mincut_v_nodes[k].insert(make_pair(i, LP));
                ns_mincut_rev_nodes[k].insert(make_pair(a, i));
            }
        }

        for (auto& Arc : subG.arcs()) {
            NODE u = Arc.first;
            NODE v = Arc.second;
            arc = ns_mincut_capgraph[k].addArc(ns_mincut_v_nodes[k][u].first,
                                               ns_mincut_v_nodes[k][v].first);
        }

        for (NODE i : subG.nodes()) {
            if (i == ns_root.at(k)) continue;
            ListNode new_node = ns_mincut_capgraph[k].split(
                ns_mincut_v_nodes[k][i].first, false);
            ns_mincut_v_nodes[k][i].second = new_node;
            ns_mincut_rev_nodes[k][new_node] = i;

            arc = ns_mincut_capgraph[k].addArc(ns_mincut_v_nodes[k][i].first,
                                               ns_mincut_v_nodes[k][i].second);

            ns_mincut_split_arc[k][i] = arc;
        }
    }

    return;
}

// new build method for NS mincut
void build_cap_graph_ns(std::shared_ptr<Graph> G,
                        ListDigraph::ArcMap<double>& x_capacities, INDEX k,
                        const map<INDEX, NODE>& ns_root,
                        const map<pair<NODE, INDEX>, double>& xSol) {
    SUB_Graph subG = G->get_subgraph()[k];
    pair<NODE, INDEX> pair_i_k;
    ListArc arc;

    for (NODE i : subG.nodes()) {
        if (i == ns_root.at(k)) continue;
        pair_i_k.first = i;
        pair_i_k.second = k;
        x_capacities[ns_mincut_split_arc[k][i]] = xSol.at(pair_i_k);
    }

    return;
}

/*  Strong Component separation for Steiner  */
bool separate_sc_Steiner(
    IloEnv masterEnv, const map<pair<NODE_PAIR, INDEX>, double>& xSol,
    std::shared_ptr<Graph> G,
    const map<pair<NODE_PAIR, INDEX>, IloNumVar>& edge_vars,
    vector<IloExpr>& cutLhs, vector<IloExpr>& cutRhs,
    vector<double>& violation) {
    bool ret = false;
    pair<NODE_PAIR, INDEX> pair_ij_k;

    for (auto k : G->p_set()) {
        pair_ij_k.second = k;

        /* Build Support subGraph */
        SmartDigraph support_graph;
        map<NODE, LemonNode> v_nodes;
        map<LemonNode, NODE> rev_nodes;
        SUB_Graph subG = G->get_subgraph()[k];
        build_support_graph_Steiner(support_graph, v_nodes, rev_nodes, xSol, G,
                                    k);

        /* Search for strong components */
        SmartDigraph::NodeMap<int> nodemap(support_graph);
        int components = stronglyConnectedComponents(support_graph, nodemap);

        cutLhs = vector<IloExpr>(components);
        cutRhs = vector<IloExpr>(components);
        violation = vector<double>(components);
        vector<int> cardinality(components, 0);
        vector<double> out_degree(components, 0);
        vector<double> max_node_degree(components, 0);
        // vector<NODE> max_node(components, 0);
        SmartDigraph::NodeMap<double> node_out_degree(support_graph, 0);

        for (SmartDigraph::NodeIt i(support_graph); i != INVALID; ++i) {
            LOG << "--------------- node " << rev_nodes[i] << endl;
            int comp = nodemap[i];
            if (0 == cardinality[comp]++) {
                LOG << "Initialized lhs" << endl;
                cutLhs[comp] = IloExpr(masterEnv);
            }
            for (NODE j : subG.adj_nodes_list().at(rev_nodes[i])) {
                pair_ij_k.first.first = rev_nodes[i];
                if (v_nodes.count(j) == 0 || comp != nodemap[v_nodes[j]]) {
                    pair_ij_k.first.second = j;
                    node_out_degree[i] += (xSol.at(pair_ij_k));
                    cutLhs[comp] += (edge_vars.at(pair_ij_k));
                }
            }
            LOG << "degree: " << node_out_degree[i] << endl;
            if (node_out_degree[i] >= max_node_degree[comp] + TOL) {
                LOG << "add rhs " << endl;
                max_node_degree[comp] = node_out_degree[i];
                cutRhs[comp] = IloExpr(masterEnv);
                for (NODE j : subG.adj_nodes_list().at(rev_nodes[i])) {
                    pair_ij_k.first.first = rev_nodes[i];
                    pair_ij_k.first.second = j;
                    cutRhs[comp] += (edge_vars.at(pair_ij_k));
                }
            }
        }

        for (int i = 0; i < components; ++i) {
            LOG << "component " << i << ": " << cardinality[i] << endl;
            LOG << "delta(S) " << out_degree[i] << endl;
            LOG << "max(delta(i)) " << max_node_degree[i] << endl;
            LOG << "lhs: " << cutLhs[i] << endl;
            LOG << "rhs: " << cutRhs[i] << endl;
            violation[i] = max_node_degree[i] - out_degree[i];
            if (violation[i] >= TOL) ret = true;
        }
    }
    return ret;
}

/*  Min cut separation for Steiner  */
bool seperate_min_cut_Steiner(
    IloEnv masterEnv, const map<pair<NODE_PAIR, INDEX>, double>& xSol,
    std::shared_ptr<Graph> G,
    const map<pair<NODE_PAIR, INDEX>, IloNumVar>& edge_vars,
    vector<IloExpr>& cutLhs, vector<IloExpr>& cutRhs, vector<double>& violation,
    const map<INDEX, NODE>& root,
    const map<NODE, IloNumVar>& primal_node_vars) {
    bool ret = false;
    pair<NODE_PAIR, INDEX> pair_ij_k;
    cutLhs = vector<IloExpr>();
    cutRhs = vector<IloExpr>();
    violation = vector<double>();

    for (auto k : G->p_set()) {
        LOG << "Ran min-cut..." << endl;
        pair_ij_k.second = k;

        // Build graph with x values as capacities
        SmartDigraph cap_graph;
        SmartDigraph::ArcMap<double> x_capacities(cap_graph);
        map<NODE, LemonNode> v_nodes;
        map<LemonNode, NODE> rev_nodes;
        SUB_Graph subG = G->get_subgraph()[k];
        build_cap_graph_Steiner(cap_graph, x_capacities, v_nodes, rev_nodes,
                                xSol, G, k);

        LOG << "Built graph..." << endl;

        IloExpr newCutLhs;
        IloExpr newCutRhs;
        double newViolation;
        double min_cut_value;

        // rv_cut >= delta(v)
        for (NODE q : subG.nodes()) {
            if (q == root.at(k)) continue;
            Preflow<SmartDigraph, SmartDigraph::ArcMap<double>> min_cut(
                cap_graph, x_capacities, v_nodes[root.at(k)], v_nodes[q]);
            min_cut.runMinCut();
            min_cut_value = min_cut.flowValue();

            // Compute the out degree of v
            double node_out_degree = 0.0;
            for (auto j : subG.adj_nodes_list().at(q)) {
                pair_ij_k.first.first = j;
                pair_ij_k.first.second = q;
                node_out_degree += xSol.at(pair_ij_k);
            }

            LOG << q << endl;
            LOG << "Min-cut " << min_cut_value << endl;
            LOG << "Node degree " << node_out_degree << endl;

            if (node_out_degree >= min_cut_value + TOL) {
                newCutLhs = IloExpr(masterEnv);
                newCutRhs = IloExpr(masterEnv);
                newViolation = node_out_degree - min_cut_value;
                for (SmartDigraph::NodeIt i(cap_graph); i != INVALID; ++i) {
                    if (min_cut.minCut(i))
                        for (NODE j : subG.adj_nodes_list().at(rev_nodes[i]))
                            if (v_nodes.count(j) == 0 ||
                                !min_cut.minCut(v_nodes[j])) {
                                pair_ij_k.first.first = rev_nodes[i];
                                pair_ij_k.first.second = j;
                                newCutLhs += (edge_vars.at(pair_ij_k));
                            }
                }

                for (auto j : subG.adj_nodes_list().at(q)) {
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

//-----------------------------------NS-----------------------------------
int* vis;       // whether a elem has been visted in dfs
int* CompV;     // node belong to which comp in sol
int** CompS;    // comp s contains which v
int* CompSnum;  // size of comp s
int* CompAcrk;  // comp divided for determine node cut
int* ACrk;      // whether a node is Acrk

void dfs1(NODE now, std::shared_ptr<SUB_Graph> subG, int& k,
          const map<pair<NODE, INDEX>, double>& xSol, int color) {
    if (vis[now]) return;
    vis[now] = 1;

    pair<NODE, INDEX> pair_i_k;

    for (auto i : subG->adj_nodes_list().at(now)) {
        pair_i_k.first = i;
        pair_i_k.second = k;
        if (!vis[i] && xSol.at(pair_i_k) > TOL) {
            CompV[i] = color;
            CompS[color][CompSnum[color]] = i;
            CompSnum[color]++;
            dfs1(i, subG, k, xSol, color);
        }
    }
    return;
}

void dfs2(NODE now, std::shared_ptr<SUB_Graph> subG, int RootComp, int color) {
    if (vis[now]) return;
    vis[now] = 1;

    for (auto i : subG->adj_nodes_list().at(now)) {
        if (CompV[i] == RootComp) {
            // cout << i << " is in RootComp;  ";
            continue;
        } else if (ACrk[now] && ACrk[i]) {
            // cout << i << " is in Acrk;  ";
            continue;
        } else if (!vis[i]) {
            // cout << i << " ";
            CompAcrk[i] = color;
            dfs2(i, subG, RootComp, color);
        }
    }
    return;
}

/*  Strong Component separation for NS  */
bool seperate_sc_ns(
    IloEnv masterEnv, const map<pair<NODE, INDEX>, double>& xSol,
    std::shared_ptr<Graph> G,
    const map<pair<NODE, INDEX>, IloNumVar>& partition_node_vars,
    vector<IloExpr>& cutLhs, vector<IloExpr>& cutRhs, vector<double>& violation,
    const map<INDEX, NODE>& ns_root, int& lazy_sep_opt) {
    bool ret = false;
    pair<NODE, INDEX> pair_i_k;

    int T, N, GNsize = G->nodes().size() + 1;

    vis = (int*)malloc(int(GNsize) * sizeof(int));
    CompV = (int*)malloc(int(GNsize) * sizeof(int));
    CompAcrk = (int*)malloc(int(GNsize) * sizeof(int));
    ACrk = (int*)malloc(int(GNsize) * sizeof(int));

    for (auto k : G->p_set()) {
        pair_i_k.second = k;

        // SUB_Graph subG = G->get_subgraph()[k];
        std::shared_ptr<SUB_Graph> subG =
            std::make_shared<SUB_Graph>(G->get_subgraph()[k]);
        T = subG->t_set().size() + 1;
        N = subG->nodes().size() + 1;

        CompS = (int**)malloc(int(T) * sizeof(int**));
        for (int i = 0; i < T; i++) {
            CompS[i] = (int*)malloc(int(N) * sizeof(int));
            memset(CompS[i], 0, sizeof(int) * (N));
        }
        CompSnum = (int*)malloc(int(T) * sizeof(int));
        memset(vis, 0, sizeof(int) * GNsize);
        memset(CompV, 0, sizeof(int) * GNsize);
        memset(CompSnum, 0, sizeof(int) * (T));

        // find all the component in initial graph
        // component starts from 1 to n
        int components = 0;
        for (auto t : subG->t_set()) {
            if (!vis[t]) {
                components++;
                CompV[t] = components;
                CompS[components][CompSnum[components]] = t;
                CompSnum[components]++;
                dfs1(t, subG, k, xSol, components);
            }
        }

        // cout << "For Partition " << k << endl;
        // print
        // cout << "Sol is: " << endl;
        /* for (auto i : subG->nodes()) {
            pair_i_k.first = i;
            pair_i_k.second = k;
            cout << pair_i_k << " " << xSol.at(pair_i_k) << endl;
        }
        */
        // cout << "number of components: " << components << endl;
        /* for (int i = 1; i <= components; i++) {
            cout << "For components: " << i << endl;
            for (int j = 0; j < CompSnum[i]; j++) {
                cout << CompS[i][j] << " ";
            }
            cout << endl;
        } */

        // begin find node cut
        // lazy-sep-opt == 0 : on to multi
        int root_comp = CompV[ns_root.at(k)];
        int RootComp = !lazy_sep_opt ? root_comp : 1;
        int RootCompIter = !lazy_sep_opt ? root_comp : components;

        for (; RootComp <= RootCompIter; RootComp++) {
            // find A[Crk]
            set<int> nAcrk;
            memset(ACrk, 0, sizeof(int) * GNsize);
            for (int num = 0; num < CompSnum[RootComp]; num++) {
                NODE i = CompS[RootComp][num];
                for (auto j : subG->adj_nodes_list().at(i)) {
                    if (CompV[j] != RootComp) {
                        nAcrk.insert(j);
                        ACrk[j] = 1;
                    }
                }
            }

            // print
            /*  cout << "Root Comp is: " << RootComp << endl;
             cout << "root Comp is: " << root_comp << endl;
             cout << "nAcrk are: ";
             for (auto i : nAcrk) cout << i << " ";
             cout << endl; */

            memset(vis, 0, sizeof(int) * GNsize);
            memset(CompAcrk, 0, sizeof(int) * GNsize);
            int AcNum = 0;
            // rep(i, 1, 16) cout << vis[i];
            // cout << endl;
            for (auto i : nAcrk) {
                if (!vis[i]) {
                    // cout << i << " ";
                    CompAcrk[i] = ++AcNum;
                    dfs2(i, subG, RootComp, AcNum);
                }
                // cout << endl;
            }

            int TarComp = !lazy_sep_opt ? 1 : RootComp;
            for (; TarComp <= components; TarComp++) {
                if (TarComp == RootComp) continue;

                // Perform the check procedure (whether s and t is connected
                int t = CompS[TarComp][0];
                IloExpr newCutLhs(masterEnv);
                IloExpr newCutRhs(masterEnv);
                double newCutValue = 0;
                double newViolation = 0;
                double totvalue = 1.0;

                // cout << "cut added is: ";
                set<NODE> cutset;
                for (auto s : nAcrk) {
                    // int s = nAcrk[num];
                    if (CompAcrk[s] == CompAcrk[t]) {
                        pair_i_k.first = s;
                        pair_i_k.second = k;
                        newCutLhs += partition_node_vars.at(pair_i_k);
                        // cout << s << " ";
                        cutset.insert(s);
                        newCutValue += xSol.at(pair_i_k);
                    } else {
                        continue;
                    }
                }
                // cout << endl;

                newViolation = 1.0 - newCutValue;
                IloNumVar temp_var =
                    IloNumVar(masterEnv, 1, 1, IloNumVar::Float);
                newCutRhs += temp_var;

                if (newCutValue < 1 - TOL && cutset.size() != 0) {
                    cutLhs.push_back(newCutLhs);
                    cutRhs.push_back(newCutRhs);
                    violation.push_back(newViolation);

                    LOG << "lhs: " << newCutLhs << endl;
                    LOG << "rhs: " << newCutRhs << endl;
                    LOG << "violation: " << newViolation << endl;

                    if (newViolation >= TOL) ret = true;
                }

                cutpool.AddLhs(k, cutset);
                cutpool.AddViolation(k, newViolation);
                cutset.clear();
            }
        }

        for (int i = 0; i < T; i++) free(CompS[i]);
        free(CompS);
        free(CompSnum);
        // cout << sizeof(CompS) << endl;
        // cout << endl;
    }
    free(vis);
    free(CompV);
    free(CompAcrk);
    free(ACrk);

    //---------------------------------------------------------------------------
    /*    for (auto k : G->p_set()) {
           pair_i_k.second = k;

           // Build Support subGraph
           ListDigraph support_graph;
           map<NODE, ListNode> v_nodes;
           map<ListNode, NODE> rev_nodes;
           SUB_Graph subG = G->get_subgraph()[k];
           map<INDEX, NODE_SET> T_k_set = G->t_set();
           build_support_graph_ns(support_graph, v_nodes, rev_nodes, xSol, G,
       k);

           map<NODE, bool> NodeUsedInLemon;
           for (auto i : subG.nodes()) {
               NodeUsedInLemon[i] = false;
           }

           // Search for strongly connected components
           ListDigraph::NodeMap<int> node_comp_map(support_graph);

           int components = stronglyConnectedComponents(
               support_graph,
               node_comp_map);  // return the number of SCCs and map i to its
       SCC
                                // if there is only one SCC

           vector<int> cardinality(components, 0);
           vector<double> value_comp(components, 0);
           vector<NODE_SET> comp_set(components);
           vector<bool> CompHasTerminal(components, 0);

           // Nodes in each SCC: comp_set[comp]
           int root_comp;
           for (ListDigraph::NodeIt i(support_graph); i != INVALID; ++i) {
               int comp = node_comp_map[i];
               if (cardinality[comp] == 0) cardinality[comp]++;
               if (rev_nodes[i] == ns_root.at(k)) root_comp = comp;
               comp_set[comp].insert(rev_nodes[i]);
               NodeUsedInLemon[rev_nodes[i]] = true;
               if (subG.CheckNodeIsTerminal().at(rev_nodes[i])) {
                   CompHasTerminal[comp] = 1;
               }
           }

           // Enumerate the components set
           int RootComp = !lazy_sep_opt ? root_comp : 0;
           int RootCompIter = !lazy_sep_opt ? root_comp + 1 : components;

           for (; RootComp < RootCompIter; RootComp++) {
               if (CompHasTerminal[RootComp] == 0) continue;

               // Begin to search path
               set<NODE> root_adj_nodes;
               for (auto i : comp_set[RootComp]) {
                   for (auto j : subG.adj_nodes_list().at(i)) {
                       if (!v_nodes.count(j)) {
                           root_adj_nodes.insert(j);
                       }
                   }
               }

               // Add the arc between the different node.
               UnionFind<NODE> forest(subG.nodes());
               map<NODE, bool> reached;
               for (auto& arc : subG.arcs()) {
                   NODE u = arc.first;
                   NODE v = arc.second;

                   bool u_selected = v_nodes.count(u);
                   bool v_selected = v_nodes.count(v);
                   if ((u_selected && node_comp_map[v_nodes[u]] == RootComp) ||
                       (v_selected && node_comp_map[v_nodes[v]] == RootComp))
                       continue;
                   if (root_adj_nodes.count(u) && root_adj_nodes.count(v))
                       continue;
                   reached[u] = true;
                   reached[v] = true;

                   if (forest.find_set(u) != forest.find_set(v)) {
                       forest.join(u, v);
                   }
               }

               int TarComp = !lazy_sep_opt ? 0 : RootComp;

               for (; TarComp < components; TarComp++) {
                   if (CompHasTerminal[TarComp] == 0 || RootComp == TarComp)
                       continue;

                   // Perform the check procedure (whether s and t is connected
                   auto firstElement = comp_set[TarComp].begin();
                   auto t = *firstElement;
                   IloExpr newCutLhs(masterEnv);
                   IloExpr newCutRhs(masterEnv);
                   double newCutValue = 0;
                   double newViolation = 0;
                   double totvalue = 1;

                   set<NODE> cutset;
                   for (auto s : root_adj_nodes) {
                       if (reached[s] && reached[t] &&
                           forest.find_set(s) == forest.find_set(t)) {
                           pair_i_k.second = k;
                           pair_i_k.first = s;
                           newCutLhs += (partition_node_vars.at(pair_i_k));
                           newCutValue += xSol.at(pair_i_k);  // 0

                           cutset.insert(s);

                       } else
                           continue;
                   }

                   IloNumVar temp_var =
                       IloNumVar(masterEnv, 1, 1, IloNumVar::Float);
                   newCutRhs += temp_var;

                   newViolation = 1.0 - newCutValue;

                   if (newCutValue < 1 - TOL && cutset.size() != 0) {
                       cutLhs.push_back(newCutLhs);
                       cutRhs.push_back(newCutRhs);
                       violation.push_back(newViolation);

                       LOG << "lhs: " << newCutLhs << endl;
                       LOG << "rhs: " << newCutRhs << endl;
                       LOG << "violation: " << newViolation << endl;

                       if (newViolation >= TOL) ret = true;
                   }

                   cutpool.AddLhs(k, cutset);
                   cutpool.AddViolation(k, newViolation);
                   cutset.clear();
               }
           }
       }
        */

    ret = cutLhs.size() > 0 ? 1 : 0;
    return ret;
}

//---------------------------------DINIC------------------------------------
const double inf = 999999999.0;
const double EPS = 0.0001;
struct st {
    int u, v, next;
    double w;
};
st** edge;
st* DinicEdge;
int **head, *dis, *q, *work, *use, *t;
int GNsize, N, M, T;
double min(double a, double b) { return (b - a) > EPS ? a : b; }

void DinicAddUni(int u, int v, double w, int* head, st* edge, int& t) {
    edge[t].u = u;
    edge[t].v = v;
    edge[t].w = w;
    edge[t].next = head[u];
    head[u] = t;
    t++;

    edge[t].u = v;
    edge[t].v = u;
    edge[t].w = 0.0;
    edge[t].next = head[v];
    head[v] = t++;
}
int DinicBfs(int S, int T, st* edge, int* head) {
    int rear = 0;
    memset(dis, -1, sizeof(int) * (2 * N));
    // for (int i = 0; i < 2 * N; i++)cout << dis[i] << " ";

    dis[S] = 0;
    q[rear++] = S;
    for (int i = 0; i < rear; i++) {
        for (int j = head[q[i]]; j != -1; j = edge[j].next) {
            int v = edge[j].v;
            if (edge[j].w > EPS && dis[v] == -1) {
                dis[v] = dis[q[i]] + 1;
                q[rear++] = v;
                if (v == T) return 1;
            }
        }
    }
    return 0;
}
double DinicDfs(int S, double a, int T, st* edge) {
    if (S == T) return a;
    for (int& i = work[S]; i != -1; i = edge[i].next) {
        int v = edge[i].v;
        if (edge[i].w > EPS && dis[v] == dis[S] + 1) {
            double tt = DinicDfs(v, min(a, edge[i].w), T, edge);
            if (tt > EPS) {
                edge[i].w -= tt;
                edge[i ^ 1].w += tt;
                return tt;
            }
        }
    }
    return 0;
}
double Dinic(int S, int T, int* head, st* edge) {
    double ans = 0.0;
    while (DinicBfs(S, T, edge, head)) {
        memcpy(work, head, sizeof(int) * (2 * N));
        double tt = DinicDfs(S, inf, T, edge);
        while (tt > EPS) {
            ans += tt;
            tt = DinicDfs(S, inf, T, edge);
        }
    }
    return ans;
}
void FindCut(int S, int* head, st* edge) {
    // memset(use, 0, sizeof(use));
    use[S] = 1;
    for (int i = head[S]; i != -1; i = edge[i].next) {
        int v = edge[i].v;
        if (!use[v] && edge[i].w > EPS) FindCut(v, head, edge);
    }
}

vector<map<NODE_PAIR, int>> EdgeToNode;
void PreBuildGraph(std::shared_ptr<Graph> G, const map<INDEX, NODE>& ns_root) {
    int P = G->p_set().size() + 1;
    head = (int**)malloc(int(P) * sizeof(int**));
    edge = (st**)malloc(int(P) * sizeof(st**));
    t = (int*)malloc(int(P) * sizeof(int));
    memset(t, 0, sizeof(int) * int(P));
    EdgeToNode.resize(P);

    for (auto k : G->p_set()) {
        std::shared_ptr<SUB_Graph> subG =
            std::make_shared<SUB_Graph>(G->get_subgraph()[k]);

        int N = G->nodes().size() + 1;
        int M = subG->arcs().size() + 1;
        int T = N + M + 1;
        edge[k] = (st*)malloc(int(2 * T) * sizeof(st));
        head[k] = (int*)malloc(int(2 * N) * sizeof(int));
        memset(head[k], -1, sizeof(int) * int(2 * N));

        int n = G->nodes().size();
        // cout << "root is: " << ns_root.at(k) << endl;
        for (auto i : subG->nodes()) {
            if (i == ns_root.at(k)) {
                for (auto j : subG->adj_nodes_list().at(i)) {
                    DinicAddUni(i, j, inf, head[k], edge[k], t[k]);
                }
            } else {
                for (auto j : subG->adj_nodes_list().at(i)) {
                    DinicAddUni(i + n, j, inf, head[k], edge[k], t[k]);
                }
                EdgeToNode[k][NODE_PAIR(i, i + n)] = t[k];
                DinicAddUni(i, i + n, 0.0, head[k], edge[k], t[k]);
            }
        }

        /*cout << endl << "For Partition " << k << ": " << G->nodes().size() <<
        endl << endl; cout << t[k] << endl; int idx = 0; for (auto u :
        subG->nodes()) { for (int i = head[k][u]; i != -1; i = edge[k][i].next)
        { cout << u << " to " << " " << edge[k][i].v << " is " << edge[k][i].w
                                << endl;
                        idx++;
                }
                for (int i = head[k][u + n]; i != -1; i = edge[k][i].next)
                {
                        cout << u + n << " to " << " " << edge[k][i].v << " is "
        << edge[k][i].w
                                << endl;
                        idx++;
                }
        }
        cout << idx << endl << endl;*/
    }

    return;
}
bool BuildDinicGraph(std::shared_ptr<Graph> G, std::shared_ptr<SUB_Graph> subG,
                     const map<pair<NODE, INDEX>, double>& xSol, int k,
                     st* edge, int root) {
    bool flag = false;
    int n = G->nodes().size();
    for (auto i : subG->nodes()) {
        if (i == root) continue;
        NODE_PAIR pair_i_k = NODE_PAIR(i, k);
        int E = EdgeToNode[k][NODE_PAIR(i, i + n)];
        edge[E].w = xSol.at(pair_i_k);
        if (edge[E].w < (1 - EPS) && edge[E].w > EPS) flag = true;
        // cout << i << " " << xSol.at(pair_i_k) << endl;
    }
    /*int idx = 0;
    for (auto u : subG->nodes()) {
            for (int i = head[k][u]; i != -1; i = edge[i].next) {
                    cout << edge[i].u << " to " << " " << edge[i].v << " is " <<
                            edge[i].w
                            << endl;
                    idx++;
            }
            for (int i = head[k][u + n]; i != -1; i = edge[i].next) {
                    cout << edge[i].u << " to " << " " << edge[i].v << " is " <<
                            edge[i].w
                            << endl;
                    idx++;
            }
    }
    cout << idx << endl << endl;*/
    return flag;
}

bool SingleLazySeperate(
    IloEnv masterEnv, std::shared_ptr<SUB_Graph> subG,
    const map<pair<NODE, INDEX>, double>& xSol, int k,
    const map<pair<NODE, INDEX>, IloNumVar>& partition_node_vars,
    vector<IloExpr>& cutLhs, vector<IloExpr>& cutRhs, vector<double>& violation,
    const map<INDEX, NODE>& ns_root) {
    // cout << "---------------Single Lazy Seperate----------" << endl;
    bool ret = false;
    pair<NODE, INDEX> pair_i_k;

    int T, N;

    vis = (int*)malloc(int(GNsize) * sizeof(int));
    CompV = (int*)malloc(int(GNsize) * sizeof(int));
    CompAcrk = (int*)malloc(int(GNsize) * sizeof(int));
    ACrk = (int*)malloc(int(GNsize) * sizeof(int));

    pair_i_k.second = k;

    T = subG->t_set().size() + 1;
    N = subG->nodes().size() + 1;

    CompS = (int**)malloc(int(T) * sizeof(int**));
    for (int i = 0; i < T; i++) {
        CompS[i] = (int*)malloc(int(N) * sizeof(int));
        memset(CompS[i], 0, sizeof(int) * (N));
    }
    CompSnum = (int*)malloc(int(T) * sizeof(int));
    memset(vis, 0, sizeof(int) * GNsize);
    memset(CompV, 0, sizeof(int) * GNsize);
    memset(CompSnum, 0, sizeof(int) * (T));

    // find all the component in initial graph
    // component starts from 1 to n
    int components = 0;
    for (auto t : subG->t_set()) {
        if (!vis[t]) {
            components++;
            CompV[t] = components;
            CompS[components][CompSnum[components]] = t;
            CompSnum[components]++;
            dfs1(t, subG, k, xSol, components);
        }
    }

    // cout << "For Partition " << k << endl;
    // print
    // cout << "Sol is: " << endl;
    /* for (auto i : subG->nodes()) {
        pair_i_k.first = i;
        pair_i_k.second = k;
        cout << pair_i_k << " " << xSol.at(pair_i_k) << endl;
    }
    */
    // cout << "number of components: " << components << endl;
    /* for (int i = 1; i <= components; i++) {
        cout << "For components: " << i << endl;
        for (int j = 0; j < CompSnum[i]; j++) {
            cout << CompS[i][j] << " ";
        }
        cout << endl;
    } */

    // begin find node cut
    // lazy-sep-opt == 0 : on to multi
    int lazy_sep_opt = 0;
    int root_comp = CompV[ns_root.at(k)];
    int RootComp = !lazy_sep_opt ? root_comp : 1;
    int RootCompIter = !lazy_sep_opt ? root_comp : components;

    for (; RootComp <= RootCompIter; RootComp++) {
        // find A[Crk]
        set<int> nAcrk;
        memset(ACrk, 0, sizeof(int) * GNsize);
        for (int num = 0; num < CompSnum[RootComp]; num++) {
            NODE i = CompS[RootComp][num];
            for (auto j : subG->adj_nodes_list().at(i)) {
                if (CompV[j] != RootComp) {
                    nAcrk.insert(j);
                    ACrk[j] = 1;
                }
            }
        }

        // print
        /*  cout << "Root Comp is: " << RootComp << endl;
         cout << "root Comp is: " << root_comp << endl;
         cout << "nAcrk are: ";
         for (auto i : nAcrk) cout << i << " ";
         cout << endl; */

        memset(vis, 0, sizeof(int) * GNsize);
        memset(CompAcrk, 0, sizeof(int) * GNsize);
        int AcNum = 0;
        // rep(i, 1, 16) cout << vis[i];
        // cout << endl;
        for (auto i : nAcrk) {
            if (!vis[i]) {
                // cout << i << " ";
                CompAcrk[i] = ++AcNum;
                dfs2(i, subG, RootComp, AcNum);
            }
            // cout << endl;
        }

        int TarComp = !lazy_sep_opt ? 1 : RootComp;
        for (; TarComp <= components; TarComp++) {
            if (TarComp == RootComp) continue;

            // Perform the check procedure (whether s and t is connected
            int t = CompS[TarComp][0];
            IloExpr newCutLhs(masterEnv);
            IloExpr newCutRhs(masterEnv);
            double newCutValue = 0;
            double newViolation = 0;
            double totvalue = 1.0;

            // cout << "cut added is: ";
            set<NODE> cutset;
            for (auto s : nAcrk) {
                // int s = nAcrk[num];
                if (CompAcrk[s] == CompAcrk[t]) {
                    pair_i_k.first = s;
                    pair_i_k.second = k;
                    newCutLhs += partition_node_vars.at(pair_i_k);
                    // cout << s << " ";
                    cutset.insert(s);
                    newCutValue += xSol.at(pair_i_k);
                } else {
                    continue;
                }
            }
            // cout << endl;

            newViolation = 1.0 - newCutValue;
            IloNumVar temp_var = IloNumVar(masterEnv, 1, 1, IloNumVar::Float);
            newCutRhs += temp_var;

            if (newCutValue < 1 - TOL && cutset.size() != 0) {
                cutLhs.push_back(newCutLhs);
                cutRhs.push_back(newCutRhs);
                violation.push_back(newViolation);

                LOG << "lhs: " << newCutLhs << endl;
                LOG << "rhs: " << newCutRhs << endl;
                LOG << "violation: " << newViolation << endl;

                if (newViolation >= TOL) ret = true;
            }

            cutpool.AddLhs(k, cutset);
            cutpool.AddViolation(k, newViolation);
            cutset.clear();
        }
    }

    for (int i = 0; i < T; i++) free(CompS[i]);
    free(CompS);
    free(CompSnum);
    // cout << sizeof(CompS) << endl;
    // cout << endl;

    free(vis);
    free(CompV);
    free(CompAcrk);
    free(ACrk);

    //---------------------------------------------------------------------------
    /*    for (auto k : G->p_set()) {
           pair_i_k.second = k;

           // Build Support subGraph
           ListDigraph support_graph;
           map<NODE, ListNode> v_nodes;
           map<ListNode, NODE> rev_nodes;
           SUB_Graph subG = G->get_subgraph()[k];
           map<INDEX, NODE_SET> T_k_set = G->t_set();
           build_support_graph_ns(support_graph, v_nodes, rev_nodes, xSol, G,
       k);

           map<NODE, bool> NodeUsedInLemon;
           for (auto i : subG.nodes()) {
               NodeUsedInLemon[i] = false;
           }

           // Search for strongly connected components
           ListDigraph::NodeMap<int> node_comp_map(support_graph);

           int components = stronglyConnectedComponents(
               support_graph,
               node_comp_map);  // return the number of SCCs and map i to its
       SCC
                                // if there is only one SCC

           vector<int> cardinality(components, 0);
           vector<double> value_comp(components, 0);
           vector<NODE_SET> comp_set(components);
           vector<bool> CompHasTerminal(components, 0);

           // Nodes in each SCC: comp_set[comp]
           int root_comp;
           for (ListDigraph::NodeIt i(support_graph); i != INVALID; ++i) {
               int comp = node_comp_map[i];
               if (cardinality[comp] == 0) cardinality[comp]++;
               if (rev_nodes[i] == ns_root.at(k)) root_comp = comp;
               comp_set[comp].insert(rev_nodes[i]);
               NodeUsedInLemon[rev_nodes[i]] = true;
               if (subG.CheckNodeIsTerminal().at(rev_nodes[i])) {
                   CompHasTerminal[comp] = 1;
               }
           }

           // Enumerate the components set
           int RootComp = !lazy_sep_opt ? root_comp : 0;
           int RootCompIter = !lazy_sep_opt ? root_comp + 1 : components;

           for (; RootComp < RootCompIter; RootComp++) {
               if (CompHasTerminal[RootComp] == 0) continue;

               // Begin to search path
               set<NODE> root_adj_nodes;
               for (auto i : comp_set[RootComp]) {
                   for (auto j : subG.adj_nodes_list().at(i)) {
                       if (!v_nodes.count(j)) {
                           root_adj_nodes.insert(j);
                       }
                   }
               }

               // Add the arc between the different node.
               UnionFind<NODE> forest(subG.nodes());
               map<NODE, bool> reached;
               for (auto& arc : subG.arcs()) {
                   NODE u = arc.first;
                   NODE v = arc.second;

                   bool u_selected = v_nodes.count(u);
                   bool v_selected = v_nodes.count(v);
                   if ((u_selected && node_comp_map[v_nodes[u]] == RootComp) ||
                       (v_selected && node_comp_map[v_nodes[v]] == RootComp))
                       continue;
                   if (root_adj_nodes.count(u) && root_adj_nodes.count(v))
                       continue;
                   reached[u] = true;
                   reached[v] = true;

                   if (forest.find_set(u) != forest.find_set(v)) {
                       forest.join(u, v);
                   }
               }

               int TarComp = !lazy_sep_opt ? 0 : RootComp;

               for (; TarComp < components; TarComp++) {
                   if (CompHasTerminal[TarComp] == 0 || RootComp == TarComp)
                       continue;

                   // Perform the check procedure (whether s and t is connected
                   auto firstElement = comp_set[TarComp].begin();
                   auto t = *firstElement;
                   IloExpr newCutLhs(masterEnv);
                   IloExpr newCutRhs(masterEnv);
                   double newCutValue = 0;
                   double newViolation = 0;
                   double totvalue = 1;

                   set<NODE> cutset;
                   for (auto s : root_adj_nodes) {
                       if (reached[s] && reached[t] &&
                           forest.find_set(s) == forest.find_set(t)) {
                           pair_i_k.second = k;
                           pair_i_k.first = s;
                           newCutLhs += (partition_node_vars.at(pair_i_k));
                           newCutValue += xSol.at(pair_i_k);  // 0

                           cutset.insert(s);

                       } else
                           continue;
                   }

                   IloNumVar temp_var =
                       IloNumVar(masterEnv, 1, 1, IloNumVar::Float);
                   newCutRhs += temp_var;

                   newViolation = 1.0 - newCutValue;

                   if (newCutValue < 1 - TOL && cutset.size() != 0) {
                       cutLhs.push_back(newCutLhs);
                       cutRhs.push_back(newCutRhs);
                       violation.push_back(newViolation);

                       LOG << "lhs: " << newCutLhs << endl;
                       LOG << "rhs: " << newCutRhs << endl;
                       LOG << "violation: " << newViolation << endl;

                       if (newViolation >= TOL) ret = true;
                   }

                   cutpool.AddLhs(k, cutset);
                   cutpool.AddViolation(k, newViolation);
                   cutset.clear();
               }
           }
       }
        */

    ret = cutLhs.size() > 0 ? 1 : 0;
    return ret;
}

/*  Min cut seperation for NS  */
bool seperate_min_cut_ns(
    IloEnv masterEnv, const map<pair<NODE, INDEX>, double>& xSol,
    std::shared_ptr<Graph> G,
    const map<pair<NODE, INDEX>, IloNumVar>& partition_node_vars,
    vector<IloExpr>& cutLhs, vector<IloExpr>& cutRhs, vector<double>& violation,
    const map<INDEX, NODE>& ns_root) {
    // cout << " ------------------- use seperate_min_cut_ns// " << endl;
    bool ret = false;
    pair<NODE, INDEX> pair_i_k;
    cutLhs = vector<IloExpr>();
    cutRhs = vector<IloExpr>();
    violation = vector<double>();
    double a, b;

    for (auto k : G->p_set()) {
        std::shared_ptr<SUB_Graph> subG =
            std::make_shared<SUB_Graph>(G->get_subgraph()[k]);

        GNsize = G->nodes().size() + 1;
        N = G->nodes().size() + 1;
        M = subG->arcs().size() + 1;
        T = N + M + 1;
        dis = (int*)malloc(int(2 * GNsize) * sizeof(int));
        work = (int*)malloc(int(2 * GNsize) * sizeof(int));
        use = (int*)malloc(int(2 * GNsize) * sizeof(int));
        q = (int*)malloc(int(2 * T) * sizeof(int));
        DinicEdge = (st*)malloc(int(2 * T) * sizeof(st));

        bool flag = BuildDinicGraph(G, subG, xSol, k, edge[k], ns_root.at(k));
        if (!flag) {
            // SingleLazySeperate(masterEnv, subG, xSol, k,
            // partition_node_vars,cutLhs, cutRhs, violation, ns_root);
            continue;
        }

        /*for (auto i : subG->nodes()) {
                pair_i_k = NODE_PAIR(i, k);
                cout << i << " " << xSol.at(pair_i_k) << endl;
        }*/

        /*int idx = 0;
        for (auto u : subG->nodes()) {
                for (int i = head[k][u]; i != -1; i = edge[k][i].next) {
                        cout << u << " to " << " " << edge[k][i].v << " is " <<
        edge[k][i].w
                                << endl;
                        idx++;
                }
                for (int i = head[k][u + GNsize-1]; i != -1; i =
        edge[k][i].next) { cout << u + GNsize-1 << " to " << " " << edge[k][i].v
        << " is " << edge[k][i].w
                                << endl;
                        idx++;
                }
        }
        cout << idx << endl << endl;*/

        IloExpr newCutLhs;
        IloExpr newCutRhs;
        double newViolation;
        double min_cut_value;

        int root = ns_root.at(k);
        for (auto _q : subG->t_set()) {
            if (_q == root) continue;
            memcpy(DinicEdge, edge[k], sizeof(st) * int(2 * T));
            /* rep(i, 1, T) cout << edge[k][i].v << " " << edge[k][i].next <<
               "->"
                              << DinicEdge[i].v << DinicEdge[i].next << endl; */
            min_cut_value = Dinic(root, _q, head[k], DinicEdge);

            // cout << "Maxflow for Dinic " << _q << " is " << min_cut_value<<
            // endl;

            a = 0.0;

            if (min_cut_value < 1 - EPS) {
                newCutLhs = IloExpr(masterEnv);
                newCutRhs = IloExpr(masterEnv);
                newViolation = 1 - min_cut_value;
                set<NODE> cutset;

                memset(use, 0, sizeof(int) * (2 * GNsize));

                FindCut(root, head[k], DinicEdge);

                int n = G->nodes().size();
                for (auto i : subG->nodes()) {
                    if (i == root) continue;
                    if ((use[i] && !use[i + n]) || (use[i + n] && !use[i])) {
                        pair_i_k = NODE_PAIR(i, k);
                        newCutLhs += (partition_node_vars.at(pair_i_k));
                        cutset.insert(i);

                        // a += xSol.at(pair_i_k);
                    }
                }

                IloNumVar temp_var = IloNumVar(masterEnv, 1, 1, IloNumVar::Int);
                newCutRhs += temp_var;

                cutLhs.push_back(newCutLhs);
                cutRhs.push_back(newCutRhs);
                violation.push_back(newViolation);

                cutpool.AddLhs(k, cutset);
                cutpool.AddViolation(k, newViolation);
                // cout << "cut value is: " << a << endl;
                // cout << "Partition " << k << endl;
                // cout << cutset;
                cutset.clear();
            }
        }

        free(dis);
        free(work);
        free(use);
        free(q);
        free(DinicEdge);
    }

    // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    const double INF = 10000000000000.0;
    /*
        for (auto k : G->p_set()) {
            LOG << "For partation: " << k << endl;
            pair_i_k.second = k;

            ListDigraph::ArcMap<double> x_capacities(ns_mincut_capgraph[k],
       INF); SUB_Graph subG = G->get_subgraph()[k];

            build_cap_graph_ns(G, x_capacities, k, ns_root, xSol);

            LOG << "Built cap graph..." << endl;

            IloExpr newCutLhs;
            IloExpr newCutRhs;
            double newViolation;
            double min_cut_value;

            for (auto q : subG.t_set()) {
                if (q == ns_root.at(k)) continue;

                Preflow<ListDigraph, ListDigraph::ArcMap<double>> min_cut(
                    ns_mincut_capgraph[k], x_capacities,
                    ns_mincut_v_nodes[k][ns_root.at(k)].first,
                    ns_mincut_v_nodes[k][q].first);
                min_cut.runMinCut();
                min_cut_value = min_cut.flowValue();
                b = min_cut_value;

                cout << "Min-cut for Lemon " << q << " is " << min_cut_value
                     << endl;

                if (min_cut_value < 1 - TOL) {
                    newCutLhs = IloExpr(masterEnv);
                    newCutRhs = IloExpr(masterEnv);
                    newViolation = 1 - min_cut_value;
                    set<NODE> cutset;

                    for (NODE i : subG.nodes()) {
                        if (i == ns_root.at(k)) continue;
                        ListNode a = ns_mincut_v_nodes[k][i].first;
                        ListNode b = ns_mincut_v_nodes[k][i].second;

                        if ((min_cut.minCut(a) && !min_cut.minCut(b)) ||
                            (!min_cut.minCut(a) && min_cut.minCut(b))) {
                            LOG << "find cur arc: " << i << endl;
                            pair_i_k.first = i;
                            newCutLhs += (partition_node_vars.at(pair_i_k));

                            cutset.insert(i);
                        }
                    }
                    IloNumVar temp_var = IloNumVar(masterEnv, 1, 1,
       IloNumVar::Int); newCutRhs += temp_var;

                    // cutLhs.push_back(newCutLhs);
                    // cutRhs.push_back(newCutRhs);
                    // violation.push_back(newViolation);

                    cutpool.AddLhs(k, cutset);
                    cutpool.AddViolation(k, newViolation);
                    cout << "Partition " << k << endl;
                    cout << cutset;
                    cutset.clear();

                    LOG << "node " << q << endl;
                    LOG << "cut " << cutLhs.size() << endl;
                    LOG << "flowValue " << min_cut_value << endl;
                    LOG << "lhs: " << newCutLhs << endl;
                    LOG << "rhs: " << newCutRhs << endl;
                }
            }
        }
     */
    ret = cutLhs.size() > 0 ? 1 : 0;
    return ret;
}

bool seperate_min_cut_ns_cutpool(
    IloEnv masterEnv, const map<pair<NODE, INDEX>, double>& xSol,
    std::shared_ptr<Graph> G,
    const map<pair<NODE, INDEX>, IloNumVar>& partition_node_vars,
    vector<IloExpr>& cutLhs, vector<IloExpr>& cutRhs, vector<double>& violation,
    const map<INDEX, NODE>& ns_root) {
    // cout << " ------------------------- use
    // seperate_min_cut_ns_cutpool----------------" << endl;
    const double INF = 10000000000000.0;
    bool ret = false;
    pair<NODE, INDEX> pair_i_k;
    map<INDEX, NODE_SET> V_k_set = G->v_set();
    map<INDEX, NODE_SET> T_k_set = G->t_set();

    cutLhs = vector<IloExpr>();
    cutRhs = vector<IloExpr>();
    violation = vector<double>();
    for (auto k : G->p_set()) {
        SUB_Graph subG = G->get_subgraph()[k];
        vector<NODE_SET> Kth_CutPool = cutpool.cutPoolLhs()[k];
        for (auto s : Kth_CutPool) {
            // calculate sth-constraint value
            double cutValue = 0.0;
            for (auto i : s) {
                pair_i_k.first = i;
                pair_i_k.second = k;

                cutValue += xSol.at(pair_i_k);
            }
            if (cutValue < 1.0 - TOL) {
                IloExpr newCutLhs = IloExpr(masterEnv);
                IloExpr newCutRhs = IloExpr(masterEnv);
                double newViolation = 1 - cutValue;

                for (auto i : s) {
                    pair_i_k.first = i;
                    pair_i_k.second = k;

                    newCutLhs += partition_node_vars.at(pair_i_k);
                }
                IloNumVar temp_var = IloNumVar(masterEnv, 1, 1, IloNumVar::Int);
                newCutRhs += temp_var;

                cutLhs.push_back(newCutLhs);
                cutRhs.push_back(newCutRhs);
                violation.push_back(newViolation);
            }
        }
    }

    ret = cutLhs.size() > 0 ? 1 : 0;
    return ret;
}