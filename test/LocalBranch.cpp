#include "LocalBranch.h"

#include <stdlib.h>

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

#include "callback.h"
#include "separation.h"
#include "smp.h"
#include "type.h"

using namespace std;

LBSolver::LBSolver(IloEnv env, std::shared_ptr<Graph> g_ptr,
                   SmpForm _formulation, int _callbackOption, bool _relax,
                   bool _ns_sep_opt, int _LB_CallTime, int _LB_InterationTime,
                   int _R, int _BCSolNum, int _BCTime, double _epsilon_lazy,
                   double _epsilon_user, int _max_cuts_lazy, int _max_cuts_user,
                   IloNumVarArray _x_vararray,
                   IloNumVarArray _x_vararray_primal,
                   map<NODE, IloNumVar> _primal_node_vars,
                   map<pair<NODE, INDEX>, IloNumVar> _partition_node_vars,
                   map<INDEX, NODE> _ns_root,
                   map<pair<NODE, INDEX>, int> _x_varindex_ns,
                   map<NODE, int> _x_varindex_ns_primal) {
    model = IloModel(env);
    objective = IloObjective();

    formulation = _formulation;
    G = g_ptr;
    callbackOption = _callbackOption;
    relax = _relax;
    ns_sep_opt = _ns_sep_opt;
    max_cuts_lazy = _max_cuts_lazy;
    max_cuts_user = _max_cuts_user;
    tol_lazy = _epsilon_lazy;
    tol_user = _epsilon_user;

    LB_CallTime = _LB_CallTime;
    LB_InterationTime = _LB_InterationTime;
    R = _R;
    BCSolNum = _BCSolNum;
    BCTime = _BCTime;

    x_vararray = _x_vararray;
    x_vararray_primal = _x_vararray_primal;
    primal_node_vars = _primal_node_vars;
    partition_node_vars = _partition_node_vars;
    ns_root = _ns_root;
    x_varindex_ns = _x_varindex_ns;
    x_varindex_ns_primal = _x_varindex_ns_primal;

    cplex = IloCplex(model);
    cplex.setParam(IloCplex::IntSolLim, BCSolNum);
    cplex.setParam(IloCplex::TiLim, BCTime);

    switch (formulation) {
        case NONE:
        case SCF:
        case MCF:
        case STEINER:
            break;

        case NS: {
            switch (callbackOption) {
                case 0:
                    break;
                case 1:
                    cplex.use(NS_StrongComponentLazyCallback(
                        env, G, partition_node_vars, x_vararray, x_varindex_ns,
                        tol_lazy, max_cuts_lazy, formulation, ns_root,
                        x_vararray_primal, x_varindex_ns_primal));
                    break;
                case 2:
                    cplex.use(NS_CutCallback(env, G, partition_node_vars,
                                             x_vararray, x_varindex_ns,
                                             tol_user, max_cuts_user,
                                             formulation, ns_root, ns_sep_opt));
                    break;
                case 3:
                    cplex.use(NS_StrongComponentLazyCallback(
                        env, G, partition_node_vars, x_vararray, x_varindex_ns,
                        tol_lazy, max_cuts_lazy, formulation, ns_root,
                        x_vararray_primal, x_varindex_ns_primal));
                    cplex.use(NS_CutCallback(env, G, partition_node_vars,
                                             x_vararray, x_varindex_ns,
                                             tol_user, max_cuts_user,
                                             formulation, ns_root, ns_sep_opt));
                    break;
                default:
                    break;
            }
            break;
        }
        default:
            break;
    }
}

void LBSolver::Floyd(map<NODE, vector<int>>& subGnodesIdx,
                     map<int, NODE>& rev_subGnodesIdx,
                     vector<vector<int>>& distance, vector<vector<int>>& path,
                     int& idx, int k) {
    const int NodesValueINF = 1000 * G->nodes().size() + 1;
    SUB_Graph subG = G->get_subgraph()[k];

    // Add nodes
    idx = 0;  // Total Nodes Number
    for (auto i : subG.nodes()) {
        subGnodesIdx[i].resize(2);

        subGnodesIdx[i][0] = ++idx;
        rev_subGnodesIdx[idx] = i;

        subGnodesIdx[i][1] = ++idx;
        rev_subGnodesIdx[idx] = i;
    }
    distance.resize(idx + 1);
    path.resize(idx + 1);
    for (int i = 1; i <= idx; i++) {
        distance[i].resize(idx, NodesValueINF);
        path[i].resize(idx, -1);
    }

    // Add arc
    for (auto i : subG.nodes()) {
        int u = subGnodesIdx[i][0];
        int v = subGnodesIdx[i][1];
        distance[u][v] = (int)(subG.node_value().at(i));
        distance[v][u] = (int)(subG.node_value().at(i));
    }
    for (auto arc : subG.arcs()) {
        int u = subGnodesIdx[arc.first][1];
        int v = subGnodesIdx[arc.second][0];
        distance[u][v] = 0;
        distance[v][u] = 0;
    }

    // Run Floyd
    for (int K = 1; K <= idx; K++) {
        for (int i = 1; i <= idx; i++) {
            for (int j = 1; j <= idx; j++) {
                if (distance[i][j] < NodesValueINF &&
                    distance[K][j] < NodesValueINF &&
                    distance[i][j] > distance[i][K] + distance[K][j]) {
                    distance[i][j] = distance[i][K] + distance[K][j];
                    path[i][j] = k;
                }
            }
        }
    }

    return;
}

void LBSolver::GenerateInitialSolution(int k, map<NODE, bool>& xPartSol) {
    SUB_Graph subG = G->get_subgraph()[k];
    map<NODE, vector<int>> subGnodesIdx;
    map<int, NODE> rev_subGnodesIdx;
    vector<vector<int>> distance;
    vector<vector<int>> path;
    int idx = 0;

    Floyd(subGnodesIdx, rev_subGnodesIdx, distance, path, idx, k);

    for (auto i : subG.nodes()) {
        xPartSol[i] = false;
    }

    auto t1 = subG.t_set().begin();
    auto t2 = subG.t_set().begin();
    t2++;
    while (t2 != subG.t_set().end()) {
        int s = *t1, t = *t2;
        xPartSol[s] = 1;
        xPartSol[t] = 1;
        while (path[s][t] != -1) {
            xPartSol[path[s][t]] = 1;
            s = path[s][t];
        }
        t1++, t2++;
    }

    return;
}

void LBSolver::LocalBranchSearch() {
    for (int lbtime = 1; lbtime < LB_CallTime; lbtime++) {
        for (auto k : G->p_set()) {
            SUB_Graph subG = G->get_subgraph()[k];
            map<NODE, bool> xPartSol;
            GenerateInitialSolution(k, xPartSol);
        }
    }
    return;
}

void LBSolver::LocalBranch() {}