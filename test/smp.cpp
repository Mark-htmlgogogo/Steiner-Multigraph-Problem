#include "smp.h"

#include <stdlib.h>

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

#include "callback.h"
#define SPACING 9
#define LSPACING 16
#define LOG \
    if (false) cerr
#define TOL 0.001

using namespace std;

CutPool cutpool;

/*******************************************************/
/************* Build the model *************************/
/*******************************************************/

/* Constructor */

SmpSolver::SmpSolver(IloEnv env, std::shared_ptr<Graph> g_ptr,
                     SmpForm formulation_, double epsilon_lazy_,
                     double epsilon_user_, int time_limit_, int max_cuts_lazy_,
                     int max_cuts_user_, int callbackOption_, bool relax_,
                     bool ns_sep_opt_, string filename_, int LB_CP_Option_,
                     int lazy_sep_opt_, int MIPDisplayLevel_) {
    /* Initialize Cplex Sturctures */
    model = IloModel(env);
    objective = IloObjective();

    /* settings */
    formulation = formulation_;
    G = g_ptr;
    time_limit = time_limit_;
    max_cuts_lazy = max_cuts_lazy_;
    max_cuts_user = max_cuts_user_;
    tol_lazy = epsilon_lazy_;
    tol_user = epsilon_user_;
    callbackOption = callbackOption_;
    relax = relax_;
    ns_sep_opt = ns_sep_opt_;
    filename = filename_;
    LB_CP_Option = LB_CP_Option_;
    fianlsolveflag = 0;
    lazy_sep_opt = lazy_sep_opt_;
    MIPDisplayLevel = MIPDisplayLevel_;

    /* Add x_i variables: primal_node_vars */
    char var_name[255];
    set<NODE> T = G->t_total();
    int idx = 0;
    x_vararray_primal = IloNumVarArray(env);
    for (auto& node : G->nodes()) {
        IloNumVar var;
        snprintf(var_name, 255, "x_%d", node);
        // set x_i = 1 if i belongs to T, others to {0,1}
        if (T.find(node) != T.end())  // Constraint(2)
        {
            if (relax)
                var = IloNumVar(env, 1, 1, IloNumVar::Float, var_name);
            else
                var = IloNumVar(env, 1, 1, IloNumVar::Bool, var_name);
        } else {
            if (relax)
                var = IloNumVar(env, 0, 1, IloNumVar::Float, var_name);
            else
                var = IloNumVar(env, 0, 1, IloNumVar::Bool, var_name);
        }
        primal_node_vars[node] = var;
        x_vararray_primal.add(var);
        x_varindex_ns_primal[node] = idx++;
        // printInfo(var);
    }

    LOG << "Building problem...";
    cout << "Begin to build problem..." << endl;

    switch (formulation) {
        case SCF:
            build_problem_scf();
            break;
        case MCF:
            build_problem_mcf_terminal();
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
    try {
        cplex = IloCplex(model);
    } catch (IloException e) {
        cout << e << endl;
    }
    // cplex = IloCplex(model);  // create a ILOG CPLEX algorithm and extract a
    // model
    cplex.setParam(IloCplex::MIPDisplay, MIPDisplayLevel);  // set display level
    if (formulation > 0) cplex.setParam(IloCplex::AdvInd, 1);  // start value: 1
    cplex.setParam(IloCplex::EpGap, 1e-09);  // set MIP gap tolerance
    cplex.setParam(IloCplex::Threads, 8);  // set the number of parallel threads
    cplex.setParam(IloCplex::TreLim,
                   12288);  // set the limit of tree memory in megabytes
    cplex.setParam(IloCplex::TiLim,
                   time_limit);  // set time limit in secs, default:1200
    cplex.setParam(IloCplex::EpInt, 1e-06);  // set integrality tolerance

    /* Implement user cut or lazy constraint here*/
    switch (formulation) {
        case NONE:
        case SCF:
        case MCF:
            break;
        case STEINER: {
            switch (callbackOption) {
                case 0:
                    break;
                case 1:
                    cplex.use(StrongComponentLazyCallback(
                        env, G, edge_vars, x_vararray, x_varindex_steiner,
                        tol_lazy, max_cuts_lazy, formulation, Steiner_root,
                        primal_node_vars));
                    break;
                case 2:
                    cplex.use(SmpCutCallback(env, G, edge_vars, x_vararray,
                                             x_varindex_steiner, tol_user,
                                             max_cuts_user, formulation,
                                             Steiner_root, primal_node_vars));
                    break;
                case 3:
                    cplex.use(StrongComponentLazyCallback(
                        env, G, edge_vars, x_vararray, x_varindex_steiner,
                        tol_lazy, max_cuts_lazy, formulation, Steiner_root,
                        primal_node_vars));
                    cplex.use(SmpCutCallback(env, G, edge_vars, x_vararray,
                                             x_varindex_steiner, tol_user,
                                             max_cuts_user, formulation,
                                             Steiner_root, primal_node_vars));
                    break;
                default:
                    break;
            }
            break;
        }
        case NS: {
            switch (callbackOption) {
                case 0:
                    break;
                case 1:
                    cplex.use(NS_StrongComponentLazyCallback(
                        env, G, partition_node_vars, x_vararray, x_varindex_ns,
                        tol_lazy, max_cuts_lazy, formulation, ns_root,
                        x_vararray_primal, x_varindex_ns_primal, lazy_sep_opt));
                    break;
                case 2:
                    cplex.use(NS_CutCallback(
                        env, G, partition_node_vars, x_vararray, x_varindex_ns,
                        tol_user, max_cuts_user, formulation, ns_root,
                        ns_sep_opt, LB_CP_Option, fianlsolveflag,
                        lazy_sep_opt));
                    break;
                case 3:
                    cplex.use(NS_StrongComponentLazyCallback(
                        env, G, partition_node_vars, x_vararray, x_varindex_ns,
                        tol_lazy, max_cuts_lazy, formulation, ns_root,
                        x_vararray_primal, x_varindex_ns_primal, lazy_sep_opt));
                    cplex.use(NS_CutCallback(
                        env, G, partition_node_vars, x_vararray, x_varindex_ns,
                        tol_user, max_cuts_user, formulation, ns_root,
                        ns_sep_opt, LB_CP_Option, fianlsolveflag,
                        lazy_sep_opt));
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

/* Update the objective function */
void SmpSolver::update_problem(const const map<NODE, double>& obj_coeff,
                               SmpForm formulation) {
    /* Build objective function */
    IloEnv env = model.getEnv();
    IloExpr totalCost(env);
    double objcoeff = 0.0;
    IloNumVar var;

    for (auto& node : G->nodes()) {
        var = primal_node_vars[node];
        objcoeff = obj_coeff.at(node);
        totalCost += var * objcoeff;
    }

    // change the objective of NS: \sigma\sigma x_i^k
    /* if (formulation == NS) {
        double M = 0;
        IloExpr sigma_vars(env);
        for (auto k : G->p_set()) {
            SUB_Graph subG = G->get_subgraph()[k];
            pair<NODE, INDEX> pair_i_k;
            for (auto i : subG.nodes()) {
                pair_i_k.first = i;
                pair_i_k.second = k;
                sigma_vars += partition_node_vars[pair_i_k];
                M += 1.0;
            }
        }
        totalCost += ((1.0 / (M + 1.0)) * sigma_vars);
    } */

    objective = IloObjective(env, totalCost, IloObjective::Minimize);
    model.add(objective);
}

void SmpSolver::solve() {
    LOG << "--------------------- SOLVING SMP " << endl;

    double start_time = cplex.getCplexTime();
    double start_ticks = cplex.getDetTime();

    // cout << "Number of constraints: " << cplex.getNrows() << endl;
    // cout << "Number of variables(int):   " << cplex.getNintVars() << endl;
    // cout << "Number of variables(binary):   " << cplex.getNbinVars() << endl;
    // cout << "Number of variables(cols):   " << cplex.getNcols() << endl;

    try {
        cplex.solve();
    } catch (IloException e) {
        cout << e << endl;
    }

    elapsed_time = cplex.getCplexTime() - start_time;
    elapsed_ticks = cplex.getDetTime() - start_ticks;
    cout << "Solution status \t= \t" << cplex.getStatus() << endl;

    // cout << "Elapsed ticks \t= \t" << elapsed_ticks << endl;
    cout << "Gap \t= \t" << cplex.getMIPRelativeGap() << endl;
    cout << "Elapsed time \t= \t" << elapsed_time << endl;

    cout << "Objectvie value \t= \t" << cplex.getObjValue() << endl;
    cout << "Number of nodes \t= \t" << cplex.getNnodes() << endl;
    cout << "Number of cuts \t= \t" << cplex.getNcuts(IloCplex::CutUser)
         << endl;

    /*for (auto var : primal_node_vars)
            cout << var.second.getName() << "\t" << cplex.getValue(var.second)
            << endl;*/
    /*for (auto var : source_node_vars)
            cout << var.second.getName() << "\t" << cplex.getValue(var.second)
            << endl;*/
    /*for (auto var : partition_node_vars)
            cout << var.second.getName() << "\t" << cplex.getValue(var.second)
            << endl;*/
    /*for (auto var : partition_flow_vars)
            cout << var.second.getName() << "\t" << cplex.getValue(var.second)
            << endl;
    for (auto var : multi_flow_vars)
            cout << var.second.getName() << "\t" << cplex.getValue(var.second)
            << endl;
    for (auto var : edge_vars)
            cout << var.second.getName() << "\t" << cplex.getValue(var.second)
            << endl;*/

    print_to_file();
}

void SmpSolver::clear() {
    LOG << "Clearing cuts and lazy constraints" << endl;
    cplex.clearUserCuts();
    cplex.clearLazyConstraints();
}

/* Single commodity flow formulation */
void SmpSolver::build_problem_scf() {
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
    for (auto k : G->p_set()) {
        // cout <<  "k = " << k << endl;
        V_k_size = static_cast<int>(V_k_set[k].size());
        subG = G->get_subgraph()[k];

        // Add  x_i^k (binary), for each node in G[V_k]:
        pair_i_k.second = k;
        for (auto i : subG.nodes()) {
            IloNumVar var;
            snprintf(var_name, 255, "x_%d^%d", i, k);
            if (T_k_set[k].find(i) != T_k_set[k].end()) {
                if (relax)
                    var = IloNumVar(env, 1, 1, IloNumVar::Float, var_name);
                else
                    var = IloNumVar(env, 1, 1, IloNumVar::Int, var_name);
            }

            else {
                if (relax)
                    var = IloNumVar(env, 0, 1, IloNumVar::Float, var_name);
                else
                    var = IloNumVar(env, 0, 1, IloNumVar::Int, var_name);
            }
            pair_i_k.first = i;
            partition_node_vars[pair_i_k] = var;
            // printInfo(var);
        }

        // Add x_s_k (float [0, |V_k_set|]), for each k in P_index_set:
        snprintf(var_name, 255, "x_s%d", k);
        temp_var = IloNumVar(env, 0, V_k_size, IloNumVar::Float, var_name);
        source_node_vars[k] = temp_var;
        // printInfo(temp_var);

        // for each T_k, choose a root r_k
        auto firstElement = T_k_set.at(k).begin();
        root[k] = *firstElement;
        // cout << "r" << k << " = " << root[k] << endl;

        // Add y_ij^k (float), for arc (s_k, r_k) and each arc in G[V_k]:
        pair_ij_k.second = k;
        pair_source_root.first = 0;
        pair_source_root.second = root[k];
        pair_ij_k.first = pair_source_root;
        // Add y_(s_k,r_k)
        snprintf(var_name, 255, "y_s%dr%d", k, k);
        temp_var = IloNumVar(env, 0.0, IloInfinity, IloNumVar::Float, var_name);
        partition_flow_vars[pair_ij_k] = temp_var;
        // printInfo(temp_var);

        // Add  y_ij^k
        for (auto& arc : subG.arcs()) {
            IloNumVar var;
            snprintf(var_name, 255, "y_%d%d^%d", arc.first, arc.second, k);
            var = IloNumVar(env, 0.0, IloInfinity, IloNumVar::Float, var_name);
            pair_ij_k.first = arc;
            partition_flow_vars[pair_ij_k] = var;
            // printInfo(var);
        }
    }

    /*******************/
    /* Add constraints */
    /*******************/

    cout << "Begin to add constraints..." << endl;

    // Constraint��5��: x_i >= x_i^k k \in P, i \in v_k\t_k
    for (auto i : G->v_total()) {
        // cout << "constraint(5)" << endl;
        if (std::find(G->t_total().begin(), G->t_total().end(), i) !=
            G->t_total().end())
            continue;
        for (auto k : G->nodes_of_v().at(i)) {
            pair_i_k.first = i;
            pair_i_k.second = k;
            snprintf(con_name, 255, "%s >= %s (5)",
                     primal_node_vars[i].getName(),
                     partition_node_vars[pair_i_k].getName());
            model.add(primal_node_vars[i] >= partition_node_vars[pair_i_k])
                .setName(con_name);
            // cout << con_name << endl;
        }
    }

    // Constraints (6) - (9)
    for (auto k : G->p_set()) {
        // Set variables's k part
        V_k_size = static_cast<int>(V_k_set.at(k).size());
        subG = G->get_subgraph()[k];
        pair_i_k.second = k;
        pair_ij_k.second = k;
        pair_source_root.second = root[k];

        // Constraint (6): flow conservation at s_k
        //  y_(s_k, r_k) + x_s_k = |V_k|
        // cout << "constraint(6)" << endl;
        pair_source_root.first = 0;
        pair_source_root.second = root[k];
        pair_ij_k.first = pair_source_root;

        snprintf(con_name, 255, "%s + %s >= %d (6)",
                 partition_flow_vars[pair_ij_k].getName(),
                 source_node_vars[k].getName(), V_k_size);
        model
            .add(partition_flow_vars[pair_ij_k] + source_node_vars[k] >=
                 V_k_size)
            .setName(con_name);
        // cout << con_name << endl;

        snprintf(con_name, 255, "%s + %s <= %d (6)",
                 partition_flow_vars[pair_ij_k].getName(),
                 source_node_vars[k].getName(), V_k_size);
        model
            .add(partition_flow_vars[pair_ij_k] + source_node_vars[k] <=
                 V_k_size)
            .setName(con_name);
        // cout << con_name << endl;

        // Constraint (7): y_ij^k <= |V_k| * x_j^k
        // cout << "constraint(7)" << endl;
        for (auto& arc : subG.arcs()) {
            pair_ij_k.first = arc;
            pair_i_k.first = arc.second;

            snprintf(con_name, 255, "%s <= %d * %s (7)",
                     partition_flow_vars[pair_ij_k].getName(), V_k_size,
                     partition_node_vars[pair_i_k].getName());
            model
                .add(partition_flow_vars[pair_ij_k] <=
                     V_k_size * partition_node_vars[pair_i_k])
                .setName(con_name);
            cout << con_name << endl;
        }

        // Constraint (8): flow conservation at x_j^k
        // (k-flow)->j = x_j_k + y->(k-flow)

        // cout << "constraint(8)" << endl;
        for (auto j : subG.nodes()) {
            IloExpr flow_in_j(env);
            IloExpr flow_out_j(env);
            pair_ij_k.first.second = j;
            char flow_in_j_name[255];
            memset(flow_in_j_name, 0, 255);
            for (auto i : subG.adj_nodes_list().at(j)) {
                // cout << "in" << j << ": " << i << endl;
                pair_ij_k.first.first = i;
                flow_in_j += partition_flow_vars[pair_ij_k];
                strcat(flow_in_j_name,
                       partition_flow_vars[pair_ij_k].getName());
                strcat(flow_in_j_name, " + ");
            }
            // consider r_k and (s_k, r_k)
            if (j == root[k]) {
                // cout << "in" << j << ": "
                //<< "s" << k << endl;
                pair_ij_k.first = pair_source_root;
                flow_in_j += partition_flow_vars[pair_ij_k];
                strcat(flow_in_j_name,
                       partition_flow_vars[pair_ij_k].getName());
            }
            flow_in_j.setName(flow_in_j_name);

            pair_ij_k.first.first = j;
            char flow_out_j_name[255];
            memset(flow_out_j_name, 0, 255);
            for (auto i : subG.adj_nodes_list().at(j)) {
                // cout << "out" << j << ": " << i << endl;
                pair_ij_k.first.second = i;
                flow_out_j += partition_flow_vars[pair_ij_k];
                strcat(flow_out_j_name,
                       partition_flow_vars[pair_ij_k].getName());
                strcat(flow_out_j_name, " + ");
            }
            flow_out_j.setName(flow_out_j_name);

            pair_i_k.first = j;

            snprintf(con_name, 255, "%s >= %s + %s (8)", flow_in_j.getName(),
                     flow_out_j.getName(),
                     partition_node_vars[pair_i_k].getName());
            model.add(flow_in_j >= flow_out_j + partition_node_vars[pair_i_k])
                .setName(con_name);
            // cout << con_name << endl;

            snprintf(con_name, 255, "%s <= %s + %s (8)", flow_in_j.getName(),
                     flow_out_j.getName(),
                     partition_node_vars[pair_i_k].getName());
            model.add(flow_in_j <= flow_out_j + partition_node_vars[pair_i_k])
                .setName(con_name);
            // cout << con_name << endl;
        }

        // Constraint (9): total system flow = flow of (0, r_k)
        // sum(x_j^k) = y_(s_k,r_k)
        // cout << "constraint(9)" << endl;
        IloExpr selected_node_size(env);
        for (auto i : subG.nodes()) {
            pair_i_k.first = i;
            selected_node_size += partition_node_vars[pair_i_k];
            // cout << partition_node_vars[pair_i_k].getName() << " + ";
        }
        pair_ij_k.first = pair_source_root;

        // cout << " = " << partition_flow_vars[pair_ij_k].getName() << endl;

        snprintf(con_name, 255, "selected_node_size_%d", k);
        selected_node_size.setName(con_name);

        snprintf(con_name, 255, "%s >= %s (9)", selected_node_size.getName(),
                 partition_flow_vars[pair_ij_k].getName());
        model.add(selected_node_size >= partition_flow_vars[pair_ij_k])
            .setName(con_name);
        // cout << con_name << endl;

        snprintf(con_name, 255, "%s <= %s (9)", selected_node_size.getName(),
                 partition_flow_vars[pair_ij_k].getName());
        model.add(selected_node_size <= partition_flow_vars[pair_ij_k])
            .setName(con_name);
        // cout << con_name << endl;
    }
}

/* Multi-commodity flow formulation */
void SmpSolver::build_problem_mcf() {
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
    for (auto k : G->p_set()) {
        V_k_size = static_cast<int>(V_k_set[k].size());
        subG = G->get_subgraph()[k];
        int subG_nodesNum = subG.nodes().size();

        // Add  x_i^k (binary), for each node in G[V_k]:
        pair_i_k.second = k;
        for (auto i : subG.nodes()) {
            IloNumVar var;
            if (subG.CheckNodeIsTerminal().at(i)) {
                // if (T_k_set[k].find(i) != T_k_set[k].end()) {
                if (relax)
                    var = IloNumVar(env, 1, 1, IloNumVar::Float, var_name);
                else
                    var = IloNumVar(env, 1, 1, IloNumVar::Int, var_name);
            } else {
                if (relax)
                    var = IloNumVar(env, 0, 1, IloNumVar::Float, var_name);
                else
                    var = IloNumVar(env, 0, 1, IloNumVar::Int, var_name);
            }
            pair_i_k.first = i;
            partition_node_vars[pair_i_k] = var;
            // printInfo(var);
        }

        // For each T_k, choose a root r_k
        auto firstElement = T_k_set[k].begin();
        root[k] = *firstElement;
        // cout << "r" << k << " = " << root[k] << endl;

        // Add y_ij_km (float)
        pair_ij_km.second.first = k;
        for (auto m : subG.nodes()) {
            if (m == root[k]) continue;
            pair_ij_km.second.second = m;
            for (auto arc : subG.arcs()) {
                IloNumVar var;
                var = IloNumVar(env, 0.0, IloInfinity, IloNumVar::Float,
                                var_name);
                pair_ij_km.first = arc;
                multi_flow_vars[pair_ij_km] = var;
                // printInfo(var);
            }
        }
    }

    /*******************/
    /* Add constraints */
    /*******************/

    cout << "Begin to add constraints..." << endl;

    // Add constraint (13-1)
    for (auto i : G->v_total()) {
        // cout << "constraint(13-1)" << endl;
        if (std::find(G->t_total().begin(), G->t_total().end(), i) !=
            G->t_total().end())
            continue;
        for (auto k : G->nodes_of_v().at(i)) {
            pair_i_k.first = i;
            pair_i_k.second = k;
            model.add(primal_node_vars[i] >= partition_node_vars[pair_i_k]);
            // cout << con_name << endl;
        }
    }

    // Add constraint (14-21)
    // cout << "Add constraint (14-16)" << endl;
    for (auto k : G->p_set()) {
        // cout << endl
        //<< "For partition: " << k << endl;
        V_k_size = static_cast<int>(V_k_set[k].size());
        subG = G->get_subgraph()[k];
        pair_ij_km.second.first = k;

        // Add constraint (15-18)
        for (auto m : subG.nodes()) {
            if (m == root[k]) continue;
            pair_ij_km.second.second = m;

            // Partition for con_15： y_ir^k_km
            IloRange con_15 = IloRange(env, 0, 0);
            for (auto i : subG.adj_nodes_list().at(root[k])) {
                pair_ij_km.first.first = i;
                pair_ij_km.first.second = root[k];
                con_15.setLinearCoef(multi_flow_vars[pair_ij_km], 1);
            }
            model.add(con_15);
            // cout << "con_15_name: " << con_15_name << endl;

            // Partition for con_16 ： y_im_km
            pair_i_k.first = m;
            pair_i_k.second = k;
            IloNumVar coeff = partition_node_vars[pair_i_k];
            IloExpr flow_in_m(env);
            for (auto i : subG.adj_nodes_list().at(m)) {
                pair_ij_km.first.first = i;
                pair_ij_km.first.second = m;
                flow_in_m += multi_flow_vars[pair_ij_km];
            }
            // cout << con_name << endl;
            model.add(flow_in_m <= coeff).setName(con_name);

            // cout << con_name << endl;
            model.add(flow_in_m >= coeff).setName(con_name);

            // Partition for con_17： y_mj_km
            IloRange con_17 = IloRange(env, 0, 0);
            for (auto j : subG.adj_nodes_list().at(m)) {
                pair_ij_km.first.first = m;
                pair_ij_km.first.second = j;
                con_17.setLinearCoef(multi_flow_vars[pair_ij_km], 1);
            }
            model.add(con_17);
            // cout << "con_17_name: " << con_17_name << endl;

            // Partition for con_18: sigma{y_iq_km} = sigma{y_kj_km}
            IloExpr flow_in_q(env);
            IloExpr flow_out_q(env);
            for (auto q : subG.nodes()) {
                if (q == m || q == root[k]) continue;
                for (auto i : subG.adj_nodes_list().at(q)) {
                    pair_ij_km.first.first = i;
                    pair_ij_km.first.second = q;
                    flow_in_q += multi_flow_vars[pair_ij_km];
                }

                for (auto j : subG.adj_nodes_list().at(q)) {
                    pair_ij_km.first.first = q;
                    pair_ij_km.first.second = j;
                    flow_out_q += multi_flow_vars[pair_ij_km];
                }
            }
            // cout << con_name << endl;
            model.add(flow_in_q <= flow_out_q).setName(con_name);

            // cout << con_name << endl;
            model.add(flow_in_q >= flow_out_q).setName(con_name);
        }

        // Add constraint (19-20)
        pair_i_k.second = k;
        for (auto m : subG.nodes()) {
            if (m == root[k]) continue;
            pair_ij_km.second.second = m;
            for (auto arc : subG.arcs()) {
                pair_ij_km.first.first = arc.first;
                pair_ij_km.first.second = arc.second;

                // constraint 18
                pair_i_k.first = arc.first;
                // cout << con_name << endl;
                model.add(multi_flow_vars[pair_ij_km] <=
                          partition_node_vars[pair_i_k]);

                // constraint 19
                pair_i_k.first = arc.second;
                // cout << con_name << endl;
                model.add(multi_flow_vars[pair_ij_km] <=
                          partition_node_vars[pair_i_k]);
            }
        }
    }

    cout << "end" << endl;
}

void SmpSolver::build_problem_mcf_terminal() {
    /*****************/
    /* Add variables */
    /*****************/

    cout << "Begin to build problem mcf terminal..." << endl;
    cout << "Begin to add variables..." << endl;

    IloEnv env = model.getEnv();

    char var_name[255];
    char con_name[255];
    int V_k_size;
    SUB_Graph subG;
    map<INDEX, NODE_SET> V_k_set = G->v_set();
    map<INDEX, NODE_SET> T_k_set = G->t_set();
    pair<NODE, INDEX> pair_i_k;
    pair<NODE_PAIR, NODE_PAIR> pair_ij_kt;
    pair<NODE, NODE_PAIR> z_v_kt;
    map<INDEX, NODE> root;
    IloNumVar temp_var;

    for (auto k : G->p_set()) {
        subG = G->get_subgraph()[k];
        z_v_kt.second.first = k;

        // For each T_k, choose a root r_k
        auto firstElement = T_k_set[k].begin();
        root[k] = *firstElement;
        // cout << "r" << k << " = " << root[k] << endl;

        // Add z_v_kt, float
        for (auto v : subG.nodes()) {
            z_v_kt.first = v;
            for (auto t : T_k_set[k]) {
                if (t == root[k]) continue;
                z_v_kt.second.second = t;
                IloNumVar var;
                snprintf(var_name, 255, "z_%d_%d%d", v, k, t);
                var = IloNumVar(env, 0.0, 1.0, IloNumVar::Float, var_name);
                path_flow_vars[z_v_kt] = var;
                model.add(var);
                // printInfo(var);
            }
        }

        // Add y_ij_km (float)
        pair_ij_kt.second.first = k;
        for (auto t : subG.t_set()) {
            if (t == root[k]) continue;
            pair_ij_kt.second.second = t;
            for (auto arc : subG.arcs()) {
                IloNumVar var;
                snprintf(var_name, 255, "y_%d,%d^%d,%d", arc.first, arc.second,
                         k, t);
                var = IloNumVar(env, 0.0, IloInfinity, IloNumVar::Float,
                                var_name);
                pair_ij_kt.first = arc;
                multi_flow_vars[pair_ij_kt] = var;
                model.add(var);
                // printInfo(var);
            }
        }
    }
    /*******************/
    /* Add constraints */
    /*******************/

    cout << "Begin to add constraints..." << endl;

    for (auto k : G->p_set()) {
        z_v_kt.second.first = k;
        pair_ij_kt.second.first = k;
        subG = G->get_subgraph()[k];

        for (auto t : T_k_set[k]) {
            // cout << "t = " << t << endl;
            if (t == root[k]) continue;
            z_v_kt.second.second = t;
            pair_ij_kt.second.second = t;

            // Add cons (2) and (5)
            IloRange con_2 = IloRange(env, 1, 1);
            IloRange con_5 = IloRange(env, 1, 1);
            z_v_kt.first = root[k];
            con_2.setLinearCoef(path_flow_vars[z_v_kt], 1);
            model.add(con_2);
            z_v_kt.first = t;
            con_5.setLinearCoef(path_flow_vars[z_v_kt], 1);
            model.add(con_5);

            // Add cons (1)
            for (auto v : V_k_set[k]) {
                if (T_k_set[k].find(v) != T_k_set[k].end()) continue;
                z_v_kt.first = v;
                snprintf(con_name, 255, "%s >= %s (1)",
                         primal_node_vars[v].getName(),
                         path_flow_vars[z_v_kt].getName());
                model.add(primal_node_vars[v] >= path_flow_vars[z_v_kt]);
                // cout << con_name << endl;
            }

            // Add cons (3)
            string con_3_name = "";
            IloRange con_3 = IloRange(env, 0, 0);
            for (auto i : subG.adj_nodes_list().at(root[k])) {
                pair_ij_kt.first.first = i;
                pair_ij_kt.first.second = root[k];
                con_3.setLinearCoef(multi_flow_vars[pair_ij_kt], 1);
                con_3_name =
                    con_3_name + " + " + multi_flow_vars[pair_ij_kt].getName();
            }
            con_3.setName(con_3_name.c_str());
            model.add(con_3);
            // cout << "cons (3): " << con_3_name << " = " << 0 << endl;

            // Add cons (4)
            string con_4_name = "";
            IloRange con_4 = IloRange(env, 1, 1);
            for (auto i : subG.adj_nodes_list().at(root[k])) {
                pair_ij_kt.first.first = root[k];
                pair_ij_kt.first.second = i;
                con_4.setLinearCoef(multi_flow_vars[pair_ij_kt], 1);
                con_4_name =
                    con_4_name + " + " + multi_flow_vars[pair_ij_kt].getName();
            }
            con_4.setName(con_4_name.c_str());
            model.add(con_4);
            // cout << "cons (4): " << con_4_name << " = " << 1 << endl;

            // Add cons(6) and cons(7)
            string con_6_name = "";
            string con_7_name = "";
            IloRange con_6 = IloRange(env, 1, 1);
            IloRange con_7 = IloRange(env, 0, 0);

            for (auto i : subG.adj_nodes_list().at(t)) {
                pair_ij_kt.first.first = i;
                pair_ij_kt.first.second = t;
                con_6.setLinearCoef(multi_flow_vars[pair_ij_kt], 1);
                con_6_name =
                    con_6_name + " + " + multi_flow_vars[pair_ij_kt].getName();
                model.add(con_6);
                // cout << "cons(6)" << con_6_name << " = " << 1 << endl;

                pair_ij_kt.first.first = t;
                pair_ij_kt.first.second = i;
                con_7.setLinearCoef(multi_flow_vars[pair_ij_kt], 1);
                con_7_name =
                    con_7_name + " + " + multi_flow_vars[pair_ij_kt].getName();
                model.add(con_7);
                // cout << "cons(7)" << con_7_name << " = " << 0 << endl;
            }

            // Add cons(8) and cons(9)
            for (auto v : subG.nodes()) {
                if (v == root[k] || v == t) continue;

                string con_8_name = "";
                string con_9_name = "";
                IloExpr flow_in_v(env);
                IloExpr flow_out_v(env);

                for (auto i : subG.adj_nodes_list().at(v)) {
                    pair_ij_kt.first.first = i;
                    pair_ij_kt.first.second = v;
                    flow_in_v += multi_flow_vars[pair_ij_kt];
                    con_8_name = con_8_name + " + " +
                                 multi_flow_vars[pair_ij_kt].getName();

                    pair_ij_kt.first.first = v;
                    pair_ij_kt.first.second = i;
                    flow_out_v += multi_flow_vars[pair_ij_kt];
                    con_9_name = con_9_name + " + " +
                                 multi_flow_vars[pair_ij_kt].getName();
                }
                z_v_kt.first = v;

                model.add(flow_in_v >= path_flow_vars[z_v_kt]);
                model.add(flow_in_v <= path_flow_vars[z_v_kt]);

                model.add(flow_out_v >= path_flow_vars[z_v_kt]);
                model.add(flow_out_v <= path_flow_vars[z_v_kt]);

                // cout << "cons (8): " << con_8_name << " = " <<
                // path_flow_vars[z_v_kt].getName() << endl; cout << "cons (9):
                // "
                // << con_9_name << " = " << path_flow_vars[z_v_kt].getName() <<
                // endl;
            }
        }
    }
    cout << " end " << endl;
}

/* Steiner forest formulation */
void SmpSolver::build_problem_steiner() {
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

    // Add varaible y_ij_k
    int idx = 0;
    x_vararray = IloNumVarArray(env);
    for (auto k : G->p_set()) {
        pair_ij_k.second = k;
        subG = G->get_subgraph()[k];
        for (auto& arc : subG.arcs()) {
            IloNumVar var;
            snprintf(var_name, 255, "y_%d%d_%d", arc.first, arc.second, k);
            pair_ij_k.first = arc;
            if (relax)
                var = IloNumVar(env, 0, 1, IloNumVar::Float, var_name);
            else
                var = IloNumVar(env, 0, 1, IloNumVar::Int, var_name);
            edge_vars[pair_ij_k] = var;
            x_vararray.add(var);
            x_varindex_steiner[pair_ij_k] = idx++;
            model.add(var);
            // printInfo(var);
        }

        // For each T_k, choose a root r_k
        auto firstElement = T_k_set[k].begin();
        Steiner_root[k] = *firstElement;
        // cout << "r" << k << " = " << Steiner_root[k] << endl;
    }

    /*******************/
    /* Add constraints */
    /*******************/

    // cout << "Begin to Add the Constraint..." << endl;
    for (auto k : G->p_set()) {
        // cout << "For part " << k << "..." << endl;
        subG = G->get_subgraph()[k];

        // Partition for con_23: x_j >= sigma{y_ij_k}
        for (auto j : subG.nodes()) {
            if (std::find(G->t_total().begin(), G->t_total().end(), j) !=
                G->t_total().end())
                continue;
            string con_23_name = "";
            IloExpr sigma_vars(env);
            pair_ij_k.second = k;
            for (auto i : subG.adj_nodes_list().at(j)) {
                pair_ij_k.first.first = i;
                pair_ij_k.first.second = j;
                sigma_vars += edge_vars[pair_ij_k];
                con_23_name =
                    con_23_name + " + " + edge_vars[pair_ij_k].getName();
            }
            snprintf(con_name, 255, "con_23: %s >= %s",
                     primal_node_vars[j].getName(), con_23_name.c_str());
            // cout << con_name << endl;
            model.add(primal_node_vars[j] >= sigma_vars);
        }

        // Partition for con_24: sigma{y_it_k = 1}
        for (auto t : T_k_set[k]) {
            pair_ij_k.second = k;

            string con_24_name = "";
            IloRange con_24 = IloRange(env, 1, 1);
            for (auto i : subG.adj_nodes_list().at(t)) {
                pair_ij_k.first.first = i;
                pair_ij_k.first.second = t;
                con_24.setLinearCoef(edge_vars[pair_ij_k], 1);
                con_24_name =
                    con_24_name + " + " + edge_vars[pair_ij_k].getName();
            }
            con_24.setName(con_24_name.c_str());
            model.add(con_24);
            // cout << "con_24: " << con_24_name << "= 1" << endl;
        }

        // Partition for con_25: sigma{y_iq_k <= 1}
        for (auto q : V_k_set[k]) {
            if (std::find(T_k_set[k].begin(), T_k_set[k].end(), q) !=
                T_k_set[k].end())
                continue;
            pair_ij_k.second = k;

            string con_25_name = "";
            IloRange con_25 = IloRange(env, 0, 1);
            for (auto i : subG.adj_nodes_list().at(q)) {
                pair_ij_k.first.first = i;
                pair_ij_k.first.second = q;
                con_25.setLinearCoef(edge_vars[pair_ij_k], 1);
                con_25_name =
                    con_25_name + " + " + edge_vars[pair_ij_k].getName();
            }
            con_25.setName(con_25_name.c_str());
            model.add(con_25);
            // cout << "con_25: 0 <= " << con_25_name << " <= 1" << endl;
        }

        // Partition for con_26: y_ij_k + y_ji_k <= 1;
        for (auto arc : subG.arcs()) {
            string con_26_name = "";
            pair_ij_k.second = k;
            NODE i = arc.first, j = arc.second;
            if (i == Steiner_root[k] || j == Steiner_root[k]) continue;
            IloRange con_26 = IloRange(env, 0, 1);

            pair_ij_k.first = make_pair(i, j);
            con_26.setLinearCoef(edge_vars[pair_ij_k], 1);
            con_26_name = con_26_name + edge_vars[pair_ij_k].getName();

            pair_ij_k.first = make_pair(j, i);
            con_26.setLinearCoef(edge_vars[pair_ij_k], 1);
            con_26_name = con_26_name + " + " + edge_vars[pair_ij_k].getName();

            con_26.setName(con_26_name.c_str());
            // cout << "con_26: " << con_26_name << " = 1" << endl;
            model.add(con_26);
        }
        // cout << endl;
    }
    cout << "end" << endl;
}

/* node separator formulation */
void SmpSolver::build_problem_ns() {
    /*****************/
    /* Add variables */
    /*****************/

    cout << "Begin to build problem NS..." << endl;
    cout << "Begin to add variables..." << endl;
    IloEnv env = model.getEnv();

    char var_name[255];
    char con_name[255];
    int V_k_size;
    SUB_Graph subG;
    map<INDEX, NODE_SET> V_k_set = G->v_set();
    map<INDEX, NODE_SET> T_k_set = G->t_set();
    pair<NODE, INDEX> pair_i_k;
    pair<NODE, INDEX> pair_j_k;
    pair<NODE_PAIR, INDEX> pair_ij_k;
    IloNumVar temp_var;

    // Add  x_i^k (binary), for each node in G[V_k]:
    int idx = 0;
    x_vararray = IloNumVarArray(env);
    for (auto k : G->p_set()) {
        V_k_size = static_cast<int>(V_k_set[k].size());
        subG = G->get_subgraph()[k];
        pair_i_k.second = k;

        // Pre determine the value of x_i_k
        map<NODE, NODE_SET> ValidTerminals;
        // Node i valid to which terminals
        // -1 indocates i has no adj terminals
        // -2 means i has adj terminals but invalid for all
        map<NODE, NODE_SET> UsefulNodes;
        for (auto i : subG.nodes()) {
            if (subG.AdjTerminalNodes().at(i).size() == 0) {
                ValidTerminals[i].insert(-1);
            } else {
                bool flag = 0;
                for (auto t : subG.AdjTerminalNodes().at(i)) {
                    for (auto j : subG.adj_nodes_list().at(i)) {
                        if (j == t) {
                            continue;
                        }
                        if (std::find(subG.adj_nodes_list().at(t).begin(),
                                      subG.adj_nodes_list().at(t).end(),
                                      j) == subG.adj_nodes_list().at(t).end()) {
                            // NODE i has one adj nodes j non exists for
                            // terminal t adj nodes ->
                            // NODE i is valid for terminal t
                            ValidTerminals[i].insert(t);
                            UsefulNodes[t].insert(i);
                            flag = 1;
                        }
                    }
                }
                if (flag == 0) {
                    ValidTerminals[i].insert(-2);
                }
            }
        }

        for (auto i : subG.nodes()) {
            IloNumVar var;
            snprintf(var_name, 255, "x_%d^%d", i, k);
            if (T_k_set[k].find(i) != T_k_set[k].end()) {
                if (relax)
                    var = IloNumVar(env, 1, 1, IloNumVar::Float, var_name);
                else
                    var = IloNumVar(env, 1, 1, IloNumVar::Int, var_name);
            } else if (ValidTerminals[i].size() == 1 &&
                       *ValidTerminals[i].begin() == -2) {
                if (relax)
                    var = IloNumVar(env, 0, 0, IloNumVar::Float, var_name);
                else
                    var = IloNumVar(env, 0, 0, IloNumVar::Int, var_name);
            } else {
                if (relax)
                    var = IloNumVar(env, 0, 1, IloNumVar::Float, var_name);
                else
                    var = IloNumVar(env, 0, 1, IloNumVar::Int, var_name);
            }
            pair_i_k.first = i;
            partition_node_vars[pair_i_k] = var;
            x_vararray.add(var);
            x_varindex_ns[pair_i_k] = idx++;
            model.add(var);
            // printInfo(var);
        }

        // For each T_k, choose a root r_k
        auto firstElement = T_k_set[k].begin();
        ns_root[k] = *firstElement;

        // Add cons: sigma{x_j_k} >= 1, or 2x_i_k
        pair_i_k.second = k;
        pair_j_k.second = k;

        for (auto i : subG.nodes()) {
            pair_i_k.first = i;
            string cons32_left = "";
            IloExpr sigma_vars(env);

            if (subG.CheckNodeIsTerminal().at(i)) {
                for (auto j : subG.adj_nodes_list().at(i)) {
                    if (!subG.CheckNodeIsTerminal().at(j) &&
                        UsefulNodes[i].find(j) == UsefulNodes[i].end()) {
                        continue;
                    } else {
                        pair_j_k.first = j;
                        sigma_vars += partition_node_vars[pair_j_k];
                        cons32_left = cons32_left + " + " +
                                      partition_node_vars[pair_j_k].getName();
                    }
                }
                model.add(sigma_vars >= 1);
                // cout << "constraint(32): " << cons32_left << ">= 1" << endl;
            } else {
                if (ValidTerminals[i].size() == 1 &&
                    *ValidTerminals[i].begin() == -2) {
                    continue;
                }
                for (auto j : subG.adj_nodes_list().at(i)) {
                    pair_j_k.first = j;
                    sigma_vars += partition_node_vars[pair_j_k];
                    cons32_left = cons32_left + " + " +
                                  partition_node_vars[pair_j_k].getName();
                }
                model.add(sigma_vars >= 2 * partition_node_vars[pair_i_k]);
                // cout << "constraint(32): " << cons32_left << ">= 2*"
                //<< partition_node_vars[pair_i_k].getName() << endl;
            }
        }
    }

    /*******************/
    /* Add constraints */
    /*******************/
    cout << "Begin to Add the Constraint..." << endl;

    // Begin to add cons 29: x_i >= x_i_k
    for (auto i : G->v_total()) {
        if (std::find(G->t_total().begin(), G->t_total().end(), i) !=
            G->t_total().end())
            continue;
        for (auto k : G->nodes_of_v().at(i)) {
            pair_i_k.first = i;
            pair_i_k.second = k;
            snprintf(con_name, 255, "%s >= %s (5)",
                     primal_node_vars[i].getName(),
                     partition_node_vars[pair_i_k].getName());
            model.add(primal_node_vars[i] >= partition_node_vars[pair_i_k])
                .setName(con_name);
            // cout << "constraint(29):  " << con_name << endl;
        }
    }

    // Add cut pool constraint
    /*int CutPoolSize = cutpool.cutPoolLhs().size();
    for (int i = 0; i < CutPoolSize; i++) {
        model.add(cutpool.cutPoolLhs()[i] >= 1);
    }*/

    generate_ns_mincut_graph(G, ns_root);

    return;
}

/* frequent used function */
void SmpSolver::printInfo(IloNumVar var) {
    printf("add %s: (%f, %f), type: %d\n", var.getName(), var.getLb(),
           var.getUb(), var.getType());
}

void SmpSolver::print_to_file() {
    // begin to write the information into the file
    string store = filename;
    string graph_id = "";
    string newfilename = "";
    if (isdigit(store[store.size() - 6]))
        graph_id = store[store.size() - 6] + store[store.size() - 5];
    else
        graph_id = store[store.size() - 5];
    while (store[store.size() - 1] != '\\') {
        newfilename.push_back(*store.rbegin());
        store.pop_back();
    }
    std::reverse(newfilename.begin(), newfilename.end());
    switch (formulation) {
        case SCF: {
            store = store + "1_SCF";
            break;
        }
        case MCF: {
            store = store + "1_MCF";
            break;
        }
        case STEINER: {
            store = store + "1_STEINER";
            break;
        }
        case NS: {
            store = store + "1_NS";
            break;
        }
    }
    if (relax) store = store + "_relax";
    store = store + ".txt";

    //[Gap] [time] [Status] [Value] [Nodes number] [User number]
    ofstream flow(store, ios::app);
    flow.setf(ios::left, ios::adjustfield);
    flow << setw(LSPACING) << newfilename;  // graph number
    flow << setw(SPACING) << elapsed_time;
    flow << setw(SPACING) << cplex.getNnodes();
    flow << setw(SPACING) << cplex.getNcuts(IloCplex::CutUser);
    flow << setw(SPACING) << cplex.getMIPRelativeGap();
    flow << setw(SPACING) << cplex.getObjValue();
    flow << setw(SPACING) << cplex.getStatus();

    // flow << setw(SPACING) << formulation ;
    // flow << setw(SPACING) << callbackOption ;
    // flow << setw(SPACING) << ns_sep_opt ;
    // flow << setw(SPACING) << time_limit ;
    // flow << setw(SPACING) << max_cuts_lazy;
    // flow << setw(SPACING) << tol_lazy;
    // flow << setw(SPACING) << max_cuts_user;
    // flow << setw(SPACING) << tol_user;
    switch (callbackOption) {
        case 0:
            flow << setw(LSPACING) << "NULL";
            break;
        case 1:
            switch (lazy_sep_opt) {
                case 0:
                    flow << setw(LSPACING + 6) << "L(1-m)";
                    flow << "(" << max_cuts_lazy << ", " << tol_lazy
                         << setw(SPACING + 9) << ")";
                    break;
                case 1:
                    flow << setw(LSPACING + 6) << "L(m-m)";
                    flow << "(" << max_cuts_lazy << ", " << tol_lazy
                         << setw(SPACING + 9) << ")";
                    break;
            }

            break;
        case 2:
            flow << setw(SPACING) << "U";
            break;
        case 3:
            switch (lazy_sep_opt) {
                case 0:
                    flow << "L(1-m)";
                    break;
                case 1:
                    flow << "L(m-m)";
                    break;
            }

            switch (ns_sep_opt) {
                case 1:
                    flow << setw(LSPACING) << "U(SCC-MinCut)";
                    break;
                case 0:
                    flow << setw(LSPACING) << "U(MinCut)";
                    break;
            }

            flow << "(" << max_cuts_lazy << ", " << tol_lazy << ";"
                 << max_cuts_user << "," << tol_user << setw(SPACING) << ")";

            break;
        default:
            break;
    }
    flow << endl;
}

/*******************************************************/
/************************** LB *************************/
/*******************************************************/

LBSolver::LBSolver(IloEnv env, std::shared_ptr<Graph> g_ptr,
                   SmpForm _formulation, int _callbackOption, bool _relax,
                   bool _ns_sep_opt, int _LB_MaxRestarts, int _LB_MaxIter,
                   int _Rmin, int _Rmax, int _BCSolNum, int _BCTime,
                   double _epsilon_lazy, double _epsilon_user,
                   int _max_cuts_lazy, int _max_cuts_user, string _filename,
                   int _MIPDisplayLevel, int LB_CP_Option_, int lazy_sep_opt_) {
    LBmodel = IloModel(env);
    LBobjective = IloObjective();

    formulation = _formulation;
    G = g_ptr;
    callbackOption = _callbackOption;
    relax = _relax;
    ns_sep_opt = _ns_sep_opt;

    LB_MaxRestarts = _LB_MaxRestarts;
    LB_MaxIter = _LB_MaxIter;
    Rmin = _Rmin;
    Rmax = _Rmax;
    BCSolNum = _BCSolNum;
    BCTime = _BCTime;

    max_cuts_lazy = _max_cuts_lazy;
    max_cuts_user = _max_cuts_user;
    tol_lazy = _epsilon_lazy;
    tol_user = _epsilon_user;
    filename = _filename;
    MIPDisplayLevel = _MIPDisplayLevel;
    LB_CP_Option = LB_CP_Option_;
    lazy_sep_opt = lazy_sep_opt_;

    /* Add x_i variables: primal_node_vars */
    char var_name[255];
    set<NODE> T = G->t_total();
    int idx = 0;
    x_vararray_primal = IloNumVarArray(env);
    for (auto& node : G->nodes()) {
        IloNumVar var;
        snprintf(var_name, 255, "x_%d", node);
        // set x_i = 1 if i belongs to T, others to {0,1}
        if (T.find(node) != T.end())  // Constraint(2)
        {
            if (relax)
                var = IloNumVar(env, 1, 1, IloNumVar::Float, var_name);
            else
                var = IloNumVar(env, 1, 1, IloNumVar::Bool, var_name);
        } else {
            if (relax)
                var = IloNumVar(env, 0, 1, IloNumVar::Float, var_name);
            else
                var = IloNumVar(env, 0, 1, IloNumVar::Bool, var_name);
        }
        primal_node_vars[node] = var;
        x_vararray_primal.add(var);
        x_varindex_ns_primal[node] = idx++;
        // printInfo(var);
    }

    cout << "Begin to build Local Branch problem..." << endl;

    build_problem_ns_simplifer();

    LBcplex = IloCplex(LBmodel);
    LBcplex.setParam(IloCplex::IntSolLim, BCSolNum);
    LBcplex.setParam(IloCplex::MIPDisplay, MIPDisplayLevel);
    LBcplex.setParam(IloCplex::TiLim, BCTime);
    LBcplex.setParam(IloCplex::Threads, 0);

    fianlsolveflag = 0;

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
                    LBcplex.use(NS_StrongComponentLazyCallback(
                        env, G, partition_node_vars, x_vararray, x_varindex_ns,
                        tol_lazy, max_cuts_lazy, formulation, ns_root,
                        x_vararray_primal, x_varindex_ns_primal, lazy_sep_opt));
                    break;
                case 2:
                    LBcplex.use(NS_CutCallback(
                        env, G, partition_node_vars, x_vararray, x_varindex_ns,
                        tol_user, max_cuts_user, formulation, ns_root,
                        ns_sep_opt, LB_CP_Option, fianlsolveflag,
                        lazy_sep_opt));
                    break;
                case 3:
                    LBcplex.use(NS_StrongComponentLazyCallback(
                        env, G, partition_node_vars, x_vararray, x_varindex_ns,
                        tol_lazy, max_cuts_lazy, formulation, ns_root,
                        x_vararray_primal, x_varindex_ns_primal, lazy_sep_opt));
                    LBcplex.use(NS_CutCallback(
                        env, G, partition_node_vars, x_vararray, x_varindex_ns,
                        tol_user, max_cuts_user, formulation, ns_root,
                        ns_sep_opt, LB_CP_Option, fianlsolveflag,
                        lazy_sep_opt));
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

/* Final build problem */
void LBSolver::build_problem_ns_final() {
    /*****************/
    /* Add variables */
    /*****************/
    primal_node_vars.clear();
    partition_node_vars.clear();
    ns_root.clear();
    x_varindex_ns.clear();
    x_varindex_ns_primal.clear();

    cout << "Begin to build Final NS Model..." << endl;
    cout << "Begin to add variables..." << endl;

    IloEnv env = FLBmodel.getEnv();

    /* Add x_i variables: primal_node_vars */
    char var_name[255];
    set<NODE> T = G->t_total();
    int idx = 0;
    x_vararray_primal = IloNumVarArray(env);
    for (auto& node : G->nodes()) {
        IloNumVar var;
        snprintf(var_name, 255, "x_%d", node);
        // set x_i = 1 if i belongs to T, others to {0,1}
        if (T.find(node) != T.end())  // Constraint(2)
        {
            if (relax)
                var = IloNumVar(env, 1, 1, IloNumVar::Float, var_name);
            else
                var = IloNumVar(env, 1, 1, IloNumVar::Bool, var_name);
        } else {
            if (relax)
                var = IloNumVar(env, 0, 1, IloNumVar::Float, var_name);
            else
                var = IloNumVar(env, 0, 1, IloNumVar::Bool, var_name);
        }
        primal_node_vars[node] = var;
        x_vararray_primal.add(var);
        x_varindex_ns_primal[node] = idx++;
        // printInfo(var);
    }

    char con_name[255];
    int V_k_size;
    SUB_Graph subG;
    map<INDEX, NODE_SET> V_k_set = G->v_set();
    map<INDEX, NODE_SET> T_k_set = G->t_set();
    pair<NODE, INDEX> pair_i_k;
    pair<NODE, INDEX> pair_j_k;
    pair<NODE_PAIR, INDEX> pair_ij_k;
    IloNumVar temp_var;

    // Add  x_i^k (binary), for each node in G[V_k]:
    idx = 0;
    x_vararray = IloNumVarArray(env);
    for (auto k : G->p_set()) {
        V_k_size = static_cast<int>(V_k_set[k].size());
        subG = G->get_subgraph()[k];
        pair_i_k.second = k;

        // Pre determine the value of x_i_k
        map<NODE, NODE_SET> ValidTerminals;
        // Node i valid to which terminals
        // -1 indocates i has no adj terminals
        // -2 means i has adj terminals but invalid for all
        map<NODE, NODE_SET> UsefulNodes;
        for (auto i : subG.nodes()) {
            if (subG.AdjTerminalNodes().at(i).size() == 0) {
                ValidTerminals[i].insert(-1);
            } else {
                bool flag = 0;
                for (auto t : subG.AdjTerminalNodes().at(i)) {
                    for (auto j : subG.adj_nodes_list().at(i)) {
                        if (j == t) {
                            continue;
                        }
                        if (std::find(subG.adj_nodes_list().at(t).begin(),
                                      subG.adj_nodes_list().at(t).end(),
                                      j) == subG.adj_nodes_list().at(t).end()) {
                            // NODE i has one adj nodes j non exists for
                            // terminal t adj nodes ->
                            // NODE i is valid for terminal t
                            ValidTerminals[i].insert(t);
                            UsefulNodes[t].insert(i);
                            flag = 1;
                        }
                    }
                }
                if (flag == 0) {
                    ValidTerminals[i].insert(-2);
                }
            }
        }

        for (auto i : subG.nodes()) {
            IloNumVar var;
            snprintf(var_name, 255, "x_%d^%d", i, k);
            if (T_k_set[k].find(i) != T_k_set[k].end()) {
                if (relax)
                    var = IloNumVar(env, 1, 1, IloNumVar::Float, var_name);
                else
                    var = IloNumVar(env, 1, 1, IloNumVar::Int, var_name);
            } else if (ValidTerminals[i].size() == 1 &&
                       *ValidTerminals[i].begin() == -2) {
                if (relax)
                    var = IloNumVar(env, 0, 0, IloNumVar::Float, var_name);
                else
                    var = IloNumVar(env, 0, 0, IloNumVar::Int, var_name);
            } else {
                if (relax)
                    var = IloNumVar(env, 0, 1, IloNumVar::Float, var_name);
                else
                    var = IloNumVar(env, 0, 1, IloNumVar::Int, var_name);
            }
            pair_i_k.first = i;
            partition_node_vars[pair_i_k] = var;
            x_vararray.add(var);
            x_varindex_ns[pair_i_k] = idx++;
            FLBmodel.add(var);
            // printInfo(var);
        }

        // For each T_k, choose a root r_k
        auto firstElement = T_k_set[k].begin();
        ns_root[k] = *firstElement;

        // Add cons: sigma{x_j_k} >= 1, or 2x_i_k
        pair_i_k.second = k;
        pair_j_k.second = k;

        for (auto i : subG.nodes()) {
            pair_i_k.first = i;
            string cons32_left = "";
            IloExpr sigma_vars(env);

            if (subG.CheckNodeIsTerminal().at(i)) {
                for (auto j : subG.adj_nodes_list().at(i)) {
                    if (!subG.CheckNodeIsTerminal().at(j) &&
                        UsefulNodes[i].find(j) == UsefulNodes[i].end()) {
                        continue;
                    } else {
                        pair_j_k.first = j;
                        sigma_vars += partition_node_vars[pair_j_k];
                        cons32_left = cons32_left + " + " +
                                      partition_node_vars[pair_j_k].getName();
                    }
                }
                FLBmodel.add(sigma_vars >= 1);
                // cout << "constraint(32): " << cons32_left << ">= 1" << endl;
            } else {
                if (ValidTerminals[i].size() == 1 &&
                    *ValidTerminals[i].begin() == -2) {
                    continue;
                }
                for (auto j : subG.adj_nodes_list().at(i)) {
                    pair_j_k.first = j;
                    sigma_vars += partition_node_vars[pair_j_k];
                    cons32_left = cons32_left + " + " +
                                  partition_node_vars[pair_j_k].getName();
                }
                FLBmodel.add(sigma_vars >= 2 * partition_node_vars[pair_i_k]);
                // cout << "constraint(32): " << cons32_left << ">= 2*"
                //<< partition_node_vars[pair_i_k].getName() << endl;
            }
        }
    }

    /*******************/
    /* Add constraints */
    /*******************/
    cout << "Begin to Add the Constraint..." << endl;

    // Begin to add cons 29: x_i >= x_i_k
    for (auto i : G->v_total()) {
        if (std::find(G->t_total().begin(), G->t_total().end(), i) !=
            G->t_total().end())
            continue;
        for (auto k : G->nodes_of_v().at(i)) {
            pair_i_k.first = i;
            pair_i_k.second = k;
            snprintf(con_name, 255, "%s >= %s (5)",
                     primal_node_vars[i].getName(),
                     partition_node_vars[pair_i_k].getName());
            FLBmodel.add(primal_node_vars[i] >= partition_node_vars[pair_i_k])
                .setName(con_name);
            // cout << "constraint(29):  " << con_name << endl;
        }
    }

    generate_ns_mincut_graph(G, ns_root);

    return;
}

/* node separator formulation */
void LBSolver::build_problem_ns_simplifer() {
    /*****************/
    /* Add variables */
    /*****************/

    cout << "Begin to build problem NS..." << endl;
    cout << "Begin to add variables..." << endl;
    IloEnv env = LBmodel.getEnv();

    char var_name[255];
    char con_name[255];
    int V_k_size;
    SUB_Graph subG;
    map<INDEX, NODE_SET> V_k_set = G->v_set();
    map<INDEX, NODE_SET> T_k_set = G->t_set();
    pair<NODE, INDEX> pair_i_k;
    pair<NODE, INDEX> pair_j_k;
    pair<NODE_PAIR, INDEX> pair_ij_k;
    IloNumVar temp_var;

    // Add  x_i^k (binary), for each node in G[V_k]:
    int idx = 0;
    x_vararray = IloNumVarArray(env);
    for (auto k : G->p_set()) {
        V_k_size = static_cast<int>(V_k_set[k].size());
        subG = G->get_subgraph()[k];
        pair_i_k.second = k;

        // Pre determine the value of x_i_k
        map<NODE, NODE_SET> ValidTerminals;
        // Node i valid to which terminals
        // -1 indocates i has no adj terminals
        // -2 means i has adj terminals but invalid for all
        map<NODE, NODE_SET> UsefulNodes;
        for (auto i : subG.nodes()) {
            if (subG.AdjTerminalNodes().at(i).size() == 0) {
                ValidTerminals[i].insert(-1);
            } else {
                bool flag = 0;
                for (auto t : subG.AdjTerminalNodes().at(i)) {
                    for (auto j : subG.adj_nodes_list().at(i)) {
                        if (j == t) {
                            continue;
                        }
                        if (std::find(subG.adj_nodes_list().at(t).begin(),
                                      subG.adj_nodes_list().at(t).end(),
                                      j) == subG.adj_nodes_list().at(t).end()) {
                            // NODE i has one adj nodes j non exists for
                            // terminal t adj nodes ->
                            // NODE i is valid for terminal t
                            ValidTerminals[i].insert(t);
                            UsefulNodes[t].insert(i);
                            flag = 1;
                        }
                    }
                }
                if (flag == 0) {
                    ValidTerminals[i].insert(-2);
                }
            }
        }

        for (auto i : subG.nodes()) {
            IloNumVar var;
            snprintf(var_name, 255, "x_%d^%d", i, k);
            if (T_k_set[k].find(i) != T_k_set[k].end()) {
                if (relax)
                    var = IloNumVar(env, 1, 1, IloNumVar::Float, var_name);
                else
                    var = IloNumVar(env, 1, 1, IloNumVar::Int, var_name);
            } else if (ValidTerminals[i].size() == 1 &&
                       *ValidTerminals[i].begin() == -2) {
                if (relax)
                    var = IloNumVar(env, 0, 0, IloNumVar::Float, var_name);
                else
                    var = IloNumVar(env, 0, 0, IloNumVar::Int, var_name);
            } else {
                if (relax)
                    var = IloNumVar(env, 0, 1, IloNumVar::Float, var_name);
                else
                    var = IloNumVar(env, 0, 1, IloNumVar::Int, var_name);
            }
            pair_i_k.first = i;
            partition_node_vars[pair_i_k] = var;
            x_vararray.add(var);
            x_varindex_ns[pair_i_k] = idx++;
            LBmodel.add(var);
            // printInfo(var);
        }

        // For each T_k, choose a root r_k
        auto firstElement = T_k_set[k].begin();
        ns_root[k] = *firstElement;

        // Add cons: sigma{x_j_k} >= 1, or 2x_i_k
        pair_i_k.second = k;
        pair_j_k.second = k;

        for (auto i : subG.nodes()) {
            pair_i_k.first = i;
            string cons32_left = "";
            IloExpr sigma_vars(env);

            if (subG.CheckNodeIsTerminal().at(i)) {
                for (auto j : subG.adj_nodes_list().at(i)) {
                    if (!subG.CheckNodeIsTerminal().at(j) &&
                        UsefulNodes[i].find(j) == UsefulNodes[i].end()) {
                        continue;
                    } else {
                        pair_j_k.first = j;
                        sigma_vars += partition_node_vars[pair_j_k];
                        cons32_left = cons32_left + " + " +
                                      partition_node_vars[pair_j_k].getName();
                    }
                }
                LBmodel.add(sigma_vars >= 1);
                // cout << "constraint(32): " << cons32_left << ">= 1" << endl;
            } else {
                if (ValidTerminals[i].size() == 1 &&
                    *ValidTerminals[i].begin() == -2) {
                    continue;
                }
                for (auto j : subG.adj_nodes_list().at(i)) {
                    pair_j_k.first = j;
                    sigma_vars += partition_node_vars[pair_j_k];
                    cons32_left = cons32_left + " + " +
                                  partition_node_vars[pair_j_k].getName();
                }
                LBmodel.add(sigma_vars >= 2 * partition_node_vars[pair_i_k]);
                // cout << "constraint(32): " << cons32_left << ">= 2*"
                //<< partition_node_vars[pair_i_k].getName() << endl;
            }
        }
    }

    /*******************/
    /* Add constraints */
    /*******************/
    cout << "Begin to Add the Constraint..." << endl;

    // Begin to add cons 29: x_i >= x_i_k
    for (auto i : G->v_total()) {
        if (std::find(G->t_total().begin(), G->t_total().end(), i) !=
            G->t_total().end())
            continue;
        for (auto k : G->nodes_of_v().at(i)) {
            pair_i_k.first = i;
            pair_i_k.second = k;
            snprintf(con_name, 255, "%s >= %s (5)",
                     primal_node_vars[i].getName(),
                     partition_node_vars[pair_i_k].getName());
            LBmodel.add(primal_node_vars[i] >= partition_node_vars[pair_i_k])
                .setName(con_name);
            // cout << "constraint(29):  " << con_name << endl;
        }
    }

    generate_ns_mincut_graph(G, ns_root);

    return;
}

void LBSolver::update_LB_problem() {
    /* Build objective function */
    IloEnv env = LBmodel.getEnv();
    IloExpr totalCost(env);
    double objcoeff = 0.0;
    IloNumVar var;

    for (auto& node : G->nodes()) {
        var = primal_node_vars[node];
        objcoeff = G->nodes_value().at(node);
        totalCost += var * objcoeff;
    }

    LBobjective = IloObjective(env, totalCost, IloObjective::Minimize);
    LBmodel.add(LBobjective);
}

void LBSolver::update_LB_problem_final() {
    /* Build objective function */
    IloEnv env = FLBmodel.getEnv();
    IloExpr totalCost(env);
    double objcoeff = 0.0;
    IloNumVar var;

    for (auto& node : G->nodes()) {
        var = primal_node_vars[node];
        objcoeff = G->nodes_value().at(node);
        totalCost += var * objcoeff;
    }

    FLBobjective = IloObjective(env, totalCost, IloObjective::Minimize);
    FLBmodel.add(FLBobjective);
}

void LBSolver::Floyd(map<NODE, int>& subGnodesIdx,
                     map<int, NODE>& rev_subGnodesIdx,
                     vector<vector<int>>& distance, vector<vector<int>>& path,
                     int& idx, int k) {
    const int NodesValueINF = 0x3f3f3f3f;
    SUB_Graph subG = G->get_subgraph()[k];

    // Add nodes
    idx = 0;  // Total Nodes Number
    for (auto i : subG.nodes()) {
        idx++;
        subGnodesIdx[i] = idx;
        rev_subGnodesIdx[idx] = i;
    }
    distance.resize(idx + 2);
    path.resize(idx + 2);
    for (int i = 1; i <= idx + 1; i++) {
        distance[i].resize(idx + 1, NodesValueINF);
        path[i].resize(idx + 1, -1);
    }

    // Add arc
    for (auto i : subG.nodes()) {
        int u = subGnodesIdx[i];
        distance[u][u] = 0;
    }
    for (auto arc : subG.arcs()) {
        int u = subGnodesIdx[arc.first];
        int v = subGnodesIdx[arc.second];
        distance[u][v] = subG.node_value().at(arc.second);
        distance[v][u] = subG.node_value().at(arc.first);
    }

    if (k == 1) {
        int u = subGnodesIdx[8];
        int v = subGnodesIdx[13];
        cout << distance[u][v] << endl;
    }

    // Run Floyd
    for (int K = 1; K <= idx; K++) {
        for (int i = 1; i <= idx; i++) {
            for (int j = 1; j <= idx; j++) {
                if (distance[i][j] > distance[i][K] + distance[K][j]) {
                    distance[i][j] = distance[i][K] + distance[K][j];
                    path[i][j] = K;
                }
            }
        }
    }

    return;
}

//
// void LBSolver::GenerateInitialSolution(int k) {
//    SUB_Graph subG = G->get_subgraph()[k];
//    map<NODE, int> subGnodesIdx;
//    map<int, NODE> rev_subGnodesIdx;
//    vector<vector<int>> distance;
//    vector<vector<int>> path;
//    int idx = 0;
//
//    Floyd(subGnodesIdx, rev_subGnodesIdx, distance, path, idx, k);
//
//    const int NodesValueINF = 1000 * G->nodes().size() + 1;
//
//    if (k == 1) {
//        cout << distance[subGnodesIdx[2]][subGnodesIdx[19]] << endl;
//        int ns = subGnodesIdx[2];
//        int nt = subGnodesIdx[19];
//        while (path[ns][nt] != -1) {
//            xPartSol[k][rev_subGnodesIdx[path[ns][nt]]] = 1;
//            xPrimalSol[rev_subGnodesIdx[path[ns][nt]]] = 1;
//            nt = path[ns][nt];
//            cout << rev_subGnodesIdx[nt] << " ";
//        }
//        cout << endl;
//    }
//
//    for (auto i : subG.nodes()) {
//        xPartSol[k][i] = false;
//    }
//
//    /*auto t1 = subG.t_set().begin();
//    auto t2 = subG.t_set().begin();
//    t2++;
//    while (t2 != subG.t_set().end()) {
//        int s = *t1, t = *t2;
//        xPartSol[k][s] = 1;
//        xPartSol[k][t] = 1;
//        xPrimalSol[s] = 1;
//        xPrimalSol[t] = 1;
//
//        int ns = subGnodesIdx[s][0];
//        int nt = subGnodesIdx[t][1];
//        while (path[ns][nt] != -1) {
//            xPartSol[k][rev_subGnodesIdx[path[ns][nt]]] = 1;
//            xPrimalSol[rev_subGnodesIdx[path[ns][nt]]] = 1;
//            nt = path[ns][nt];
//        }
//        t1++, t2++;
//    }*/
//
//    int s = *(subG.t_set().begin());
//    xPartSol[k][s] = 1;
//    xPrimalSol[s] = 1;
//    for (auto t : subG.t_set()) {
//        if (s == t) continue;
//        xPartSol[k][t] = 1;
//        xPrimalSol[t] = 1;
//
//        int ns = subGnodesIdx[s];
//        int nt = subGnodesIdx[t];
//        while (path[ns][nt] != -1) {
//            xPartSol[k][rev_subGnodesIdx[path[ns][nt]]] = 1;
//            xPrimalSol[rev_subGnodesIdx[path[ns][nt]]] = 1;
//            nt = path[ns][nt];
//        }
//    }
//
//    return;
//}
//

void LBSolver::GenerateInitialSolution(int k) {
    SUB_Graph subG = G->get_subgraph()[k];

    lemon::SmartDigraph g;
    lemon::SmartDigraph::ArcMap<double> costMap(g);
    lemon::SmartDigraph::NodeMap<int> nodeMap(g);

    map<NODE, lemon::SmartDigraph::Node> v_nodes;
    map<lemon::SmartDigraph::Node, NODE> rev_nodes;

    // defining the type of the Dijkstra Class
    using SptSolver = lemon::Dijkstra<lemon::SmartDigraph,
                                      lemon::SmartDigraph::ArcMap<double>>;

    lemon::SmartDigraph::Node currentNode;
    for (auto i : subG.nodes()) {
        currentNode = g.addNode();
        v_nodes[i] = currentNode;
        rev_nodes[currentNode] = i;
    }

    lemon::SmartDigraph::Arc currentArc;
    for (auto arc : subG.arcs()) {
        int sourceIndex = arc.first;
        int targetIndex = arc.second;
        // cout << arc.first << " " << arc.second << endl;
        lemon::SmartDigraph::Node sourceNode = v_nodes[sourceIndex];
        lemon::SmartDigraph::Node targetNode = v_nodes[targetIndex];

        currentArc = g.addArc(sourceNode, targetNode);
        costMap[currentArc] = 1.0 * subG.node_value().at(targetIndex);
    }

    for (auto i : HeuristicPool) {
        if (v_nodes.count(i)) {
            for (SmartDigraph::InArcIt a(g, v_nodes[i]); a != INVALID; ++a) {
                costMap[a] = 0.0;
            }
        }
    }

    /*for (SmartDigraph::NodeIt n(g); n != INVALID; ++n) {
            std::cout << rev_nodes[n] << std::endl;
    }
    cout << endl;
    for (SmartDigraph::ArcIt a(g); a != INVALID; ++a) {
            cout << rev_nodes[g.source(a)] <<" "<< rev_nodes[g.target(a)] <<
    endl;
    }*/

    // add source
    /*auto firstElement = subG.t_set().begin();
    int s = *firstElement;*/
    auto r = rand() % subG.t_set().size();  // not _really_ random
    auto it = std::begin(subG.t_set());
    std::advance(it, r);
    int s = *it;
    lemon::SmartDigraph::Node startN = v_nodes[s];  // strat point
    NODE_SET TmpHeuristic;                          // tmp delate nodes

    for (auto t : subG.t_set()) {
        if (s == t) continue;

        // set the selected nodes weight to 0
        for (auto x : TmpHeuristic) {
            for (SmartDigraph::InArcIt a(g, v_nodes[x]); a != INVALID; ++a) {
                costMap[a] = 0.0;
            }
        }
        TmpHeuristic.clear();

        lemon::SmartDigraph::Node endN = v_nodes[t];
        SptSolver spt(g, costMap);
        spt.run(startN, endN);

        std::vector<lemon::SmartDigraph::Node> path;
        for (lemon::SmartDigraph::Node v = endN; v != startN;
             v = spt.predNode(v)) {
            if (v != lemon::INVALID && spt.reached(v)) {
                path.push_back(v);
            }
        }
        path.push_back(startN);

        for (auto p = path.rbegin(); p != path.rend(); ++p) {
            int m = rev_nodes[*p];
            // std::cout << m << std::endl;
            xPartSol[k][m] = 1;
            xPrimalSol[m] = 1;
            HeuristicPool.insert(m);
            TmpHeuristic.insert(m);
        }
    }

    // for (auto t : subG.t_set()) {
    //    if (s == t) continue;
    //    lemon::SmartDigraph::Node endN = v_nodes[t];
    //
    //    std::vector<lemon::SmartDigraph::Node> path;
    //    for (lemon::SmartDigraph::Node v = endN; v != startN;
    //         v = spt.predNode(v)) {
    //        if (v != lemon::INVALID && spt.reached(v)) {
    //            path.push_back(v);
    //        }
    //    }
    //    path.push_back(startN);
    //
    //    for (auto p = path.rbegin(); p != path.rend(); ++p) {
    //        int m = rev_nodes[*p];
    //        // std::cout << m << std::endl;
    //        xPartSol[k][m] = 1;
    //        xPrimalSol[m] = 1;
    //    }
    //}
    return;
}

void LBSolver::LocalBranchSearch() {
    Final_Obj = INF;
    TOT_LB_TIME = 0.0;
    TOT_TIME = 0.0;
    LocalBranchTime = 0;

    for (int lbtime = 1; lbtime <= LB_MaxRestarts; lbtime++) {
        xPartSol.clear();
        xPrimalSol.clear();
        HeuristicPool.clear();

        // generate initial solution
        for (auto k : G->p_set()) {
            GenerateInitialSolution(k);
        }

        CheckSolution();

        // calculate objvalue
        int ObjValue = 0;
        for (auto i : xPrimalSol) {
            ObjValue += i.second * G->nodes_value().at(i.first);
        }

        double gap = 0.0;
        LocalBranch(ObjValue, gap);

        // change optimal solution
        if (ObjValue < Final_Obj) {
            Final_xPrimalSol = xPrimalSol;
            Final_xPartSol = xPartSol;
            Final_Obj = ObjValue;
            Final_gap = gap;
        }
    }

    FinalSolve();

    return;
}

void LBSolver::LocalBranch(int& ObjValue, double& gap) {
    IloEnv env = LBmodel.getEnv();

    pair<NODE, INDEX> pair_i_k;
    pair<NODE, INDEX> pair_j_k;

    int Iter = 1, R = Rmin, Rdelta = (int)(2 * (Rmax - Rmin) / LB_MaxIter);
    int MIPStartIndex = 0;
    while (Iter++ <= LB_MaxIter && R <= Rmax) {
        IloConstraintArray cons_array(env);

        // Add asymmetric constraint
        for (auto k : G->p_set()) {
            SUB_Graph subG = G->get_subgraph()[k];
            IloExpr sigma_vars(env);
            for (auto i : xPartSol[k]) {
                if (i.second == 0) continue;
                pair_i_k.first = i.first;
                pair_i_k.second = k;
                sigma_vars += (1 - partition_node_vars[pair_i_k]);
            }
            IloConstraint cons = sigma_vars <= R;
            LBmodel.add(cons);
            cons_array.add(cons);
        }

        // add cutpool constraints
        for (auto k : G->p_set()) {
            int CutPoolSize = cutpool.cutPoolLhs()[k].size();
            for (auto s : cutpool.cutPoolLhs()[k]) {
                IloExpr sigma_vars(env);
                for (auto i : s) {
                    pair_i_k.first = i;
                    pair_i_k.second = k;
                    sigma_vars += partition_node_vars[pair_i_k];
                }
                LBmodel.add(sigma_vars >= 1);
            }
        }

        // Warm start
        IloNumVarArray VarArray(env);
        IloNumArray NumArray(env);
        for (auto k : G->p_set()) {
            pair_i_k.second = k;
            for (auto i : xPartSol[k]) {
                pair_i_k.first = i.first;
                VarArray.add(partition_node_vars[pair_i_k]);
                if (i.second) {
                    NumArray.add(1);
                } else {
                    NumArray.add(0);
                }
            }
        }

        for (auto i : xPrimalSol) {
            VarArray.add(primal_node_vars[i.first]);
            if (i.second) {
                NumArray.add(1);
            } else {
                NumArray.add(0);
            }
        }
        LBcplex.addMIPStart(VarArray, NumArray);
        VarArray.end();
        NumArray.end();

        // Cplex Info
        LBcplex.setParam(IloCplex::RandomSeed, Iter);

        double start_time = LBcplex.getCplexTime();
        double start_ticks = LBcplex.getDetTime();

        LBcplex.solve();

        double elapsed_time = LBcplex.getCplexTime() - start_time;
        double elapsed_ticks = LBcplex.getDetTime() - start_ticks;

        TOT_LB_TIME += elapsed_time;
        TOT_TIME += elapsed_time;

        cout << "Solution status \t= \t" << LBcplex.getStatus() << endl;
        try {
            cout << "Objectvie value \t= \t" << LBcplex.getObjValue() << endl;
        } catch (IloException e) {
            cout << e << endl;
        }
        cout << "Elapsed time \t= \t" << elapsed_time << endl;

        // update Sol Value
        if (LBcplex.getObjValue() < ObjValue) {
            IloNumArray val = IloNumArray(env, partition_node_vars.size());
            IloNumArray val_primal = IloNumArray(env, G->nodes().size());
            LBcplex.getValues(val, x_vararray);
            LBcplex.getValues(val_primal, x_vararray_primal);

            SUB_Graph subG;
            for (auto k : G->p_set()) {
                subG = G->get_subgraph()[k];
                pair_i_k.second = k;
                for (auto i : subG.nodes()) {
                    pair_i_k.first = i;
                    xPartSol[k][i] = val[x_varindex_ns[pair_i_k]];
                }
            }
            for (auto i : G->nodes()) {
                xPrimalSol[i] = val_primal[x_varindex_ns_primal[i]];
            }
            ObjValue = LBcplex.getObjValue();
            gap = LBcplex.getMIPRelativeGap();

            // reset R
            R = Rmin;

        } else {
            R += Rdelta;
        }

        // Remove asymmetric constraint
        for (int i = 0; i < G->p_set().size(); i++) {
            LBmodel.remove(cons_array[i]);
        }
        if (LBcplex.getNMIPStarts() > 0) {
            LBcplex.deleteMIPStarts(0, LBcplex.getNMIPStarts());
        }
    }
    LocalBranchTime += (Iter - 1);
    return;
}

void LBSolver::CheckSolution() {
    cout << "-------------- Heuristic Solution Status "
            "--------------"
         << endl;
    for (auto k : G->p_set()) {
        SUB_Graph subG = G->get_subgraph()[k];
        UnionFind<NODE> forest(subG.nodes());
        for (auto arc : subG.arcs()) {
            int u = arc.first;
            int v = arc.second;
            if (xPartSol[k][u] && xPartSol[k][v]) {
                if (forest.find_set(u) != forest.find_set(v)) {
                    forest.join(u, v);
                }
            }
        }
        bool inconnect = false;
        int s = *(subG.t_set().begin());
        for (auto t : subG.t_set()) {
            if (forest.find_set(s) != forest.find_set(t)) {
                inconnect = true;
            }
        }
        if (inconnect) {
            cout << "partiton: " << k << " not connect" << endl;
            for (auto i : xPartSol[k]) {
                cout << i.first << " " << i.second << endl;
            }
        } else {
            cout << "partition: " << k << " terminal all connect" << endl;
        }
    }
    cout << "--------------------------- END "
            "-------------------------"
         << endl;
    return;
}

void LBSolver::FinalSolve() {
    cout << endl
         << "--------------------------------------- Begin Final Solve "
            "WithValue:  "
         << Final_Obj << "  -------------------------------" << endl;

    IloEnv nenv;
    FLBmodel = IloModel(nenv);
    FLBobjective = IloObjective();

    build_problem_ns_final();
    update_LB_problem_final();

    FLBcplex = IloCplex(FLBmodel);
    fianlsolveflag = 1;

    switch (callbackOption) {
        case 0:
            break;
        case 1:
            FLBcplex.use(NS_StrongComponentLazyCallback(
                nenv, G, partition_node_vars, x_vararray, x_varindex_ns,
                tol_lazy, max_cuts_lazy, formulation, ns_root,
                x_vararray_primal, x_varindex_ns_primal, lazy_sep_opt));
            break;
        case 2:
            FLBcplex.use(NS_CutCallback(
                nenv, G, partition_node_vars, x_vararray, x_varindex_ns,
                tol_user, max_cuts_user, formulation, ns_root, ns_sep_opt,
                LB_CP_Option, fianlsolveflag, lazy_sep_opt));
            break;
        case 3:
            FLBcplex.use(NS_StrongComponentLazyCallback(
                nenv, G, partition_node_vars, x_vararray, x_varindex_ns,
                tol_lazy, max_cuts_lazy, formulation, ns_root,
                x_vararray_primal, x_varindex_ns_primal, lazy_sep_opt));
            FLBcplex.use(NS_CutCallback(
                nenv, G, partition_node_vars, x_vararray, x_varindex_ns,
                tol_user, max_cuts_user, formulation, ns_root, ns_sep_opt,
                LB_CP_Option, fianlsolveflag, lazy_sep_opt));
            break;
        default:
            break;
    }

    pair<NODE, INDEX> pair_i_k;
    pair<NODE, INDEX> pair_j_k;

    // 0 use LB cutpool constraints as initial constraint
    if (!LB_CP_Option) {
        for (auto k : G->p_set()) {
            for (auto s : cutpool.cutPoolLhs()[k]) {
                IloExpr sigma_vars(nenv);
                for (auto i : s) {
                    pair_i_k.first = i;
                    pair_i_k.second = k;
                    sigma_vars += partition_node_vars[pair_i_k];
                }
                FLBmodel.add(sigma_vars >= 1);
            }
        }
    }

    // Warm start
    IloNumVarArray VarArray(nenv);
    IloNumArray NumArray(nenv);
    for (auto k : G->p_set()) {
        pair_i_k.second = k;
        for (auto i : Final_xPartSol[k]) {
            pair_i_k.first = i.first;
            VarArray.add(partition_node_vars[pair_i_k]);
            if (i.second) {
                NumArray.add(1);
            } else {
                NumArray.add(0);
            }
        }
    }
    for (auto i : Final_xPrimalSol) {
        VarArray.add(primal_node_vars[i.first]);
        if (i.second) {
            NumArray.add(1);
        } else {
            NumArray.add(0);
        }
    }
    FLBcplex.addMIPStart(VarArray, NumArray);

    /*IloExpr sigma_vars(nenv);
    for (auto i : G->nodes()) {
            sigma_vars += G->nodes_value().at(i)*primal_node_vars[i];
    }
    FLBmodel.add(sigma_vars <= Final_Obj);
    */
    VarArray.end();
    NumArray.end();

    // FLBcplex.setParam(IloCplex::Param::Preprocessing::Reduce, 1);
    // FLBcplex.setParam(IloCplex::Param::Preprocessing::Symmetry, 0);
    FLBcplex.setParam(IloCplex::TiLim, 3600);
    if (formulation > 0)
        FLBcplex.setParam(IloCplex::AdvInd, 1);  // start value: 1
    FLBcplex.setParam(IloCplex::EpGap, 1e-09);   // set MIP gap tolerance
    FLBcplex.setParam(IloCplex::Threads, 0);
    FLBcplex.setParam(IloCplex::TreLim, 12288);
    FLBcplex.setParam(IloCplex::EpInt, 1e-06);  // set integrality tolerance
    FLBcplex.setParam(IloCplex::MIPDisplay, MIPDisplayLevel);

    double start_time = FLBcplex.getCplexTime();
    double start_ticks = FLBcplex.getDetTime();

    FLBcplex.solve();

    double elapsed_time = FLBcplex.getCplexTime() - start_time;
    double elapsed_ticks = FLBcplex.getDetTime() - start_ticks;

    TOT_TIME += elapsed_time;
    FINAL_SOLVE_TIME = elapsed_time;

    cout << "Solution status \t= \t" << FLBcplex.getStatus() << endl;
    cout << "Objectvie value \t= \t" << FLBcplex.getObjValue() << endl;
    cout << "Final Elapsed time \t= \t" << FINAL_SOLVE_TIME << endl;
    cout << "LB Elapsed time \t= \t" << TOT_LB_TIME << endl;
    cout << "TOT Elapsed time \t= \t" << TOT_TIME << endl;
    cout << "Use Local Branch Time \t=\t" << LocalBranchTime << endl << endl;

    print_to_file();
    return;
}

void LBSolver::print_to_file() {
    // begin to write the information into the file
    string store = filename;
    string graph_id = "";
    string newfilename = "";
    if (isdigit(store[store.size() - 6]))
        graph_id = store[store.size() - 6] + store[store.size() - 5];
    else
        graph_id = store[store.size() - 5];
    while (store[store.size() - 1] != '\\') {
        newfilename.push_back(*store.rbegin());
        store.pop_back();
    }
    std::reverse(newfilename.begin(), newfilename.end());
    switch (formulation) {
        case SCF: {
            store = store + "1_SCF";
            break;
        }
        case MCF: {
            store = store + "1_MCF";
            break;
        }
        case STEINER: {
            store = store + "1_STEINER";
            break;
        }
        case NS: {
            store = store + "1_NS";
            break;
        }
    }
    if (relax) store = store + "_relax";
    store = store + ".txt";

    //[Gap] [time] [Status] [Value] [Nodes number] [User number]
    ofstream flow(store, ios::app);
    flow.setf(ios::left, ios::adjustfield);
    flow << setw(LSPACING) << newfilename;
    flow << setw(LSPACING) << TOT_TIME;
    flow << setw(LSPACING) << TOT_LB_TIME;
    flow << setw(LSPACING) << FINAL_SOLVE_TIME;
    flow << setw(LSPACING) << LocalBranchTime;
    flow << setw(LSPACING) << Final_gap;
    flow << setw(LSPACING) << Final_Obj;
    flow << setw(LSPACING) << FLBcplex.getMIPRelativeGap();
    flow << setw(LSPACING) << FLBcplex.getObjValue();
    flow << setw(LSPACING) << FLBcplex.getNnodes();
    flow << setw(LSPACING) << FLBcplex.getNcuts(IloCplex::CutUser);
    flow << setw(LSPACING) << FLBcplex.getStatus();
    // flow << setw(LSPACING) << formulation ;
    // flow << setw(SPACING) << callbackOption ;
    // flow << setw(SPACING) << ns_sep_opt ;
    // flow << setw(SPACING) << time_limit ;
    // flow << setw(SPACING) << max_cuts_lazy;
    // flow << setw(SPACING) << tol_lazy;
    // flow << setw(SPACING) << max_cuts_user;
    // flow << setw(SPACING) << tol_user;
    switch (callbackOption) {
        case 0:
            flow << setw(LSPACING) << "NULL";
            break;
        case 1:
            switch (lazy_sep_opt) {
                case 0:
                    flow << setw(LSPACING + 6) << "L(1-m)";
                    flow << "(" << max_cuts_lazy << ", " << tol_lazy
                         << setw(SPACING + 9) << ")";
                    break;
                case 1:
                    flow << setw(LSPACING + 6) << "L(m-m)";
                    flow << "(" << max_cuts_lazy << ", " << tol_lazy
                         << setw(LSPACING + 9) << ")";
                    break;
            }

            break;
        case 2:
            flow << setw(LSPACING) << "U";
            break;
        case 3:
            switch (lazy_sep_opt) {
                case 0:
                    flow << "L(1-m)";
                    break;
                case 1:
                    flow << "L(m-m)";
                    break;
            }

            switch (ns_sep_opt) {
                case 1:
                    flow << setw(LSPACING) << "U(SCC-MinCut)";
                    break;
                case 0:
                    flow << setw(LSPACING) << "U(MinCut)";
                    break;
            }

            flow << "(" << max_cuts_lazy << ", " << tol_lazy << ";"
                 << max_cuts_user << "," << tol_user << setw(SPACING) << ")";

            break;
        default:
            break;
    }

    flow << setw(LSPACING) << LB_MaxRestarts;
    flow << setw(LSPACING) << LB_MaxIter;
    flow << setw(LSPACING) << BCTime;

    switch (LB_CP_Option) {
        case 0:
            flow << setw(LSPACING) << "USE AS INITIAL SOLUTION";
            break;
        case 1:
            flow << setw(LSPACING) << "USE AS USER CUT";
            break;
        default:
            break;
    }

    flow << endl;
}
