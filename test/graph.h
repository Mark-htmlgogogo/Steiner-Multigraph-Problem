/*
  file : graph.h
*/

#pragma once
#include <algorithm>
#include <iostream>
#include <map>
#include <set>
#include <vector>

#include "type.h"

class SUB_Graph {
   public:
    SUB_Graph(){};

    /*  Add Info.  */
    void Add_node(map<INDEX, NODE_SET>& v_set, INDEX p) {
        Nodes.assign(v_set[p].begin(), v_set[p].end());
    }

    void Add_Arc(map<NODE, vector<NODE>>& adj_nodes_list,
                 vector<NODE_PAIR>& arcs, INDEX p,
                 map<INDEX, NODE_SET>& v_set) {
        /* for (auto i : Nodes) {
             for (auto j : adj_nodes_list[i]) {
                 // Adj_nodes[i].push_back(j);
                 if (i != j &&
                     std::find(Nodes.begin(), Nodes.end(), j) != Nodes.end() &&
                     std::find(arcs.begin(), arcs.end(), NODE_PAIR(i, j)) !=
                         arcs.end()) {
                     Adj_nodes[i].push_back(j);
                     Arcs.push_back(NODE_PAIR(i, j));
                 }
             }
         } */

        for (auto a : arcs) {
            int u = a.first, v = a.second;
            if (v_set[p].count(u) && v_set[p].count(v)) {
                Adj_nodes[u].push_back(v);
                Arcs.push_back(NODE_PAIR(u, v));
            }
        }
    }

    void Add_T_terminal(map<INDEX, NODE_SET> t_set, INDEX p) {
        for (auto i : t_set[p]) {
            T_Terminal.insert(i);
        }
    }

    void Add_Node_Value(map<NODE, double> nodes_value, INDEX p) {
        for (auto i : Nodes) {
            Node_Value[i] = nodes_value[i];
        }
    }

    void NodeIsTerminal(map<INDEX, NODE_SET> v_set, map<INDEX, NODE_SET> t_set,
                        INDEX p) {
        for (auto i : v_set[p]) {
            if (t_set[p].find(i) == t_set[p].end()) {
                node_is_terminal[i] = 0;
            } else {
                node_is_terminal[i] = 1;
            }
        }
    }

    void AddNodeAdj(INDEX p) {
        for (auto i : Nodes) {
            Adj_Terminal_nodes[i] = vector<NODE>();
            Adj_General_nodes[i] = vector<NODE>();
            for (auto j : Adj_nodes[i]) {
                if (node_is_terminal[j]) {
                    Adj_Terminal_nodes[i].push_back(j);
                } else {
                    Adj_General_nodes[i].push_back(j);
                }
            }
        }
    }

    const vector<NODE>& nodes() const { return Nodes; }
    const vector<NODE_PAIR>& arcs() const { return Arcs; }
    const NODE_SET& t_set() const { return T_Terminal; }
    const map<NODE, double>& node_value() const { return Node_Value; }
    const map<NODE, vector<NODE>>& adj_nodes_list() const { return Adj_nodes; }
    const map<NODE, vector<NODE>>& AdjTerminalNodes() const {
        return Adj_Terminal_nodes;
    }
    const map<NODE, vector<NODE>>& AdjGeneralNodes() const {
        return Adj_General_nodes;
    }
    const map<NODE, bool>& CheckNodeIsTerminal() const {
        return node_is_terminal;
    }

    /*  Print Info.  */
    void print_graph() {
        cout << "NODEs are: " << Nodes << endl;
        cout << "Arcs are: " << endl;
        for (auto i : Arcs) cout << i << endl;
        cout << "Terminals are: " << T_Terminal << endl;
        cout << "Node_Values are: " << endl;
        for (auto i : Node_Value) {
            cout << i.first << " " << i.second << endl;
        }
        cout << "Terminals are: " << endl;
        for (auto i : node_is_terminal) {
            if (!i.second) continue;
            cout << i.first << " " << i.second << endl;
        }
        cout << "Nodes adjacent node info: " << endl;
        for (auto p : Adj_Terminal_nodes) {
            cout << "For Nodes " << p.first << ": " << endl;
            cout << "Terminal nodes are: ";
            for (auto i : p.second) {
                cout << " " << i;
            }
            cout << endl;
            cout << "General nodes are: ";
            for (auto i : Adj_General_nodes[p.first]) {
                cout << " " << i;
            }
            cout << endl;
        }
        cout << endl;
    }

   private:
    /*  NODEs and ARCs  */
    vector<NODE> Nodes;
    vector<NODE_PAIR> Arcs;
    map<NODE, vector<NODE>> Adj_nodes;
    map<NODE, vector<NODE>> Adj_Terminal_nodes;
    map<NODE, vector<NODE>> Adj_General_nodes;
    /*  Terminals  */
    NODE_SET T_Terminal;
    /*  NODE Value  */
    map<NODE, double> Node_Value;
    map<NODE, bool> node_is_terminal;
};

class Graph {
   public:
    Graph(){};

    bool Has_node(NODE n) {
        return std::find(Nodes.begin(), Nodes.end(), n) != Nodes.end();
    }

    /*  add node to graph */
    void Add_node(NODE n, double val) {
        if (!Has_node(n)) Nodes.push_back(n);
        if (Adj_nodes.count(n) == 0) Adj_nodes[n] = vector<NODE>();
        if (Node_Value.count(n) == 0) Node_Value[n] = val;
    }

    /*  add arc to graph */
    void Add_arc(NODE u, NODE v) {
        Adj_nodes[u].push_back(v);
        Adj_nodes[v].push_back(u);
        Arcs.push_back(NODE_PAIR(u, v));
        Arcs.push_back(NODE_PAIR(v, u));
    }

    /*  add index to graph partation */
    void Add_index_set(INDEX part) {
        P_index_set.insert(part);
        V_sub_graph[part] = NODE_SET();
    }

    /*  add sub_graphs */
    void Add_sub_graphs(INDEX part, NODE n) {
        V_sub_graph[part].insert(n);
        Node_Belong_to_V[n].insert(part);
        V_Total.insert(n);
    }

    /*  add terminal on every sub_graph */
    void Add_terminals(INDEX part, NODE n) {
        T_terminal[part].insert(n);
        Node_Belong_to_T[n].insert(part);
        T_Total.insert(n);
    }

    /*  Interface */
    const vector<NODE>& nodes() const { return Nodes; }
    const vector<NODE_PAIR>& arcs() const { return Arcs; }

    const set<INDEX>& p_set() const { return P_index_set; }            // P
    const NODE_SET& v_total() const { return V_Total; }                // V
    const NODE_SET& t_total() const { return T_Total; }                // T
    const map<INDEX, NODE_SET>& v_set() const { return V_sub_graph; }  // V_k
    const map<INDEX, NODE_SET>& t_set() const { return T_terminal; }   // T_k

    const map<NODE, double>& nodes_value() const {
        return Node_Value;
    }  // node weight
    const map<NODE, NODE_SET>& nodes_of_v() const { return Node_Belong_to_V; }
    const map<NODE, NODE_SET>& nodes_of_t() const { return Node_Belong_to_T; }
    const map<NODE, vector<NODE>>& adj_nodes_list() const { return Adj_nodes; }
    const vector<SUB_Graph>& get_subgraph() const { return sub_g; };

    unsigned nodes_num() { return Nodes.size(); }
    unsigned undi_arcs_num() { return Arcs.size() / 2; }
    unsigned di_arcs_num() { return Arcs.size(); }

    /* addition function */
    void Generate_sub_graph_P() {
        sub_g.push_back(SUB_Graph());  // contain the relative p_index
                                       // indentical
        for (auto i : P_index_set) {
            SUB_Graph subg;
            subg.Add_node(V_sub_graph, i);
            subg.Add_Arc(Adj_nodes, Arcs, i, V_sub_graph);
            subg.Add_T_terminal(T_terminal, i);
            subg.Add_Node_Value(Node_Value, i);
            subg.NodeIsTerminal(V_sub_graph, T_terminal, i);
            subg.AddNodeAdj(i);
            sub_g.push_back(subg);
        }
    }

    bool Check_Graph_Logic() {
        cout << "Check _Graph_Logic" << endl;
        // check amount of V
        NODE_SET v_num;
        for (auto i : P_index_set)
            for (auto j : V_sub_graph[i]) {
                v_num.insert(j);
            }
        if (v_num.size() > Nodes.size()) {
            cout << "Number of v_num wrong: " << v_num.size() << " and "
                 << Nodes.size() << endl;
            return false;
        }
        // check amount of T in every V
        v_num.clear();
        bool flag = true;
        for (auto i : P_index_set) {
            if (T_terminal[i].size() > V_sub_graph[i].size()) {
                flag = false;
                cout << "Terminal number in " << i
                     << "is greater than the partiton nodes number" << endl;
            } else
                cout << "Terminal number in " << i << "is correct" << endl;
        }
        if (!flag) return false;
        // check nodes compared to different partation
        flag = true;
        for (auto i : Node_Belong_to_V) {
            for (auto j : i.second) {
                if (V_sub_graph[j].find(i.first) == V_sub_graph[j].end()) {
                    cout << "NODE " << i.first << "is wrong with partation"
                         << endl;
                    return false;
                }
            }
            cout << "NODE " << i.first << "partation is correct" << endl;
        }
        // check nodes compared to different terminal
        flag = true;
        for (auto i : Node_Belong_to_T) {
            for (auto j : i.second) {
                if (T_terminal[j].find(i.first) == T_terminal[j].end()) {
                    cout << "NODE " << i.first << "is wrong with terminial"
                         << endl;
                    return false;
                }
            }
            cout << "NODE " << i.first << "terminal belonging is correct"
                 << endl;
        }
        cout << endl;
    }

    void Print_Graph() {
        cout << "Print_Graph" << endl;
        cout << "FOR NODEs and ARC Info. " << endl;
        cout << "Nodes number is: " << Nodes.size() << endl;
        for (auto i : Nodes) {
            cout << i << " ";
        }
        cout << endl;
        cout << "Arcs number is: " << Arcs.size() << endl;
        for (auto i : Arcs) {
            cout << i << endl;
        }
        cout << "FOR Partation Info." << endl;
        cout << "Total index number is: " << P_index_set.size() << endl;
        for (auto i : P_index_set) {
            cout << i << " ";
        }
        cout << endl;
        cout << "Total sub_graph number is: " << V_sub_graph.size() << endl;
        cout << "Total nodes number stored in V_Total is: " << V_Total.size()
             << endl;
        for (auto i : V_sub_graph) {
            cout << "For index set: " << i.first << endl;
            for (auto j : i.second) cout << j << " ";
            cout << endl;
        }
        cout << "Total terminial number is: " << T_terminal.size() << endl;
        cout << "Total nodes number stored in T_Total is: " << T_Total.size()
             << endl;
        for (auto i : T_terminal) {
            cout << "For index set: " << i.first << endl;
            for (auto j : i.second) cout << j << " ";
            cout << endl;
        }

        cout << "FOR NODE Info." << endl;
        for (auto i : Nodes) {
            cout << "FOR NODE: " << i << endl;
            cout << "NODE value is: " << Node_Value[i] << endl;
            cout << "This NODE belong to which Vs: ";
            cout << Node_Belong_to_V[i];
            cout << endl;
            cout << "This NODE belong to which Ts: ";
            cout << Node_Belong_to_T[i];
            cout << endl;
            cout << "The Adjency NODEs are: ";
            for (auto j : Adj_nodes[i]) {
                cout << j << " ";
            }
            cout << endl;
            cout << endl;
        }

        cout << "FOR Sub Graph Info. " << endl;
        int index = 0;
        for (auto i : sub_g) {
            cout << "FOR V[" << index++ << "]: " << endl;
            i.print_graph();
        }
    }

   private:
    /*  Node and Arc number */
    vector<NODE> Nodes;
    vector<NODE_PAIR> Arcs;
    /*  Partition information */
    set<INDEX> P_index_set;
    NODE_SET V_Total;
    NODE_SET T_Total;
    map<INDEX, NODE_SET> V_sub_graph;
    map<INDEX, NODE_SET> T_terminal;
    /*  NODE informatino */
    map<NODE, double> Node_Value;
    map<NODE, INDEX_SET> Node_Belong_to_V;  // node belong which V part
    map<NODE, INDEX_SET> Node_Belong_to_T;  // node belong which T part
    map<NODE, vector<NODE>> Adj_nodes;
    /*  Sub graph -> Vi */
    vector<SUB_Graph> sub_g;  // partation of V[i]
};
