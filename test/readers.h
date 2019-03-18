#pragma once
//asdad
#include "type.h"
#include "graph.h"

#include <cstdio>
#include <iostream>
#include <fstream>
#include <memory>
#include <map>
#include <string>
#include <sstream>

using namespace std;

class Reader {
public:
	static bool read_graph(std::string filename, std::shared_ptr<Graph>G) {
		ifstream infile(filename);
		if (!infile.good()) {
			cout << "Can't read file: " << filename << endl;
			return false;
		}

		size_t lastdot = filename.find_last_of(".");
		string suffix = filename.substr(lastdot + 1, string::npos);

		if (suffix == "txt")
			return read_graphfile(infile, G);

		return false;
	};

private:
	static bool read_graphfile(std::ifstream& infile, std::shared_ptr<Graph> G) {
		// READ NODE
		int n; infile >> n;
		for (int i = 1; i <= n; i++) {
			int o; infile >> o;
			double val; infile >> val;
			G->Add_node(o, double(val));
		}
		// READ ARC
		int m; infile >> m;
		for (int i = 1; i <= m; i++) {
			int u, v; infile >> u >> v;
			G->Add_arc(u, v);
		}
		// READ PARTATION
		int p; infile >> p;
		for (int i = 1; i <= p; i++) {
			int o; infile >> o;
			G->Add_index_set(o);
			int v_size; infile >> v_size;
			for (int j = 1; j <= v_size; j++) {
				int tmp;
				infile >> tmp;
				G->Add_sub_graphs(o, tmp);
			}
			int t_size; infile >> t_size;
			for (int j = 1; j <= t_size; j++) {
				int tmp;
				infile >> tmp;
				G->Add_terminals(o, tmp);
			}
		}
		// CHECK AND PRINT GRAPH
		G->Generate_sub_graph_P();
		//G->Print_Graph();
		//G->Check_Graph_Logic();
		return true;
	}
};