/*
  file: type.h
*/

#ifndef __TYPE_H__
#define __TYPE_H__

#include <cstdio>
#include <cstring>
#include <iostream>
#include <map>
#include <set>
#include <utility>
#include <vector>

#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wcast-qual"
#pragma GCC diagnostic ignored "-Wredundant-decls"
#include <ilcplex/ilocplex.h>
#pragma GCC diagnostic pop
#pragma GCC diagnostic pop
#pragma GCC diagnostic pop

using namespace std;

typedef long long LL;

typedef int NODE;
typedef int INDEX;
typedef pair<int, int> NODE_PAIR;
typedef set<NODE> NODE_SET;
typedef set<INDEX> INDEX_SET;

// typedef			pair<NODE_PAIR, INDEX>
//				TRIPLET;

// typedef			map<pair<NODE_PAIR, INDEX>, IloNumVar>
//				EDGE_VAR;

enum SmpForm { NONE = 0, SCF = 1, MCF = 2, STEINER = 3, NS = 4 };

static std::vector<std::string> Multi_Stp_Name = {"NONE", "SCF", "MCF", "SF",
                                                  "NBM"};

template <class T1>
std::ostream& operator<<(std::ostream& out, const std::set<T1>& value) {
    for (auto i : value) out << i << " ";
    out << endl;
    return out;
}

template <class T2>
std::ostream& operator<<(std::ostream& out, const std::vector<T2>& value) {
    for (auto i : value) out << i << " ";
    return out;
}

template <class T1, class T2>
std::ostream& operator<<(std::ostream& out, const std::pair<T1, T2>& value) {
    out << value.first << " to " << value.second;
    return out;
}

/** Hash function for type pair< , > */
namespace std {
template <typename a, typename b>
struct hash<pair<a, b> > {
   private:
    const hash<a> ah;
    const hash<b> bh;

   public:
    hash() : ah(), bh() {}
    size_t operator()(const pair<a, b>& p) const {
        size_t seed = ah(p.first);
        return bh(p.second) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }
};
}  // namespace std

#endif