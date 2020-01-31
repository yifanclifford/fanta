/*
    $Id: gspan.h,v 1.6 2004/05/21 05:50:13 taku-ku Exp $;
 
   Copyright (C) 2004 Taku Kudo, All rights reserved.
     This is free software with ABSOLUTELY NO WARRANTY.
  
   This program is free software; you can redistribute it and/or modify
     it under the terms of the GNU General Public License as published by
     the Free Software Foundation; either version 2 of the License, or
     (at your option) any later version.
    
   This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU General Public License for more details.
    
   You should have received a copy of the GNU General Public License
     along with this program; if not, write to the Free Software
     Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
     02111-1307, USA
*/
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <set>
#include <algorithm>
#include "stopwatch.h"
#include "utils.h"
//#include <mex.h>

namespace FANTA {

template <class T> inline void _swap (T &x, T &y) { T z = x; x = y; y = z; }

struct Edge {
	int from;
	int to;
	int elabel;
    double probability;
	unsigned int id;
	Edge(): from(0), to(0), elabel(0), id(-1), probability(1.0) {};
    friend std::ostream & operator<<(std::ostream & os, const Edge*e) {
        os << "(" << e->from << "," << e->to << "," << e->id << ")";
        return os;
    }
    bool operator < (const Edge& other) const {
        double a = probability-0.5;
        if(a<0) a *= -1;
        double b = other.probability-0.5;
        if(b<0) b *= -1;
        if(a > b) return true;
        else return false;
    }
};

class Vertex
{
public:
	typedef std::vector<Edge>::iterator edge_iterator;

	int label;
	std::vector<Edge> edge;

	void push (int from, int to, int elabel)
	{
		edge.resize (edge.size()+1);
		edge[edge.size()-1].from = from;
		edge[edge.size()-1].to = to;
		edge[edge.size()-1].elabel = elabel;
		return;
	}
    
    void push(int from, int to, int elabel, double p)
    {
        push(from, to, elabel);
        edge[edge.size()-1].probability = p;
    }
};

class Graph: public std::vector<Vertex> {
private:
	unsigned int edge_size_;
public:
	typedef std::vector<Vertex>::iterator vertex_iterator;
    //yifan
    std::vector<Edge> edges;
	Graph (bool _directed)
	{
		directed = _directed;
	};
	bool directed;

	//  int y; // class label
	unsigned int edge_size ()   { return edge_size_; }
	unsigned int vertex_size () { return (unsigned int)size(); } // wrapper
	void buildEdge ();
	std::istream &read (std::istream &); // read
	std::ostream &write (std::ostream &); // write
    std::ostream &write_pattern(std::ostream &); // write pattern
//	void read (const mxArray* graph);
//	void mexprint (void);
//	mxArray* writemex (void);
	void check (void);
	Graph(): edge_size_(0), directed(false) {};
    
    //yifan
    void remove(int fromlabel, int elabel, int tolabel);
};

class DFS {
public:
	int from;
	int to;
	int fromlabel;
	int elabel;
	int tolabel;
	friend bool operator == (const DFS &d1, const DFS &d2)
	{
		return (d1.from == d2.from && d1.to == d2.to && d1.fromlabel == d2.fromlabel
			&& d1.elabel == d2.elabel && d1.tolabel == d2.tolabel);
	}
	friend bool operator != (const DFS &d1, const DFS &d2) { return (! (d1 == d2)); }
	DFS(): from(0), to(0), fromlabel(0), elabel(0), tolabel(0) {};
};

typedef std::vector<int> RMPath;

struct DFSCode: public std::vector <DFS> {
private:
	RMPath rmpath;
public:
	const RMPath& buildRMPath ();

	/* Convert current DFS code into a graph.
	 */
	bool toGraph (Graph &);

	/* Clear current DFS code and build code from the given graph.
	 */
	void fromGraph (Graph &g);

	/* Return number of nodes in the graph.
	 */
	unsigned int nodeCount (void);
    
    bool check_label(int num);
    
	void push (int from, int to, int fromlabel, int elabel, int tolabel)
	{
		resize (size() + 1);
		DFS &d = (*this)[size()-1];

		d.from = from;
		d.to = to;
		d.fromlabel = fromlabel;
		d.elabel = elabel;
		d.tolabel = tolabel;
	}
	void pop () { resize (size()-1); }
	std::ostream &write (std::ostream &); // write
};

struct PDFS {
	unsigned int id;	// ID of the original input graph
	Edge        *edge;
	PDFS        *prev;
	PDFS(): id(0), edge(0), prev(0) {};
};

class History: public std::vector<Edge*> {
private:
	std::vector<int> edge;
	std::vector<int> vertex;

public:
	bool hasEdge   (unsigned int id) { return (bool)edge[id]; }
	bool hasVertex (unsigned int id) { return (bool)vertex[id]; }
	void build     (Graph &, PDFS *);
	History() {};
	History (Graph& g, PDFS *p) { build (g, p); }

};

class Projected: public std::vector<PDFS> {
public:
	void push (int id, Edge *edge, PDFS *prev)
	{
		resize (size() + 1);
		PDFS &d = (*this)[size()-1];
		d.id = id; d.edge = edge; d.prev = prev;
	}
    
    void push (Edge *edge, PDFS *prev)
    {
        resize(size()+1);
        PDFS &d = (*this)[size()-1];
        d.edge = edge; d.prev = prev;
    }
};

typedef std::vector <Edge*> EdgeList;

bool  get_forward_pure   (Graph&, Edge *,  int,    History&, EdgeList &);
bool  get_forward_rmpath (Graph&, Edge *,  int,    History&, EdgeList &);
bool  get_forward_root   (Graph&, Vertex&, EdgeList &);
Edge *get_backward       (Graph&, Edge *,  Edge *, History&);
    
    
class Fanta {

private:

	typedef std::map<int, std::map <int, std::map <int, Projected> > >           Projected_map3;
	typedef std::map<int, std::map <int, Projected> >                            Projected_map2;
	typedef std::map<int, Projected>                                             Projected_map1;
	typedef std::map<int, std::map <int, std::map <int, Projected> > >::iterator Projected_iterator3;
	typedef std::map<int, std::map <int, Projected> >::iterator                  Projected_iterator2;
	typedef std::map<int, Projected>::iterator                                   Projected_iterator1;
	typedef std::map<int, std::map <int, std::map <int, Projected> > >::reverse_iterator Projected_riterator3;

	DFSCode                     DFS_CODE;
	DFSCode                     DFS_CODE_IS_MIN;
	Graph                       GRAPH_IS_MIN;

	bool mex;	// Shall we output to matlab structures?
	class Pattern {
	public:
		struct Pattern* next;
//		mxArray* graph;
		double lower, upper;
//		std::map<unsigned int, unsigned int> counts;
	};

	unsigned int ID;
	double sigma;
	unsigned int maxnode;	// lower bound on node count
	unsigned int maxpat;	// upper bound on node count
	bool where;
	bool enc;
	bool directed;
	std::ostream* os;

	/* Singular vertex handling stuff
	 * [graph][vertexlabel] = count.
	 */
	std::map<unsigned int, std::map<unsigned int, unsigned int> > singleVertex;
	std::map<unsigned int, unsigned int> singleVertexLabel;



	bool is_min ();
	bool project_is_min (Projected &);

	std::map<unsigned int, unsigned int> support_counts (Projected &projected);
	std::istream &read (std::istream &);

public:
	Fanta (void);
    
    //yifan
private://private variables
    double epsilon, delta, phi;
    double lower, upper;
    double esup;
//    unsigned int maxtoc;
    unsigned int endpoint, stoppoint;
    unsigned int pruning;
    bool share, I;
    std::istream* is;
    Graph G;
    std::vector<Embedding> embeddings;
    std::map<int, std::map<int, std::map<int, std::vector<double> > > > edge_pros; //expected      semantic
    std::map<int, std::map<int, std::map<int, double> > > edge_pro;                //probabilistic semantic
    std::map<int,int> edge_map;                                                    //the embedding edges
    std::vector<double> ps0;
    double p0;
    
public://public functions;
    void run(unsigned int _pruning, bool _share, unsigned int _definition, unsigned int _maxpat, unsigned int _maxnode, double _epsilon, double _delta, double _sigma, unsigned int _phi ,std::istream &_is, std::ostream &_os);
private:
    //run the algorithm
    void run_exact();
    void run_expected();
    void run_probabilistic();
    //frequent edges
    void edge_exact(Projected_map3& root);
    void edge_expected(Projected_map3& root);
    void edge_probabilistic(Projected_map3& root);
    //project
    void project_exact(Projected&);
    void project_expected(Projected&, std::vector<double>&spro);
    void project_probabilistic(Projected&);
    //support computation
    int  support_exact();
    bool support_expected(std::vector<double> &pro, std::vector<double> &spro);
    bool support_probabilistic(double &p);
    //compute n samples
    double compute_exact();
    void   compute_expected(int*X, int n);
    double compute_probabilistic(int n);
    //void   simple_dfs(int cur, int n, int&x, bool*edges, bool*projects, ProjectSet&pset);
    //share compute
    //void share_expected(int*X,int n);
    void share_dfs(int cur, int n, int*X, bool*edges, bool*projects, ProjectSet&pset);
    double share_probabilistic(int n);
    void share_dfs(int cur, int n, int&x, bool*edges, bool*projects, ProjectSet&pset);
    //void share_expected2(int*X,int n);
    //evaluate
//    int evaluate_expected(std::vector<double>&pro, std::vector<double>&spro, int x);
    int evaluate_simple(std::vector<double>&pro, int x);
    int evaluate_check(std::vector<double>&pro, int x);
    int evaluate_prune(std::vector<double>&pro, std::vector<double>&spro, int x);
    int prune_probabilistic(int x);
    //projects to embeddings
    bool preprune(int&fl, int&el, int&tl);           //check if distinct edge
    void get_embeddings(Projected& projected);
    void mexAppendPattern(Graph* g);
    void report();
    int enumerate_extensions(Projected& projected, Projected_map2&backward, Projected_map3&forward);//return maxtoc
    
private://debug
    Stopwatch compute_watch,edge_watch,sup_watch;
//    long long sample_num;
    long long share_time;
    long long sample_counter;
};
};


