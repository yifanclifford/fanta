/*
 $Id: gspan.cpp,v 1.8 2004/05/21 09:27:17 taku-ku Exp $;

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
#include <iterator>
#include <stdlib.h>
#include <unistd.h>
#include <bitset>
//#include <mex.h>
#include <string.h>
#include <algorithm>
#include "fanta.h"

//yifan
#define EDGE_MAX 700000
#define MIN(a,b) (a < b) ? a:b;
#define PRINT(x,s) for(int i=0;i<s;i++) std::cout<<x[i]<<" ";\
std::cout<<std::endl;
namespace FANTA {

Fanta::Fanta(void) {

}

std::map<unsigned int, unsigned int> Fanta::support_counts(
		Projected &projected) {
	std::map<unsigned int, unsigned int> counts;

	for (Projected::iterator cur = projected.begin(); cur != projected.end();
			++cur) {
		counts[cur->id] += 1;
	}

	return (counts);
}

//yifan
void Fanta::run(unsigned int _pruning, bool _share, unsigned int _definition,
		unsigned int _maxpat, unsigned int _maxnode, double _epsilon,
		double _delta, double _sigma, unsigned int _phi, std::istream &_is,
		std::ostream &_os) {
	pruning = _pruning;
	share = _share;
	maxpat = _maxpat;
	maxnode = _maxnode;
	epsilon = _epsilon;
	delta = _delta;
	sigma = _sigma;
	phi = _phi;
	is = &_is;
	os = &_os;
	mex = false;
	ID = 0;
	sample_counter = 0;
	compute_watch.clear();
	//print
	printf("parameters:\n");
	printf("sigma = %.2f\n", sigma);
	printf("epsilon = %.2f\n", epsilon);
	printf("delta = %.2f\n", delta);
	printf("phi = %.2f\n", phi);

	switch (pruning) {
	case 0:
		printf("no ");
		break;
	case 1:
		printf("check ");
		break;
	case 2:
		printf("eager ");
		break;
	default:
		break;
	}
	printf("pruning\n");
	if (share)
		printf("share computation\n");
	//read
	printf("read the graph\n");
	G.read(*is);
//        G.write(std::cout);
	//compute
	switch (_definition) {
	case 0:        //exact
		std::cout << "exact semantic" << std::endl;
		run_exact();
		break;
	case 1:        //expected
		std::cout << "expected semantic" << std::endl;
		run_expected();
		break;
	case 2:        //probabilistic
		std::cout << "probabilistic semantic" << std::endl;
		run_probabilistic();
		break;
	default:
		break;
	}
	printf("sample=%lld\n", sample_counter);
	printf("ID=%d\n", ID);
	printf("compute time=%f\n", compute_watch.get_time());
	printf("edge time=%f\n", edge_watch.get_time());
	printf("sup time=%f\n", sup_watch.get_time());
}

class MyEdge {
public:
	double p;
	int id;
	bool operator <(const MyEdge& other) const {
		double res = p - other.p;
		if (res > 0)
			return true;
		else
			return false;
	}
};

void Fanta::get_embeddings(Projected &projected) {
	std::map<int, int> m;
	embeddings.clear();
	edge_map.clear();
	int i = 0;
	int id;
//        std::vector<MyEdge> edges;
	for (Projected::iterator cur = projected.begin(); cur != projected.end();
			++cur) {
		Embedding embedding;
		PDFS* p = &(*cur);
		embedding.vertices.push_back(p->edge->from);
		for (; p; p = p->prev) {
			embedding.vertices.push_back(p->edge->to);
			if (m.find(p->edge->id) == m.end()) {
				id = i;
				m[p->edge->id] = id;
//                    edge_map[id] = p->edge->id;
				++i;
			} else
				id = m[p->edge->id];
			embedding.edges.push_back(id);
		}
		embeddings.push_back(embedding);
	}
	id = 0;
	for (std::map<int, int>::iterator entry = m.begin(); entry != m.end();
			++entry) {
//            printf("first=%d\n",i->first);
		edge_map[id++] = entry->first;
	}
//        std::sort(edges.begin(), edges.end(), Fanta::compare);
//        std::sort(edges.begin(), edges.end());
//        for(auto& edge:edges)
//            edge_map[ m[edge.id] ] = edge.id;
}

void Fanta::edge_exact(Projected_map3 &root) {
	for (Projected_iterator3 fromlabel = root.begin(); fromlabel != root.end();
			++fromlabel)
		for (Projected_iterator2 elabel = fromlabel->second.begin();
				elabel != fromlabel->second.end(); ++elabel)
			for (Projected_iterator1 tolabel = elabel->second.begin();
					tolabel != elabel->second.end();) {
				get_embeddings(tolabel->second);
				endpoint = support_exact();
				if (endpoint < sigma) {
					G.remove(fromlabel->first, elabel->first, tolabel->first);
					tolabel = root[fromlabel->first][elabel->first].erase(
							tolabel);
					continue;
				}
				printf("(%d,%d,%d): sup=%d\n", fromlabel->first, elabel->first,
						tolabel->first, endpoint);
				++tolabel;
			}

}

void Fanta::edge_expected(Projected_map3 &root) {
	for (Projected_iterator3 fromlabel = root.begin(); fromlabel != root.end();
			++fromlabel)
		for (Projected_iterator2 elabel = fromlabel->second.begin();
				elabel != fromlabel->second.end(); ++elabel)
			for (Projected_iterator1 tolabel = elabel->second.begin();
					tolabel != elabel->second.end();) {
				printf("edge (%d,%d,%d)", fromlabel->first, elabel->first,
						tolabel->first);
//                    index_watch.resume();
				get_embeddings(tolabel->second);
//                    index_watch.pause();
				endpoint = support_exact();
//                    printf(" endpoint=%d ",endpoint);
				if (endpoint < sigma) {
					printf(" is infrequent\n");
					G.remove(fromlabel->first, elabel->first, tolabel->first);
					tolabel = root[fromlabel->first][elabel->first].erase(
							tolabel);
					continue;
				}
				std::vector<double>&pro =
						edge_pros[fromlabel->first][elabel->first][tolabel->first];
				std::vector<double> spro;
				for (int i = 0; i < endpoint; ++i)
					spro.push_back(1.0);
				compute_watch.resume();
				bool frequent = support_expected(pro, spro);
				compute_watch.pause();
				if (!frequent) {
//                        printf(" [%.2f,%.2f] ", lower, upper);
					printf(" is infrequent\n");
					G.remove(fromlabel->first, elabel->first, tolabel->first);
					tolabel = root[fromlabel->first][elabel->first].erase(
							tolabel);
					continue;
				}
				printf(" [%.2f,%.2f]\n", lower, upper);
//                    *os << "edge ("<<fromlabel->first<<","<<elabel->first<<","<<tolabel->first<<"): ["<<lower<<","<<upper<<"]"<<std::endl;
				++tolabel;
			}

	double error;
	for (std::map<int, std::map<int, std::map<int, std::vector<double> > > >::iterator fl =
			edge_pros.begin(); fl != edge_pros.end(); ++fl)
		for (std::map<int, std::map<int, std::vector<double> > >::iterator el =
				fl->second.begin(); el != fl->second.end(); ++el)
			for (std::map<int, std::vector<double> >::iterator tl =
					el->second.begin(); tl != el->second.end(); ++tl) {
				error = sigma * epsilon / (2 * tl->second.size());
				for (auto&p : tl->second)
					if (p + error < 1)
						p += error;
					else
						p = 1;
			}
}

void Fanta::edge_probabilistic(Projected_map3 &root) {
	double p;
	for (Projected_iterator3 fromlabel = root.begin(); fromlabel != root.end();
			++fromlabel)
		for (Projected_iterator2 elabel = fromlabel->second.begin();
				elabel != fromlabel->second.end(); ++elabel)
			for (Projected_iterator1 tolabel = elabel->second.begin();
					tolabel != elabel->second.end();) {
				get_embeddings(tolabel->second);
//                    endpoint = support_exact();
//                    if (endpoint < sigma)
//                    {
//                        G.remove(fromlabel->first, elabel->first, tolabel->first);
//                        tolabel = root[fromlabel->first][elabel->first].erase(tolabel);
//                        continue;
//                    }
//                    double& p = edge_pro[fromlabel->first][elabel->first][tolabel->first];
				bool frequent = support_probabilistic(p);
				if (!frequent) {
					G.remove(fromlabel->first, elabel->first, tolabel->first);
					tolabel = root[fromlabel->first][elabel->first].erase(
							tolabel);
					continue;
				}
				printf("(%d,%d,%d) [%.2f,%.2f]\n", fromlabel->first,
						elabel->first, tolabel->first, lower, upper);
				edge_pro[fromlabel->first][elabel->first][tolabel->first] = p;
				++tolabel;
			}
//        double error;
//        for (std::map < int,std::map < int, std::map<int,double > > > ::iterator fl = edge_pro.begin(); fl != edge_pro.end(); ++fl)
//            for (std::map < int, std::map<int,double > > ::iterator el = fl->second.begin(); el != fl->second.end(); ++el)
//                for (std::map<int,double >::iterator tl = el->second.begin(); tl != el->second.end(); ++tl)
//                {
//                    error = sigma * epsilon / 2;
//                    auto&p = tl->second;
//                    p += error;
//                    if(p>1) p=1;
//                }
}

void Fanta::project_exact(FANTA::Projected &projected) {
	if (DFS_CODE.size() > 1) {
		printf("pattern ");
		DFS_CODE.write(std::cout);
		get_embeddings(projected);
		endpoint = support_exact();
		if (endpoint < sigma) {
			printf(" is infrequent\n");
			return;
		}
		if (is_min() == false) {
			printf(" is not min\n");
			return;
		}
		printf(" sup=%d\n", endpoint);
	}
	report();

	Projected_map3 forward;
	Projected_map2 backward;

	int maxtoc = enumerate_extensions(projected, backward, forward);
	/* Test all extended substructures. */
	// backward
	for (Projected_iterator2 to = backward.begin(); to != backward.end(); ++to)
		for (Projected_iterator1 elabel = to->second.begin();
				elabel != to->second.end(); ++elabel) {
			DFS_CODE.push(maxtoc, to->first, -1, elabel->first, -1);
			project_exact(elabel->second);
			DFS_CODE.pop();
		}
	// forward
	for (Projected_riterator3 from = forward.rbegin(); from != forward.rend();
			++from)
		for (Projected_iterator2 elabel = from->second.begin();
				elabel != from->second.end(); ++elabel)
			for (Projected_iterator1 tolabel = elabel->second.begin();
					tolabel != elabel->second.end(); ++tolabel) {
				DFS_CODE.push(from->first, maxtoc + 1, -1, elabel->first,
						tolabel->first);
				project_exact(tolabel->second);
				DFS_CODE.pop();
			}
}

void Fanta::project_expected(FANTA::Projected &projected,
		std::vector<double> &spro) {
//        printf("size=%d\n",maxpat);
	if (!DFS_CODE.check_label(maxnode) || DFS_CODE.size() > maxpat)
		return;

	std::vector<double> pro;
	if (DFS_CODE.size() > 1) {
		std::cout << "pattern ";
		DFS_CODE.write(std::cout);
		get_embeddings(projected);
//            printf("get embedding finishes\n");
		endpoint = support_exact();
		if (endpoint < sigma) {
			std::cout << " is infrequent" << std::endl;
			return;
		}
		if (is_min() == false) {
			std::cout << " is not min" << std::endl;
			return;
		}
		bool frequent = support_expected(pro, spro);
//            PRINT(pro, pro.size());
		if (!frequent) {
			std::cout << " is infrequent" << std::endl;
			return;
		}
		printf("[%.2f,%.2f]\n", lower, upper);
	} else {        //the probabilities of basic edges are already calculated
		int fromlabel = DFS_CODE[0].fromlabel;
		int elabel = DFS_CODE[0].elabel;
		int tolabel = DFS_CODE[0].tolabel;
		pro = edge_pros[fromlabel][elabel][tolabel];
		double sup = 0;
		for (double p : pro)
			sup += p;
		lower = sup - sigma * epsilon / 2;
		upper = sup + sigma * epsilon / 2;
	}
	report();

	Projected_map2 backward;
	Projected_map3 forward;

	/* Test all extended substructures. */
	int maxtoc = enumerate_extensions(projected, backward, forward);

	double error = sigma * epsilon / (2 * stoppoint);
	for (auto&p : pro)
		if (p + error < 1)
			p += error;
		else
			p = 1;

	// backward
	for (Projected_iterator2 to = backward.begin(); to != backward.end(); ++to)
		for (Projected_iterator1 elabel = to->second.begin();
				elabel != to->second.end(); ++elabel) {
			DFS_CODE.push(maxtoc, to->first, -1, elabel->first, -1);
			project_expected(elabel->second, pro);
			DFS_CODE.pop();
		}
	// forward
	for (Projected_riterator3 from = forward.rbegin(); from != forward.rend();
			++from)
		for (Projected_iterator2 elabel = from->second.begin();
				elabel != from->second.end(); ++elabel)
			for (Projected_iterator1 tolabel = elabel->second.begin();
					tolabel != elabel->second.end(); ++tolabel) {
				DFS_CODE.push(from->first, maxtoc + 1, -1, elabel->first,
						tolabel->first);
				project_expected(tolabel->second, pro);
				DFS_CODE.pop();
			}
}

void Fanta::project_probabilistic(FANTA::Projected &projected) {
	if (!DFS_CODE.check_label(maxnode) || DFS_CODE.size() > maxpat)
		return;
	double p;
	if (DFS_CODE.size() > 1) {
		std::cout << "pattern ";
		DFS_CODE.write(std::cout);
		if (is_min() == false) {
			std::cout << " is not min" << std::endl;
			return;
		}
		get_embeddings(projected);
		bool frequent = support_probabilistic(p);
		//            PRINT(pro, pro.size());
		printf("[%.2f,%.2f]\n", lower, upper);
		if (!frequent)
			return;

		double error = sigma * epsilon / 2;
		p += error;
		if (p > 1)
			p = 1;
	} else {        //the probabilities of basic edges are already calculated
		int fromlabel = DFS_CODE[0].fromlabel;
		int elabel = DFS_CODE[0].elabel;
		int tolabel = DFS_CODE[0].tolabel;
		p = edge_pro[fromlabel][elabel][tolabel];
		double error = epsilon * sigma / 2;
		lower = (p - error > 0) ? p - error : 0;
		upper = (p + error < 1) ? p + error : 1;
	}

	report();

	Projected_map2 backward;
	Projected_map3 forward;

	int maxtoc = enumerate_extensions(projected, backward, forward);
	/* Test all extended substructures. */
	// backward
	for (Projected_iterator2 to = backward.begin(); to != backward.end(); ++to)
		for (Projected_iterator1 elabel = to->second.begin();
				elabel != to->second.end(); ++elabel) {
			DFS_CODE.push(maxtoc, to->first, -1, elabel->first, -1);
			project_probabilistic(elabel->second);
			DFS_CODE.pop();
		}
	// forward
	for (Projected_riterator3 from = forward.rbegin(); from != forward.rend();
			++from)
		for (Projected_iterator2 elabel = from->second.begin();
				elabel != from->second.end(); ++elabel)
			for (Projected_iterator1 tolabel = elabel->second.begin();
					tolabel != elabel->second.end(); ++tolabel) {
				DFS_CODE.push(from->first, maxtoc + 1, -1, elabel->first,
						tolabel->first);
				project_probabilistic(tolabel->second);
				DFS_CODE.pop();
			}
}

int Fanta::support_exact() {
	int len = (int) embeddings[0].vertices.size();
	ProjectSet pset(len);
	for (auto& embedding : embeddings)
		pset += embedding.vertices;
	return pset.support();
}

double Fanta::compute_exact() {
	int l = (int) edge_map.size();
	long long n = 1 << l;
	int maxtoc = (int) embeddings[0].vertices.size();
	ProjectSet pset(maxtoc + 1);
	double psup = 0;
	for (long long i = 0; i < n; ++i) {
		long long edge = i;
		double p = 1;
		for (int id = 0; id < l; ++id) {
			if (edge % 2 == 1)
				p *= G.edges[edge_map[id]].probability;
			else
				p *= 1 - G.edges[edge_map[id]].probability;
			edge /= 2;
		}
		std::bitset<EDGE_MAX> graph(i);
		pset.clear();
		for (auto& embedding : embeddings) {
			bool contain = true;
			for (int id : embedding.edges)
				if (!graph[id]) {
					contain = false;
					break;
				}
			if (contain)
				pset += embedding.vertices;
		}
		int sup = pset.support();
		if (sup >= phi)
			psup += p;
	}
	return psup;
}

bool Fanta::support_expected(std::vector<double> &pro,
		std::vector<double> &spro) {
	int N, N1 = 0, n;
	double v = 2 * log(2 / delta) / (epsilon * epsilon * sigma * sigma);
	I = false;
	N = ceil(v * endpoint * endpoint);
//        printf("v=%f\n", v);
//        printf("N=%d\n",N);

	int*X = new int[endpoint];
	memset(X, 0, sizeof(int) * endpoint);
	int res = 0;
	int x = endpoint;
	if (pruning > 0)
		x = 1;
	if (pruning == 2) {
		int fl, el, tl;
		I = preprune(fl, el, tl);
		if (I)
			ps0 = edge_pros[fl][el][tl];
	}

	for (; x <= endpoint; ++x) {
		N = ceil(v * x * x);
		n = N - N1;
		N1 = N;
		pro.clear();
//            if(share)
//                share_expected(X, n);
//            else
		compute_expected(X, n);
		for (int k = 0; k < endpoint; ++k)
			pro.push_back(X[k] * 1.0 / N);
		switch (pruning) {
		case 0:
			res = evaluate_simple(pro, x);
			break;
		case 1:
			res = evaluate_check(pro, x);
			break;
		case 2:
			res = evaluate_prune(pro, spro, x);
		default:
			break;
		}
		if (res > 0)
			break;
	}
	delete[] X;
	sample_counter += N1;
//        printf("N1=%d\n",N1);
	if (res == 1)
		return true;
	else
		return false;
}

bool Fanta::support_probabilistic(double &p) {
//        printf("embedding size=%ld\n",embeddings.size());

	int N = ceil(2 * log(2 / delta) / (epsilon * epsilon));
	sample_counter += N;
	double error = epsilon * sigma / 2;
//        printf("t=%f\n",compute_watch.get_time());

	compute_watch.resume();
	if (share)
		p = share_probabilistic(N);
	else
		p = compute_probabilistic(N);
	compute_watch.pause();
	lower = p - error;
	if (lower < 0)
		lower = 0;
	upper = p + error;
	if (upper > 1)
		upper = 1;
	if (upper < sigma)
		return false;
	return true;
}

int Fanta::evaluate_simple(std::vector<double> &pro, int x) {
	double sup = 0;
	if (x == endpoint) {
		for (int i = 0; i < endpoint; ++i)
			sup += pro[i];
		double error = epsilon * sigma / 2;
		lower = sup - error;
		upper = sup + error;
		if (upper < sigma)
			return 2;
		else
			return 1;
	}
	return 0;
}

int Fanta::evaluate_check(std::vector<double> &pro, int x) {
	lower = upper = 0;
	double error = epsilon * sigma / (2 * x);
	double a;
	for (int i = 0; i < endpoint; ++i) {
		a = pro[i] - error;
		if (a > 0)
			lower += a;
		a = pro[i] + error;
		if (a > 1)
			a = 1;
		upper += a;
	}
	if (upper < sigma)
		return 2;
	if (lower > sigma * (1 - epsilon))
		return 1;
	return 0;
}

int Fanta::evaluate_prune(std::vector<double> &pro, std::vector<double> &spro,
		int x) {
	lower = upper = 0;
	double error = epsilon * sigma / (2 * x);
	double a, b;
	for (int i = 0; i < endpoint; ++i) {
		a = spro[i];
		if (I)
			a *= ps0[i];
		b = pro[i] + error;
		upper += MIN(a, b)
		;
		if (pro[i] > error)
			lower += pro[i] - error;
	}
	if (upper < sigma)
		return 2;              //infrequent
	if (lower > sigma * (1 - epsilon))
		return 1;
	return 0;
}

bool Fanta::preprune(int&fromlabel0, int&elabel0, int&tolabel0) {
	if (DFS_CODE.size() == 0)
		return false; //start from null
	int maxtoc = (int) embeddings[0].vertices.size();
	int* labels = new int[maxtoc + 1];
	for (DFS& edge : DFS_CODE) {
		if (edge.fromlabel != -1)
			labels[edge.from] = edge.fromlabel;
		if (edge.tolabel != -1)
			labels[edge.to] = edge.tolabel;
	}
	auto& e0 = DFS_CODE[DFS_CODE.size() - 1];
	//find the edge
	fromlabel0 = labels[e0.from];
	elabel0 = e0.elabel;
	tolabel0 = labels[e0.to];
	if (fromlabel0 > tolabel0)
		_swap(fromlabel0, tolabel0);
	bool I = true;
	for (int i = 0; i < DFS_CODE.size() - 1; i++) {
		auto& e = DFS_CODE[i];
		int fromlabel = labels[e.from];
		int tolabel = labels[e.to];
		if (fromlabel > tolabel)
			_swap(fromlabel, tolabel);
		if (fromlabel0 == fromlabel && e0.elabel == e.elabel
				&& tolabel0 == tolabel) {
			I = false;
			break;
		}
	}
	delete[] labels;
	return I;
}

double Fanta::share_probabilistic(int n) {
	int l1 = (int) edge_map.size();
//        printf("embedding edge size=%d\n",l1);
	int l2 = (int) embeddings.size();
	bool* edges = new bool[l1];
	bool* projects = new bool[l2];
	memset(edges, 0, sizeof(bool) * l1);
	memset(projects, 0, sizeof(bool) * l2);
	int maxtoc = (int) embeddings[0].vertices.size();
	ProjectSet pset(maxtoc);
	int x = 0;
//        long long old = share_time;
	share_dfs(0, n, x, edges, projects, pset);
//        printf("share time=%lld\n",share_time-old);
	delete[] edges;
	delete[] projects;
	return x * 1.0 / n;
}

void Fanta::share_dfs(int cur, int n, int &x, bool *edges, bool *projects,
		ProjectSet &pset) {
	int l = (int) edge_map.size();
	if (n == 1 || cur == l) {
		double r;
		std::unordered_set<int> contains;
		bool contain1, contain2;
		edge_watch.resume();
		for (int e = cur; e < l; ++e) {
			r = rand() * 1.0 / RAND_MAX;
			if (r < G.edges[edge_map[e]].probability)
				edges[e] = true;
			else
				edges[e] = false;
		}
		edge_watch.pause();
		l = (int) embeddings.size();
		int cost1 = 0;
		int cost2 = 0;
		bool* contain = new bool[l];
		memset(contain, 0, sizeof(bool) * l);
		sup_watch.resume();
		for (int i = 0; i < l; ++i) {
			contain1 = true;
			contain2 = projects[i];
			auto&embedding = embeddings[i];
			for (int id : embedding.edges)
				if (!edges[id]) {
					contain1 = false;
					break;
				}
			contain[i] = contain1;
			if (contain1)
				++cost2;
			if ((contain1 && !contain2) || (!contain1 && contain2))
				++cost1;
		}
		if (cost1 < cost2) {
			share_time += cost2 - cost1;
//                                printf("share\n");
			for (int i = 0; i < l; ++i) {
				auto&embedding = embeddings[i];
				if (contain[i] && !projects[i]) {
					pset += embedding.vertices;
					projects[i] = true;
				}
				if (!contain[i] && projects[i]) {
					pset -= embedding.vertices;
					projects[i] = false;
				}
			}
		} else {
//                                printf("scratch\n");
			pset.clear();
			for (int i = 0; i < l; ++i) {
				if (contain[i]) {
					auto&embedding = embeddings[i];
					pset += embedding.vertices;
					projects[i] = true;
				} else
					projects[i] = false;
			}
		}
		delete[] contain;
		int sup = pset.support();
		if (sup >= phi)
			x += n;
		sup_watch.pause();
	}

	else {
		edge_watch.resume();
		double r;
		int left = 0;
		for (int i = 0; i < n; ++i) {
			r = rand() * 1.0 / RAND_MAX;
			if (r < G.edges[edge_map[cur]].probability)
				++left;
		}
		edge_watch.pause();
		if (left > 0) {
			edges[cur] = true;
			share_dfs(cur + 1, left, x, edges, projects, pset);
		}
		if (n - left > 0) {
			edges[cur] = false;
			share_dfs(cur + 1, n - left, x, edges, projects, pset);
		}
	}

}

void Fanta::share_dfs(int cur, int n, int* X, bool* edges, bool* projects,
		ProjectSet &pset) {
	int l = (int) edge_map.size();
	if (n == 1 || cur == l) {
		double r;
//            std::unordered_set<int> contains;
		bool contain1, contain2;
		for (int e = cur; e < l; ++e) {
			r = rand() * 1.0 / RAND_MAX;
			if (r < G.edges[edge_map[e]].probability)
				edges[e] = true;
			else
				edges[e] = false;
		}
		l = (int) embeddings.size();
		int cost1 = 0;
		int cost2 = 0;
		bool* contain = new bool[l];
		memset(contain, 0, sizeof(bool) * l);
		for (int i = 0; i < l; ++i) {
			contain1 = true;
			contain2 = projects[i];
			auto&embedding = embeddings[i];
			for (int id : embedding.edges)
				if (!edges[id]) {
					contain1 = false;
					break;
				}
			contain[i] = contain1;
			if (contain1)
				++cost2;
			if ((contain1 && !contain2) || (!contain1 && contain2))
				++cost1;
		}

		if (cost1 < cost2) {
//                share_time += cost2-cost1;
//                printf("share\n");
			for (int i = 0; i < l; ++i) {
				auto&embedding = embeddings[i];
				if (contain[i] && !projects[i]) {
					pset += embedding.vertices;
					projects[i] = true;
				}
				if (!contain[i] && projects[i]) {
					pset -= embedding.vertices;
					projects[i] = false;
				}
			}
		} else {
//                printf("scratch\n");
			pset.clear();
			for (int i = 0; i < l; ++i) {
				if (contain[i]) {
					auto&embedding = embeddings[i];
					pset += embedding.vertices;
					projects[i] = true;
				} else
					projects[i] = false;
			}
		}
		delete[] contain;
		int sup = pset.support();
		for (int i = 0; i < sup; ++i)
			X[i] += n;
	} else {
		double r;
		int left = 0;
		for (int i = 0; i < n; ++i) {
			r = rand() * 1.0 / RAND_MAX;
			if (r < G.edges[edge_map[cur]].probability)
				++left;
		}
//            printf("left=%d,right=%d\n",left,n-left);
		if (left > 0) {
			edges[cur] = true;
			share_dfs(cur + 1, left, X, edges, projects, pset);
		}
		if (n - left > 0) {
			edges[cur] = false;
			share_dfs(cur + 1, n - left, X, edges, projects, pset);
		}
	}
}

void Fanta::compute_expected(int*X, int n) {
	int maxtoc = (int) embeddings[0].vertices.size();
	ProjectSet pset(maxtoc);
	double r;
	int l = (int) edge_map.size();
	bool* graph = new bool[l];
	for (int t = 0; t < n; ++t) {
		memset(graph, 0, sizeof(bool) * l);
		for (int e = 0; e < l; ++e) {
			r = rand() * 1.0 / RAND_MAX;
			if (r < G.edges[edge_map[e]].probability)
				graph[e] = 1;
		}
		pset.clear();
		for (auto& embedding : embeddings) {
			bool contain = true;
			for (int id : embedding.edges)
				if (!graph[id]) {
					contain = false;
					break;
				}
			if (contain)
				pset += embedding.vertices;
		}
		int sup = pset.support();
		for (int i = 0; i < sup; ++i)
			++X[i];
	}
	delete[] graph;
}

double Fanta::compute_probabilistic(int n) {
	int count = 0;
	int maxtoc = (int) embeddings[0].vertices.size();
	ProjectSet pset(maxtoc);
	double r;
	int l = (int) edge_map.size();
	bool* graph = new bool[l];
	for (int t = 0; t < n; ++t) {
		memset(graph, 0, sizeof(bool) * l);
		for (int e = 0; e < l; ++e) {
			r = rand() * 1.0 / RAND_MAX;
			if (r < G.edges[edge_map[e]].probability)
				graph[e] = 1;
		}
		pset.clear();
		for (auto& embedding : embeddings) {
			bool contain = true;
			for (int id : embedding.edges)
				if (!graph[id]) {
					contain = false;
					break;
				}
			if (contain)
				pset += embedding.vertices;
		}
		int sup = pset.support();
		if (sup >= phi)
			++count;
	}
	delete[] graph;
	return count * 1.0 / n;
}

void Fanta::report() {
	*os << "t # " << ID++ << " [" << lower << "," << upper << "]" << std::endl;
	//*os << "t # " << ID++ << std::endl;
	Graph g(directed);
	DFS_CODE.toGraph(g);
	g.write_pattern(*os);
}

int Fanta::enumerate_extensions(FANTA::Projected &projected,
		Projected_map2 &backward, Projected_map3 &forward) {
	/* Enumerate all possible one edge extensions of the current substructure.   */
	EdgeList edges;
	const RMPath &rmpath = DFS_CODE.buildRMPath();
	int minlabel = DFS_CODE[0].fromlabel;
	int maxtoc = DFS_CODE[rmpath[0]].to;
	for (unsigned int n = 0; n < projected.size(); ++n) {
		PDFS *cur = &projected[n];
		History history(G, cur);
		// backward
		for (int i = (int) rmpath.size() - 1; i >= 1; --i) {
			Edge *e = get_backward(G, history[rmpath[i]], history[rmpath[0]],
					history);
			if (e)
				backward[DFS_CODE[rmpath[i]].from][e->elabel].push(e, cur);
		}
		// pure forward
		if (get_forward_pure(G, history[rmpath[0]], minlabel, history, edges)) {
			for (EdgeList::iterator it = edges.begin(); it != edges.end(); ++it)
				forward[maxtoc][(*it)->elabel][G[(*it)->to].label].push(*it,
						cur);
		}
		// backtracked forward
		for (int i = 0; i < (int) rmpath.size(); ++i)
			if (get_forward_rmpath(G, history[rmpath[i]], minlabel, history,
					edges))
				for (EdgeList::iterator it = edges.begin(); it != edges.end();
						++it)
					forward[DFS_CODE[rmpath[i]].from][(*it)->elabel][G[(*it)->to].label].push(
							*it, cur);
	}
	return maxtoc;
}

void Fanta::run_exact() {
	EdgeList edges;
	Projected_map3 root;
	for (unsigned int from = 0; from < G.size(); ++from)
		if (get_forward_root(G, G[from], edges))
			for (EdgeList::iterator it = edges.begin(); it != edges.end(); ++it)
				root[G[from].label][(*it)->elabel][G[(*it)->to].label].push(*it,
						0);

	std::cout << "find frequent edges" << std::endl;
	edge_exact(root);
	std::cout << "start mining:" << std::endl;

	for (Projected_iterator3 fromlabel = root.begin(); fromlabel != root.end();
			++fromlabel)
		for (Projected_iterator2 elabel = fromlabel->second.begin();
				elabel != fromlabel->second.end(); ++elabel)
			for (Projected_iterator1 tolabel = elabel->second.begin();
					tolabel != elabel->second.end();) {
				//Build the initial two-node graph.  It will be grown recursively within project.
				DFS_CODE.push(0, 1, fromlabel->first, elabel->first,
						tolabel->first);
				project_exact(tolabel->second);
				DFS_CODE.pop();
				int tl = tolabel->first;
				tolabel = root[fromlabel->first][elabel->first].erase(tolabel);
				G.remove(fromlabel->first, elabel->first, tl);
			}
}

void Fanta::run_expected() {
//        share_time = 0;
	Stopwatch watch;
//        watch.start();
	EdgeList edges;
	Projected_map3 root;
	for (unsigned int from = 0; from < G.size(); ++from)
		if (get_forward_root(G, G[from], edges))
			for (EdgeList::iterator it = edges.begin(); it != edges.end(); ++it)
				root[G[from].label][(*it)->elabel][G[(*it)->to].label].push(*it,
						0);
//        watch.pause();
//        printf("t1: %f\n",watch.get_time());
//        watch.resume();
//        index_watch.start();
	watch.start();
	std::cout << "find frequent edges" << std::endl;
	edge_expected(root);
	watch.pause();
	printf("t1: %f\n", watch.get_time());
	watch.resume();
//        printf("index time: %f\n",index_watch.get_time());
	std::cout << "start mining:" << std::endl;

	for (Projected_iterator3 fromlabel = root.begin(); fromlabel != root.end();
			++fromlabel)
		for (Projected_iterator2 elabel = fromlabel->second.begin();
				elabel != fromlabel->second.end(); ++elabel)
			for (Projected_iterator1 tolabel = elabel->second.begin();
					tolabel != elabel->second.end(); ++tolabel) {
				//Build the initial two-node graph.  It will be grown recursively within project.
				DFS_CODE.push(0, 1, fromlabel->first, elabel->first,
						tolabel->first);
				std::vector<double> empty;
				empty.push_back(1);
				project_expected(tolabel->second, empty);
				DFS_CODE.pop();
//                    int tl = tolabel->first;
//                    tolabel = root[fromlabel->first][elabel->first].erase(tolabel);
//                    G.remove(fromlabel->first, elabel->first, tl);
			}
	watch.stop();
//        printf("share time=%lld\n",share_time);
	printf("t2: %f\n", watch.get_time());
}

void Fanta::run_probabilistic() {
//        share_time = 0;
	Stopwatch watch;
	EdgeList edges;
	Projected_map3 root;
	for (unsigned int from = 0; from < G.size(); ++from)
		if (get_forward_root(G, G[from], edges))
			for (EdgeList::iterator it = edges.begin(); it != edges.end(); ++it)
				root[G[from].label][(*it)->elabel][G[(*it)->to].label].push(*it,
						0);
	watch.start();
	std::cout << "find frequent edges" << std::endl;
	edge_probabilistic(root);
	std::cout << "start mining:" << std::endl;
	watch.pause();
	printf("t1=%.3f\n", watch.get_time());
	watch.resume();
	for (Projected_iterator3 fromlabel = root.begin(); fromlabel != root.end();
			++fromlabel)
		for (Projected_iterator2 elabel = fromlabel->second.begin();
				elabel != fromlabel->second.end(); ++elabel)
			for (Projected_iterator1 tolabel = elabel->second.begin();
					tolabel != elabel->second.end(); ++tolabel) {
				//Build the initial two-node graph.  It will be grown recursively within project.
				DFS_CODE.push(0, 1, fromlabel->first, elabel->first,
						tolabel->first);
				project_probabilistic(tolabel->second);
				DFS_CODE.pop();
//				int tl = tolabel->first;
//				tolabel = root[fromlabel->first][elabel->first].erase(tolabel);
//				G.remove(fromlabel->first, elabel->first, tl);
			}
	watch.pause();
	printf("t2=%.3f\n", watch.get_time());
	watch.stop();
	if (share)
		printf("share=%lld\n", share_time);
}

}
