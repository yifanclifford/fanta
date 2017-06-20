//
//  utils.h
//  fanta
//
//  Created by User on 15/9/29.
//  Copyright © 2015年 User. All rights reserved.
//

#ifndef utils_h
#define utils_h

#endif /* utils_h */
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <iostream>
//yifan
struct Embedding {
	std::vector<int> edges, vertices;
	friend std::ostream & operator<<(std::ostream & os, Embedding e) {
		for (int v : e.vertices)
			os << v << ' ';
		return os;
	}
};

class ProjectSet: public std::vector<std::unordered_map<int, int> > {

public:
	ProjectSet(int len);
	ProjectSet() {
		initial(2);
	}
	~ProjectSet();
public:
	void initial(int len);
//    void operator +=(Embedding& embedding);
	void operator +=(std::vector<int>& vertices);
	void operator -=(std::vector<int>& vertices);
//    void add(int from, int to);
//    void set(int id, int pos);
//    void add(Embedding&embedding, int id);

	void clear();
	int support();
	friend std::ostream & operator<<(std::ostream & os, const ProjectSet es);
};
