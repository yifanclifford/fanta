//
//  utils.cpp
//  fanta
//
//  Created by User on 15/9/29.
//  Copyright © 2015年 User. All rights reserved.
//
#include "utils.h"
#include <iostream>
#include <limits.h>


ProjectSet::ProjectSet(int size) {
	initial(size);
}

ProjectSet::~ProjectSet(void) {

}

void ProjectSet::initial(int len) {
	clear();
	for (int i = 0; i < len; i++) {
		std::unordered_map<int, int> m;
		push_back(m);
	}
}

void ProjectSet::operator+=(std::vector<int> &vertices) {
	for (int i = 0; i < size(); ++i)
		++(*this)[i][vertices[i]];
}

void ProjectSet::operator-=(std::vector<int> &vertices) {
	for (int i = 0; i < size(); ++i) {
		--(*this)[i][vertices[i]];
//        if((*this)[i][vertices[i]] <= 0)
//            (*this)[i].erase(vertices[i]);
	}
}

void ProjectSet::clear() {
	for (auto& s : (*this))
		s.clear();
}
//obtain the minimum image value
int ProjectSet::support() {
	int min = INT_MAX;
	for (auto& s : (*this)) {
		for (std::unordered_map<int, int>::iterator i = s.begin(); i != s.end();
				) {
			if (i->second <= 0)
				i = s.erase(i);
			else
				++i;
		}
	}
	for (auto& element : (*this)) {
		int size = (int) element.size();
		min = size < min ? size : min;
	}
	min = (min == INT_MAX) ? 0 : min;
	return min;
}

//std::ostream & operator << (std::ostream & os,const ProjectSet es)
//{
//    int v = 0;
//    for(auto vertices:es)
//    {
//        std::cout<<"vertex "<<v<<": ";
//        v++;
//        for(auto vertex:vertices)
//        {
//            std::cout<<vertex<<" ";
//        }
//        std::cout<<std::endl;
//    }
//    return os;
//}
