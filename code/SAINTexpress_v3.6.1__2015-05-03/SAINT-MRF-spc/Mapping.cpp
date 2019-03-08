/*
 * MAPPING.cpp
 *
 *  Created on: June 20, 2012
 *      Author: hwchoi
 */


#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <set>
#include <boost/algorithm/string.hpp>
#include "Mapping.hpp"
#include "globals.hpp"
#include "Stats_fns.hpp"

using namespace std;


vector<string> uniqueBait(const deque<BaitClass> &BDATA, size_t &nbait ) {
	// nbait = 0;

	// for(deque<BaitClass>::const_iterator m = BDATA.begin(); m != BDATA.end(); m++)
	//	if(m->get_isCtrl() == false) {
	//		deque<string>::iterator n;
	//		for(n = ubait.begin(); n != ubait.end(); n++)
	//			if( m->get_baitId() == *n ) break;
	//		if(n == ubait.end()) {
	//			ubait.push_back( m->get_baitId() );
	//			nbait++;
	//		}
	//	}
	set<string> ubait_set;
	vector<string> ubait_vec;
	BOOST_FOREACH(const BaitClass& bdata , BDATA)
		if( !bdata.get_isCtrl() ){
			ubait_vec.push_back(bdata.get_baitId());
			ubait_set.insert(bdata.get_baitId());
		}
	{
		auto tmp_iter = std::unique(ubait_vec.begin(), ubait_vec.end());
		ubait_vec.resize( tmp_iter - ubait_vec.begin() );
	}
	// if(ubait_set.size() != ubait_vec.size())
		// throw runtime_error("non consecutive duplicate elements found");
	nbait = ubait_set.size();
	// return {ubait_set.begin(), ubait_set.end()};
	return vector<string>(ubait_set.begin(), ubait_set.end());
	// return ubait_vec;
}


vector<vector<size_t> > createP2Pmap0(const vector<string>& prey_list, const string input_iRefWeb) {
	ifstream inF(input_iRefWeb);
	if(!inF.is_open()) throw runtime_error("file doesn't exist!");
	const size_t nprey = prey_list.size();
	const string eol = detect_line_ending(input_iRefWeb);
	const char delimiter = '\t';
	// skip first line
	skip_line(inF, eol);
	stringstream ss;
	vector<boost::array<string,2> > iRefWeb_inter_list;
	while(!inF.eof()){
		inF.get(*ss.rdbuf(),eol[0]);
		string p1,p2;
		getline(ss, p1, delimiter);
		getline(ss, p2, delimiter);
		iRefWeb_inter_list.push_back(boost::array<string,2>{{p1,p2}});
		ss.str("");ss.clear();
		skip_line(inF, eol);
	}
	map<string, unsigned> prey_index;
	for(unsigned i=0; i<prey_list.size(); ++i){
		try{
			prey_index.at(prey_list[i]);//TODO
			// exit(9);
		}catch(std::out_of_range&){}
		prey_index[prey_list[i]] = i;
	}
	Fastmat<unsigned> p2p_matrix(nprey, nprey);
	BOOST_FOREACH(auto& interaction, iRefWeb_inter_list)
		try {
			if(interaction[0]==interaction[1]) continue;
			unsigned idx1 = prey_index.at(interaction[0]);
			unsigned idx2 = prey_index.at(interaction[1]);
			++p2p_matrix(idx1, idx2);
			++p2p_matrix(idx2, idx1);
		} catch (std::out_of_range&){}
	vector<vector<size_t> > p2p_mapping(nprey);
	for(size_t r = 0; r < nprey; ++r)
		for(size_t c = r+1; c < nprey; ++c)
			if( p2p_matrix(r, c) >= 1 )
				p2p_mapping[r].push_back(c);
	return p2p_mapping;
}

vector<vector<size_t> > createP2Pmap(const vector<string>& prey_list,const string inputGO) {
	ifstream inF(inputGO);
	if(!inF.is_open()) throw runtime_error("file doesn't exist!");
	const size_t nprey = prey_list.size();
	const string eol = detect_line_ending(inputGO);
	const char delimiter = '\t';
	// skip first line
	skip_line(inF, eol);
	stringstream ss;
	vector<vector<string> > GO_terms_list;
	if(0)
		while(!inF.eof()){
			// skip first column
			inF.ignore(std::numeric_limits<std::streamsize>::max(), delimiter);
			// read file into stringstream until the next delimiter
			inF.get(*ss.rdbuf(), delimiter);
			vector<string> gene_vector;
			string tmp;
			while(getline(ss, tmp, ' '))
				gene_vector.push_back(tmp);
			GO_terms_list.push_back(gene_vector);
			ss.str("");ss.clear();
			skip_line(inF, eol);
			// cout<<tmp.back();
		}
	while(!inF.eof()){
		stringstream ss_line;
		// skip first column
		inF.ignore(std::numeric_limits<std::streamsize>::max(), delimiter);
		// read a line into ss_line
		inF.get(*ss_line.rdbuf(), eol[0]);
		// read ss_line into stringstream until the next delimiter
		ss_line.get(*ss.rdbuf(), delimiter);
		vector<string> gene_vector;
		string tmp;
		while(getline(ss, tmp, ' '))
			gene_vector.push_back(tmp);
		GO_terms_list.push_back(gene_vector);
		ss.str("");ss.clear();
		skip_line(inF, eol);
	}
	map<string, unsigned> prey_index;
	for(unsigned i=0; i<prey_list.size(); ++i){
		try{
			prey_index.at(prey_list[i]);//TODO
			// exit(9);
		}catch(std::out_of_range&){}
		prey_index[prey_list[i]] = i;
	}
	// prey_index.insert(make_pair(prey_list[i], i));
	Fastmat<unsigned> p2p_matrix(nprey, nprey);
	BOOST_FOREACH(vector<string>& gene_vector, GO_terms_list)
		for(unsigned i=0; i<gene_vector.size(); ++i)
			for(unsigned j=i+1; j<gene_vector.size(); ++j)
				try {
					unsigned idx1 = prey_index.at(gene_vector[i]);
					unsigned idx2 = prey_index.at(gene_vector[j]);
					++p2p_matrix(idx1, idx2);
					++p2p_matrix(idx2, idx1);
				} catch (std::out_of_range&){}
	vector<vector<size_t> > p2p_mapping(nprey);
	for(size_t r = 0; r < nprey; ++r)
		for(size_t c = r+1; c < nprey; ++c)
			if( p2p_matrix(r, c) >= 1 )
				p2p_mapping[r].push_back(c);
	return p2p_mapping;
}

