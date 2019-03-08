/*
 * UIClass.hpp: unique interaction in non-control purifications
 *
 *  Created on: June 18, 2012
 *      Author: hwchoi
 */

#ifndef UICLASS_HPP_
#define UICLASS_HPP_

#include <iostream>
#include <string>
#include <map>
#include <deque>

using namespace std;

class UIClass {

private:

	string baitId;
	string preyId;
	string preyGeneId;
	int rowId;
	// deque<int> colId;
	vector<int> colId;
	double score;
	// potentially other information such as totalQSum, ctrlQSum, number of IPs

public:

	UIClass();
	void print() const;

	string get_baitId() const;
	string get_preyId() const;
	string get_preyGeneId() const;
	int get_rowId() const;
	// deque<int> get_colId() const;
	vector<int> get_colId() const;
	double get_score() const;

	void set_baitId(string str);
	void set_preyId(string str);
	void set_preyGeneId(string str);
	void set_rowId(int r);
	void add_colId(int c);
	void set_score(double x);
};

#endif /* UICLASS_HPP_ */
