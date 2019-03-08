/*
 * baitClass.hpp
 *
 *  Created on: June 18, 2012
 *      Author: hwchoi
 */


#ifndef PREYCLASS_HPP_
#define PREYCLASS_HPP_

#include <iostream>
#include <string>
#include <map>
#include <deque>



using namespace std;



class PreyClass {

private:

	int rowId;
	string preyId;
	double preyLength;
	string preyGeneId;

public:

	PreyClass();
	void print() const;

	int get_rowId() const;
	string get_preyId() const;
	double get_preyLength() const;
	string get_preyGeneId() const;

	void set_rowId(int r);
	void set_preyId(string str);
	void set_preyLength(double x);
	void set_preyGeneId(string str);

};


#endif /* PREYCLASS_HPP_ */
