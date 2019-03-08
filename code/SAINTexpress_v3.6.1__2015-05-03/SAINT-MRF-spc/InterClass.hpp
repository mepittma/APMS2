/*
 * InterClass.hpp
 *
 *  Created on: June 18, 2012
 *      Author: hwchoi
 */


#ifndef INTERCLASS_HPP_
#define INTERCLASS_HPP_

#include <iostream>
#include <string>
#include <map>
#include <deque>



using namespace std;

typedef struct {

	string ipId;
	string baitId;
	string preyId;
	double quant;
	int rowId;
	int colId;

} interStruct;


class InterClass {

private:

	string ipId;
	string baitId;
	string preyId;
	double quant;
	int rowId;
	int colId;
public:

	bool is_ctrl;
	void print() const;

	string get_ipId() const;
	string get_baitId() const;
	string get_preyId() const;
	double get_quant() const;
	int get_rowId() const;
	int get_colId() const;

	void set_ipId(string str);
	void set_baitId(string str);
	void set_preyId(string str);
	void set_quant(double x);
	void set_rowId(int r);
	void set_colId(int c);

};


#endif /* INTERCLASS_HPP_ */
