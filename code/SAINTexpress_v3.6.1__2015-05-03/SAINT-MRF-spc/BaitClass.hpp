/*
 * BaitClass.hpp
 *
 *  Created on: June 18, 2012
 *      Author: hwchoi
 */


#ifndef BAITCLASS_HPP_
#define BAITCLASS_HPP_

#include <iostream>
#include <string>
#include <map>
#include <deque>


using namespace std;


class BaitClass {

private:

	int colId;
	string ipId;
	string baitId;
	bool isCtrl;

public:

	void print() const;

	int get_colId() const;
	string get_ipId() const;
	const string & get_baitId() const;
	bool get_isCtrl() const;

	void set_colId( int c );
	void set_ipId( string str );
	void set_baitId( string str );
	void set_isCtrl( bool ic );

};


#endif /* BAITCLASS_HPP_ */
