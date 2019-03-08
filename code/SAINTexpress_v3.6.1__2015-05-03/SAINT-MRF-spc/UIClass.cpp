/*
 * IntrxClass.cpp
 *
 *  Created on: Jun 18, 2012
 *      Author: hwchoi
 */

#include <iostream>
#include <string>
#include <fstream>
#include <vector>

#include "UIClass.hpp"

using namespace std;

UIClass::UIClass( ) {

}

void UIClass::print() const {
	cout << baitId << "\t";
	cout << preyId << "\t";
	cout << rowId << "\t";
	for(auto i = colId.begin(); i != colId.end(); i++) cout << *i << " ";
	cout << endl;
}

string UIClass::get_baitId() const {
	return baitId;
}

string UIClass::get_preyId() const {
	return preyId;
}

string UIClass::get_preyGeneId() const {
	return preyGeneId;
}

int UIClass::get_rowId() const {
	return rowId;
}

vector<int> UIClass::get_colId() const {
	return colId;
}

double UIClass::get_score() const {
	return score;
}

void UIClass::set_baitId(string str) {
	baitId = str;
}

void UIClass::set_preyId(string str) {
	preyId = str;
}

void UIClass::set_preyGeneId(string str) {
	preyGeneId = str;
}

void UIClass::set_rowId(int r) {
	rowId = r;
}

void UIClass::add_colId(int c) {
	colId.push_back( c );
}

void UIClass::set_score(double x) {
	score = x;
}
