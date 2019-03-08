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

#include "InterClass.hpp"

using namespace std;

string InterClass::get_ipId() const {
	return ipId;
}

string InterClass::get_baitId() const {
	return baitId;
}

string InterClass::get_preyId() const {
	return preyId;
}

double InterClass::get_quant() const {
	return quant;
}

int InterClass::get_rowId() const {
	return rowId;
}

int InterClass::get_colId() const {
	return colId;
}

void InterClass::set_ipId(string str) {
	ipId = str;
}

void InterClass::set_baitId(string str) {
	baitId = str;
}

void InterClass::set_preyId(string str) {
	preyId = str;
}

void InterClass::set_quant(double x) {
	quant = x;
}

void InterClass::set_rowId(int r) {
	rowId = r;
}

void InterClass::set_colId(int c) {
	colId = c;
}

void InterClass::print() const {
	cout << ipId << "\t";
	cout << baitId << "\t";
	cout << preyId << "\t";
	cout << quant << endl;
}
