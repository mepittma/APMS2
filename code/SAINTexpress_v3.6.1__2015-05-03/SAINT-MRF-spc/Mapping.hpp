/*
 * mapping.hpp
 *
 *  Created on: June 20, 2012
 *      Author: hwchoi
 */

#ifndef MAPPING_HPP_
 #define MAPPING_HPP_

#include <iostream>
#include <vector>
#include <string>

#include "PreyClass.hpp"
#include "BaitClass.hpp"
#include "InterClass.hpp"
#include "UIClass.hpp"
#include "globals.hpp"


vector<string> uniqueBait(const deque<BaitClass> &BDATA, size_t &nbait );
vector<vector<size_t> > createB2Bmap(deque<UIClass>& UIDATA, const vector<string>& ubait, size_t nbait, const Countmat& test_mat_DATA, const Countmat& ctrl_mat_DATA);
vector<vector<size_t> > createP2Pmap(const vector<string>& prey_list,const string inputGO);

#endif /* MAPPING_HPP_ */
