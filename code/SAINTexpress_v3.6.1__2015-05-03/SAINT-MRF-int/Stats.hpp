/*
 * globals.hpp
 *
 *  Created on: June 20, 2012
 *      Author: hwchoi
 */

#ifndef STATS_HPP_
#define STATS_HPP_

#include <iostream>
#include <vector>
#include <string>

#include "PreyClass.hpp"
#include "BaitClass.hpp"
#include "InterClass.hpp"
#include "UIClass.hpp"
#include "globals.hpp"
#include "Fastmat.hpp"
#include "Stats_fns.hpp"
using namespace std;

Model_data statModel(const vector<vector<size_t> > &p2p_mapping, const vector<string>& ubait, const Quantmat& test_mat_DATA, const Quantmat &ctrl_mat_DATA, const vector<size_t>& ip_idx_to_bait_no, const size_t nprey, const size_t nbait);

#endif /* STATS_HPP_ */
