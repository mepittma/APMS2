/*
 * Stats.cpp
 *
 *  Created on: June 20, 2012
 *      Author: hwchoi
 */

#include <iostream>
#include <iomanip>
#include <string>
#include <valarray>
#include <vector>
#include <boost/array.hpp>

//#include <array>
#include <algorithm>
#include <numeric>
#include <boost/algorithm/string.hpp>
#include <nlopt.hpp>
#include "Stats.hpp"
#include "globals.hpp"
#include "nlopt_wrap.hpp"
using namespace std;

pair<double, double> mean_sd_estimate(const valarray<Quant_t> data,const Quant_t t, const double sd_NA, const pair<double, double> lse) {
	const size_t detected = data.size() - std::count_if(begin(data),end(data),saint::isnan);
	double sum_detected = 0,sum_squares_detected = 0;
	for(size_t i=0; i<data.size(); i++)
		if(!saint::isnan(data[i])) {
			sum_detected += data[i];
			sum_squares_detected += data[i]*data[i];
		}
	if(detected == 0)
		return make_pair(t - log(1), sd_NA);
	const double mu = sum_detected/detected;
	if(detected >= 3) {
		const double sd = sqrt(sum_squares_detected/detected - mu*mu);
		return make_pair(mu,sd);
	}
	const double sd = exp( mu * lse.second + lse.first );
	return make_pair(mu,sd);
}

// simple linear regression: x predictor, y response
pair<double, double> slr(const vector<double>& x, const vector<double>& y) {
	const double x_bar = mean(x), y_bar = mean(y);
	double beta_numerator = 0, beta_denominator = 0;
	for(size_t i=0; i< x.size(); i++){
		double tmp = x[i] - x_bar;
		beta_numerator += tmp * (y[i]-y_bar);
		beta_denominator += tmp * tmp;
	}
	const double beta = beta_numerator/beta_denominator;
	return make_pair(y_bar - beta * x_bar, beta);
}



Model_data statModel(const vector<vector<size_t> >& p2p_mapping, const vector<string>& ubait, const Quantmat& test_mat_DATA, const Quantmat& ctrl_mat_DATA, const vector<size_t>& ip_idx_to_bait_no, const size_t nprey, const size_t nbait) {

	// bait index to ip index
	vector<vector<size_t> > bait_no_to_ip_idxes(ubait.size());
	for(size_t i = 0; i < ubait.size(); ++i)
		for(size_t j = 0; j < ip_idx_to_bait_no.size(); ++j)
			if( ip_idx_to_bait_no[j] == i )
				bait_no_to_ip_idxes[i].push_back(j);

	// a vector of number of test replicates for each prey (row). All preys have the same vector.
	vector<unsigned char> n_rep_vec(nbait);
	transform(bait_no_to_ip_idxes.begin(), bait_no_to_ip_idxes.end(),
			  n_rep_vec.begin(),
			  // [](const vector<size_t>& a) { return a.size(); });
			  container_size<vector<size_t> >);
	const size_t n_ctrl_ip = ctrl_mat_DATA.n_cols;
	// vector<Quant_t> ts_IP(ip_idx_to_bait_no.size());

	///// calculating the value to impute by averaging over all preys with missing values.
	vector<Quant_t> ctrls_with_NA;
	for(unsigned i=0; i<ctrl_mat_DATA.n_rows;++i){
		bool full=true;
		vector<Quant_t> tmp;
		for(unsigned j=0;j<n_ctrl_ip;++j){
			if(saint::isnan(ctrl_mat_DATA(i,j))) full=false;
			else tmp.push_back(ctrl_mat_DATA(i,j));
		}
		if(!full)
			ctrls_with_NA.insert(ctrls_with_NA.end(),tmp.cbegin(),tmp.cend());
	}
	const Quant_t t = accumulate(ctrls_with_NA.cbegin(),ctrls_with_NA.cend(),0.0)/ctrls_with_NA.size();
	// half of the overall mimimum intensity in log scale
	// Quant_t t = min(test_mat_DATA.min(),ctrl_mat_DATA.min());

	// initial esitmates for η and d
	vector<double> eta(nprey), d(nprey);
	valarray<double> sd_ctrl(nprey);
	extern Options opts;
	int L_user_input = opts.L;
	if(L_user_input > static_cast<int>(n_ctrl_ip))
		cerr << "L is larger than the number of control IPs.\n";
	if(L_user_input <= 0) {
		cerr << "L is must be a positive integer." << endl;
		throw runtime_error("invalid L");
	}
	int L = min<size_t>(n_ctrl_ip, L_user_input);

	// get linear regression estimates for ctrls without any observed, by using the ctrls with full observations
	vector<double> log_sd_ctrls, mean_ctrls;
	for(size_t idx = 0; idx < ctrl_mat_DATA.n_rows; idx++){
		bool full=true;
		for(size_t i=0; i< ctrl_mat_DATA.n_cols; ++i)
			if (saint::isnan(ctrl_mat_DATA(idx, i)))
				full=false;
		if(full) {
			log_sd_ctrls.push_back(log(sqrt(var_MLE((ctrl_mat_DATA(idx,__))))));
			mean_ctrls.push_back(mean(ctrl_mat_DATA(idx,__)));
		}
	}

	// least squares estimate to estimate the sd of ctrls without any observation.
	const auto lse = slr(mean_ctrls, log_sd_ctrls);

	vector<double> sorted_log_sd_ctrls=log_sd_ctrls, sorted_mean_ctrls=mean_ctrls;
	sort(sorted_mean_ctrls.begin(),sorted_mean_ctrls.end());
	sort(sorted_log_sd_ctrls.begin(),sorted_log_sd_ctrls.end());
	double median_sd = exp(sorted_log_sd_ctrls[sorted_log_sd_ctrls.size()/2]);
	double tenpercentile=sorted_mean_ctrls[0];
	// t = tenpercentile;
	// double sd_NA = exp( (t-log(1)) * lse.second + lse.first );
	double sd_NA = median_sd;
	// initialize eta
	// estimate ctrl means and variance
	for(size_t idx = 0; idx < ctrl_mat_DATA.n_rows; idx++){
		valarray<Quant_t> ctrl = ctrl_mat_DATA(idx, __);
		auto tmp = mean_sd_estimate(ctrl, t, sd_NA, lse);
		eta[idx] = tmp.first;
		// sd_ctrl[idx] = tmp.second;
		// sd_ctrl[idx] = max(sd_NA/2, tmp.second);
		sd_ctrl[idx] = max(sd_NA, tmp.second);
	}



	// matrix of valarrays on original data
	Fastmat<valarray<Quant_t> > test_mat(nprey, ubait.size());
	for(size_t i = 0; i < test_mat.n_rows; ++i)
		for(size_t j = 0; j < test_mat.n_cols; ++j){
			//test_mat(i,j) = valarray<Quant_t>(bait_no_to_ip_idxes.at(j).size());
			test_mat(i,j).resize(bait_no_to_ip_idxes.at(j).size());
			for(size_t k = 0; k < bait_no_to_ip_idxes.at(j).size(); ++k) {
				const auto ip_idx = bait_no_to_ip_idxes.at(j).at(k);
				test_mat(i,j)[k] = test_mat_DATA(i,ip_idx);
			}
		}

	// σ for Z=0
	// valarray<double> sd_false(&sd_ctrl[0], sd_ctrl.size());
	valarray<double> sd_false = sd_ctrl;
	// initialize d
	for(size_t i=0; i < nprey; ++i) {
		Quant_t sum = 0;
		size_t n_counts = 0;
		for(size_t ip = 0; ip < test_mat_DATA.n_cols; ++ip)
			if(saint::isnan(test_mat_DATA(i,ip)) ? false : eta[i] < test_mat_DATA(i,ip)) {
				sum += test_mat_DATA(i,ip);
				++n_counts;
			}
		double mu = 0;
		if( n_counts != 0)
			mu = static_cast<double>(sum)/n_counts;
		mu = max(mu, eta[i]);
		d[i] = mu - eta[i];
		if(d[i]<log(4))
			d[i] = log(4);
		// d[i] = max<double>(d[i],sd_false[i] * 2);
	}
	Fastmat<valarray<Quant_t> >& test_mat1 = test_mat;
	// initialize Z
	Fastmat<valarray1<bool> > Z(nprey, ubait.size());
	for(size_t i=0; i < nprey; ++i)
		for(size_t j=0; j < nbait; ++j)
			Z(i,j).resize(n_rep_vec.at(j));

	for(size_t i=0; i < nprey; ++i)
		for(size_t j=0; j < nbait; ++j)
			for(size_t rep = 0; rep < bait_no_to_ip_idxes.at(j).size(); ++rep)
				// Z(i,j)[rep] = saint::isnan(test_mat1(i,j)[rep])? false : (test_mat1(i,j)[rep] > (eta[i] + sd_false[i] * 2));
				Z(i,j)[rep] = saint::isnan(test_mat1(i,j)[rep])? false : (test_mat1(i,j)[rep] > (eta[i] + log(5)));

	// σ for Z=1
	valarray<double> sd_true(0.0,nprey);
	for(size_t i=0; i < nprey; ++i){
		vector<Quant_t> Z1;
		for(size_t j=0; j < nbait; ++j)
			for(size_t rep = 0; rep < bait_no_to_ip_idxes.at(j).size(); ++rep)
				if(Z(i,j)[rep])
					Z1.push_back(test_mat1(i,j)[rep]);
		if(Z1.size() >= 2){
			// BOOST_FOREACH(auto& z, Z1)
			// 	if(saint::isnan(z))
			// 		z=t;
			sd_true[i] = sqrt(var1(Z1));
		}
		else
			sd_true[i] = sd_false[i];
		sd_true[i] = max(sd_NA, sd_true[i]);
		// sd_true[i] = max(sd_false[i], sd_true[i]);
	}
	double beta0 = 0.0, beta1 = 0.0, gamma = 0.;
	return {
		// data
		nprey, nbait,
			(test_mat_DATA),
			(test_mat),
			(test_mat1),
			(n_rep_vec),
			(p2p_mapping),
			t,
			(n_ctrl_ip),
			// parameters
			// latent variable
			(Z),
			// MRF parameters
			beta0, beta1, gamma,
			eta, d, sd_true, sd_false
			};
}

double Model_data::llik_MRF_gamma_0( const std::vector<double>& x, std::vector<double>& /*grad*/) const {
	const double beta0_ = 0;
	double beta1_ = x[0];
	double loglik = 0;
	double MRFtrue = exp(beta1_/* + gamma * gsum*/);
	const double MRFfalse = exp(beta0_);
	double tmpsum=0,countsum=0;
	for(size_t i=0; i < nprey; ++i)
		for(size_t j=0; j < nbait; ++j) {
			const auto& y = test_mat1(i,j);
			if(all_of(begin(y), end(y), saint::isnan)) continue;
			const auto& k = Z(i,j);
			const auto m = n_rep_vec[j];
			// gsum is the covariate for the gamma term
			// double gsum = 0;
			// for(const auto l : b2b_mapping[j])
			// 	gsum += Z(i, l).mean();
			// geometric mean of the likelihood of the replicates
			double prod = 1;
			for(unsigned rep = 0; rep < m; rep++) {
				tmpsum += k[rep]/double(m); countsum+=1./m;
				prod *= (k[rep] ? MRFtrue : MRFfalse)/ (MRFtrue + MRFfalse);
			}
			loglik += log(prod)/m;
		}
	return loglik;
}

double Model_data::llik_MRF( const std::vector<double>& x, std::vector<double>& /*grad*/, const Fastmat<double>& gsum_mat) const {
	const double beta0_ = 0;
	double beta1_ = x[0], gamma_ = x[1];

	const double MRFfalse = exp(beta0_);

	double loglik = 0;
	for(size_t i=0; i < nprey; ++i)
		for(size_t j=0; j < nbait; ++j) {
			const auto& y = test_mat1(i,j);
			if(all_of(begin(y), end(y), saint::isnan)) continue;
			const auto& k = Z(i,j);
			const auto m = n_rep_vec[j];
			// gsum is the covariate for the gamma term
			// double gsum = 0;
			// for(const auto& l : b2b_mapping[j])
			// 	gsum += Z(i, l).mean();
			double MRFtrue = exp(beta1_ + gamma_ * gsum_mat(i,j));
			// mean of the likelihood of the replicates
			double prod = 1;
			for(unsigned rep = 0; rep < m; rep++)
				prod *= (k[rep] ? MRFtrue : MRFfalse)/ (MRFtrue + MRFfalse);
			loglik += log(prod)/m;
		}
	return loglik;
}

double Model_data::llikelihood() const {
	double loglik = 0;
	for(size_t i=0; i < nprey; ++i) {
		for(size_t j=0; j < nbait; ++j) {
			const auto& y = test_mat1(i,j);
			if(all_of(begin(y), end(y), saint::isnan)) continue;
			const auto& k = Z(i,j);
			const auto m = n_rep_vec[j];
			// gsum is the covariate for the gamma term
			double gsum = 0;
			// BOOST_FOREACH(const auto l , b2b_mapping[j])
			// 	gsum += Z(i, l).mean();
			BOOST_FOREACH(const auto l , p2p_mapping[i])
				gsum += Z(l, j).mean();
			double MRFtrue = exp(beta1 + gamma * gsum);
			double MRFfalse = exp(beta0);
			// mean of the likelihood of the replicates
			double prod = 1;
			for(unsigned rep = 0; rep < m; rep++){
				double y_rep = saint::isnan(y[rep]) ? t : static_cast<double>(y[rep]);//TODO
				// Quant_t y_rep = saint::isnan(y[rep]) ? t : y[rep];
				prod *= (k[rep] ? MRFtrue * exp(quant_log_pdf(min<double>(y_rep, (eta[i]+d[i])), (eta[i]+d[i]),sd_true[i])): MRFfalse * exp(quant_log_pdf(max<double>(y_rep, eta[i]), eta[i], sd_false[i])))/ (MRFtrue + MRFfalse);
			}
			loglik += log(prod)/m;
		}
	}
	return loglik;
}

double Model_data::loglikelihood_Z(const size_t i, const size_t j, const size_t rep, const Fastmat<vector<boost::array<double, 2> > >& pre_calc_loglik) const {
	const auto& k = Z(i,j);
	double gsum = 0;
	// BOOST_FOREACH(const auto l , b2b_mapping[j])
	// 	gsum += Z(i, l).mean();
	if(gamma!=0){
		BOOST_FOREACH(const auto l , p2p_mapping[i])
			gsum += Z(l, j).mean();
	}
	double logMRFtrue = (beta1 + gamma * gsum);
	double logMRFfalse = (beta0);
	const boost::array<double, 2>& pcl = pre_calc_loglik(i,j)[rep];
	return k[rep]? logMRFtrue + pcl[0] : logMRFfalse + pcl[1];
}

void Model_data::icm_Z() {
	// pre calculating log likelihood
	Fastmat<vector<boost::array<double, 2> > > pre_calc_loglik(nprey, nbait);
	for(size_t i=0; i < nprey; ++i)
		for(size_t j=0; j < nbait; ++j) {
			const auto& y = test_mat1(i,j);
			const auto m = n_rep_vec[j];
			pre_calc_loglik(i,j) = vector<boost::array<double, 2> >(m);
			for(unsigned rep = 0; rep < m; rep++) {
				double y_rep = saint::isnan(y[rep]) ? t : static_cast<double>(y[rep]);//TODO
				// Quant_t y_rep = saint::isnan(y[rep]) ? t : y[rep];
				pre_calc_loglik(i,j)[rep] = boost::array<double, 2>{{
						(quant_log_pdf(min<double>(y_rep, (eta[i]+d[i])), (eta[i]+d[i]),sd_true[i])),
						(quant_log_pdf(max<double>(y_rep, eta[i]), eta[i], sd_false[i]))}};
			}
		}
	for(unsigned iter_Z=0; iter_Z < 1; ++iter_Z)
		for(size_t i=0; i < nprey; ++i)
			for(size_t j=0; j < nbait; ++j) {
				const auto m = n_rep_vec[j];
				for(unsigned rep = 0; rep < m; rep++){
					double first = loglikelihood_Z(i,j,rep,pre_calc_loglik);
					Z(i,j)[rep] = !Z(i,j)[rep];
					double second = loglikelihood_Z(i,j,rep,pre_calc_loglik);
					if(first > second) Z(i,j)[rep] = !Z(i,j)[rep];
				}
			}
}

Fastmat<vector<boost::array<double, 2> > > Model_data::precalculate_densities() const {
	// precalculate densities
	Fastmat<vector<boost::array<double, 2> > > densities(nprey, nbait);
	for(size_t i=0; i < nprey; ++i)
		for(size_t j=0; j < nbait; ++j) {
			const auto& y = test_mat1(i,j);
			const auto m = n_rep_vec[j];
			densities(i,j).resize(m);
			for(unsigned rep = 0; rep < m; rep++){
				double y_rep = saint::isnan(y[rep]) ? t : static_cast<double>(y[rep]);//TODO
				// Quant_t y_rep = saint::isnan(y[rep]) ? t : y[rep];
				densities(i,j)[rep] = boost::array<double, 2>{{
						exp(quant_log_pdf(min<double>(y_rep, (eta[i]+d[i])), (eta[i]+d[i]),sd_true[i])),
						exp(quant_log_pdf(max<double>(y_rep, eta[i]), eta[i], sd_false[i]))}};
			}
		}
	return densities;
}

void Model_data::wrt_MRF_gamma_0() {
	nlopt::opt opt(nlopt::LN_COBYLA, 1);
	opt.set_lower_bounds({-15});
	opt.set_upper_bounds({ 15});

	opt.set_ftol_abs(1e-4);
	opt.set_maxeval(1e2);
	std::vector<double> x{beta1};
	double oldbeta1 = beta1;
	double oldf = llik_MRF_gamma_0(x, dummy_vector);
	double maxf = nlopt_wrap(opt, x, &Model_data::llik_MRF_gamma_0, this);
	beta1 = x[0];
	if(maxf < oldf)
		beta1 = oldbeta1;
}

void Model_data::wrt_MRF() {
	// maximize wrt βs, γ
	// precalculate gsum
	Fastmat<double> gsum_mat(nprey, nbait);
	for(size_t i=0; i < nprey; ++i)
		for(size_t j=0; j < nbait; ++j)
			// gsum_mat is the covariate for the gamma term
			// 		BOOST_FOREACH(const auto l , b2b_mapping[j])
			// 			gsum_mat(i,j) += Z(i, l).mean();
			BOOST_FOREACH(const auto l , p2p_mapping[i])
				gsum_mat(i,j) += Z(l, j).mean();
	nlopt::opt opt(nlopt::LN_COBYLA, 2);
	opt.set_lower_bounds({-15, 0});
	opt.set_upper_bounds({ 15, 10});

	// opt.set_xtol_rel(1e-5);
	opt.set_ftol_abs(1e-4);
	// opt.set_xtol_abs({1e-5,1e-5});
	opt.set_maxeval(1e2);
	std::vector<double> x{beta1, gamma};
	double oldbeta1 = beta1, oldgamma = gamma;
	double oldf = llik_MRF(x, dummy_vector, gsum_mat);
	double maxf = nlopt_wrap(opt, x, &Model_data::llik_MRF, this, gsum_mat);
	beta1 = x[0], gamma = x[1];
	if(maxf <= llik_MRF({x[0],0},dummy_vector, gsum_mat))
		gamma = 0;
	if(maxf < oldf)
		beta1 = oldbeta1, gamma = oldgamma;
}

void Model_data::wrt_d() {
	// maximize wrt d
	for(size_t i=0; i < nprey; ++i) {
		Quant_t sum = 0;
		size_t sum_count = 0;
		for(size_t j=0; j < nbait; ++j) {
			const auto& k = Z(i,j);
			const auto& y = test_mat1(i,j);
			const auto m = n_rep_vec[j];
			// mean of the likelihood of the replicates
			for(unsigned rep = 0; rep < m; rep++)
				if(k[rep]){
					// Quant_t y_rep = (saint::isnan(y[rep]) ? t : y[rep]);
					Quant_t y_rep = (y[rep]);
					sum += y_rep;
					sum_count++;
				}
		}
		if(sum_count>0)	d[i] = static_cast<double>(sum)/sum_count - eta[i];
		else d[i] = log(4);
		// else d[i] = sd_false[i] * 2;
		if(d[i] < log(4)) d[i] = log(4);
		// d[i] = max<double>(d[i],sd_false[i] * 2);
	}
}

boost::array<Fastmat<double>, 3> Model_data::calculateScore() const {
	Fastmat<double> average_score(nprey, nbait);
	Fastmat<double> min_log_odds_score(nprey, nbait);
	Fastmat<double> maximum_score(nprey, nbait);
	for(size_t i=0; i < nprey; ++i)
		for(size_t j=0; j < nbait; ++j) {
			// const auto& k = Z(i,j);
			const auto& y = test_mat1(i,j);
			const auto m = n_rep_vec[j];
			double gsum = 0;
			// {change to p2p
			// BOOST_FOREACH(const auto l , b2b_mapping[j])
			// 	gsum += Z(i, l).mean();
			// double MRFtrue = exp(beta1 + apply_MRF[i] * gamma * gsum);
			// }
			BOOST_FOREACH(const auto l , p2p_mapping[i])
				gsum += Z(l, j).mean();
			double MRFtrue = exp(beta1 + gamma * gsum);
			double MRFfalse = exp(beta0);
			// mean scoring
			average_score(i, j) = 0;
			min_log_odds_score(i, j) = 0;
			vector<double> tmp_scores(m),tmp_odds_scores(m);
			for(unsigned rep = 0; rep < m; ++rep) {
				double y_rep = saint::isnan(y[rep]) ? t : static_cast<double>(y[rep]);//TODO
				// Quant_t y_rep = saint::isnan(y[rep]) ? t : y[rep];
				double unnorm_score_true = MRFtrue * exp(quant_log_pdf( min<double>(y_rep, (eta[i]+d[i])), eta[i] + d[i], sd_true[i]));
				double unnorm_score_false = MRFfalse * exp(quant_log_pdf(max<double>(y_rep, eta[i]), eta[i], sd_false[i]));
				tmp_scores.at(rep) = unnorm_score_true / (unnorm_score_true + unnorm_score_false);
				if(saint::isnan(y[rep]))
					tmp_scores.at(rep) = 0;
				// tmp_odds_scores.at(rep)=unnorm_score_true/unnorm_score_false;
				tmp_odds_scores.at(rep)=log(unnorm_score_true)-log(unnorm_score_false);
			}
			extern Options opts;
			const size_t max_rep = opts.R;
			// const size_t max_rep = 100;
			// std::partial_sort(tmp_scores.begin(), tmp_scores.begin() + max_rep, tmp_scores.end(), greater<double>());
			std::sort(tmp_scores.begin(), tmp_scores.end(), greater<double>());
			std::sort(tmp_odds_scores.begin(), tmp_odds_scores.end(), greater<double>());
			if(m > max_rep){
				average_score(i, j) = std::accumulate(tmp_scores.begin(), tmp_scores.begin() + max_rep, 0.0 ) / ((double) max_rep);
				// min_log_odds_score(i, j) = std::accumulate(tmp_odds_scores.begin(), tmp_odds_scores.begin() + max_rep, 0.0 ) / ((double) max_rep);
			}else{
				average_score(i, j) = std::accumulate(tmp_scores.begin(), tmp_scores.end(), 0.0 ) / ((double) m);
				// min_log_odds_score(i, j) = std::accumulate(tmp_odds_scores.begin(), tmp_odds_scores.end(), 0.0 ) / ((double) m);
			}
			min_log_odds_score(i, j) = tmp_odds_scores.back();
			maximum_score(i, j) = tmp_scores.at(0);
		}
	return boost::array<Fastmat<double>, 3>{{average_score, maximum_score,min_log_odds_score}};
}
