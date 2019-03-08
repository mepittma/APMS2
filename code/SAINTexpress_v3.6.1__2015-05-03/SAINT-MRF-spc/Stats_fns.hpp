const unsigned log_factorial_max = 10000;
const double log_factorial_table[log_factorial_max + 1] =
#include "log_factorial.hpp"
	;
//LOG_FACTORIAL_TABLE;

inline double log_factorial(const unsigned short n){
	return (n > log_factorial_max)? lgamma(n+1) : log_factorial_table[n];
}


inline double GP_log_pmf(const unsigned short x, const double& lambda1, const double& lambda2) {
	if(x == 0) return -lambda1;
	// if(lambda2==0) return x * log(lambda1) - lambda1 - log_factorial(x);
	double tmp = lambda1 + x * lambda2;
	return log(lambda1) + (x - 1) * log(tmp) - tmp - log_factorial(x);
}

// generalized poisson parameterized by mean and lambda2
inline double GP_log_pmf1(const Count_t k, const double mean, const double lambda2 /*=0*/ ) {
	const double lambda1 = mean * (1 - lambda2);
	return GP_log_pmf(k, lambda1, lambda2);
}


inline double pois_log_pmf(Count_t& k, const double& lambda) {
	if(k == 0) return -lambda;
	return k * log(lambda) - lambda - log_factorial(k);
}


inline double gamma_log_pmf(const double &x, const double & k /* shape */, const double & theta /* scale */) {
	return -lgamma(k) - k * log(theta) + (k - 1) * log(x) - x/theta;
}

inline double gamma_log_pmf(const vector<double> &x_vector, const double & k /* shape */, const double & theta /* scale */) {
	double loglik = 0;
	BOOST_FOREACH(const double& x , x_vector)
		loglik += (k - 1) * log(x) - x/theta;
	return loglik + x_vector.size() * ( -lgamma(k) - k * log(theta) );
}


inline double GP_inv_cdf(const double &u, const double &lambda1, const double &lambda2) {
	size_t x = 0; double p = 0;
	while(true) {
		p += exp(GP_log_pmf(x, lambda1, lambda2));
		if( u < p ) return x;
		++x;
	}
}


inline double GP_log_pmf__(const unsigned short &x, const double &lambda1, const double &lambda2) {
	static double lambda1_old = 0, log_lambda1;
	if(x == 0) return -lambda1;
	double tmp = lambda1 + x * lambda2;
	if( lambda1_old == lambda1 ) {
		return log_lambda1 + (x - 1) * log(tmp) - tmp - log_factorial(x);
	}
	lambda1_old = lambda1;
	log_lambda1 = log(lambda1);
	return log_lambda1 + (x - 1) * log(tmp) - tmp - log_factorial(x);
}


inline double NB_log_pmf(const unsigned short &k, const double &r, const double &p) {
	return lgamma(k+r) - lgamma(r) - log_factorial(k) + r * log(1-p) + k * log(p);
}

template<typename T>
double mean(const valarray<T>& v) {
	return static_cast<double>(v.sum())/v.size();
}
// the function returns sample variance if it can be calculated, and return zero otherwise
template<typename T>
double var1(const valarray<T>& v) {
	if(v.size()<=1)
		return 0.0;
	valarray<double> tmp(v.size());
	for(size_t i=0; i<v.size();i++)
		tmp[i] = v[i];
	tmp -= mean(v);
	return (tmp*tmp).sum()/(v.size()-1.0);
}
template<typename T>
double var1(const vector<T>& v) {
	valarray<double> tmp(v.size());
	for(size_t i=0; i<v.size();i++)
		tmp[i] = v[i];
	return var1(tmp);
}
