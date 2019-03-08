#ifndef FASTMAT_HPP_
# define FASTMAT_HPP_
#include <valarray>
#include "globals.hpp"
namespace Fast_matrix {
	class Span {
	public:
		Span(){}
	};
	static Span __;
}

using Fast_matrix::__;
using Fast_matrix::Span;

template<typename T>
class Fastmat {
public:
	const size_t n_rows;
	const size_t n_cols;
	valarray<T> mat;
	
	size_t nrow()const{return n_rows;}
	size_t ncol()const{return n_cols;}
private:
	Fastmat<T>();
	const T& index_operator(const int &x, const int &y) const {

		return mat[y + x * n_cols];
	}

public:
	Fastmat<T>& operator=(const Fastmat<T>& other){
		if(other.n_rows == n_rows && other.n_cols == n_cols)
			mat = other.mat;
		else throw runtime_error("assignment operator error: dimensions mismatch");
		return *this;
	}
	
	Fastmat(const int n_rows_, const int n_cols_) :
		n_rows(n_rows_),
		n_cols(n_cols_),
		mat( n_rows * n_cols )
	{

	}
	Fastmat(const int n_rows_, const int n_cols_, const T& value) :
		n_rows(n_rows_),
		n_cols(n_cols_),
		mat(value, n_rows * n_cols) {}

	T min() const {
		return mat.min();
	}
	T max() const {
		return mat.max();
	}
	size_t size() const {
		return mat.size();
	}
	const T& operator[](const int &x) const {
		return mat[x];
	}
	T& operator[](const int &x) {
		return mat[x];
	}

	slice_array<T> operator()(const int &x, const Span &) {
		return mat[std::slice(x * n_cols, n_cols, 1)];
	}
	valarray<T> operator()(const int &x, const Span &) const {
		return mat[std::slice(x * n_cols, n_cols, 1)];
	}
	slice_array<T> operator()(const Span &, const int &y) {
		return mat[std::slice(y, n_rows, n_cols)];
	}
	valarray<T> operator()(const Span &, const int &y) const {
		return mat[std::slice(y, n_rows, n_cols)];
	}
	T& operator()(const int &x, const int &y) {
		return mat[y + x * n_cols];
	}
	const T& operator()(const int &x, const int &y) const {
		return mat[y + x * n_cols];
	}

	ostream& print(ostream& out = std::cout) const {
		// out << std::endl;
		for(size_t r = 0; r < n_rows; r++) {
			for(size_t c = 0; c < n_cols - 1; c++)
				out << (*this)(r,c) << "\t";
			out << (*this)(r, n_cols - 1) << '\n';
		}
		// return out << std::flush;
		return out;
	}

};

typedef Fastmat<Count_t> Countmat;

template<typename T>
ostream& operator<<(ostream& out, Fastmat<T> const& mat) {
	return mat.print(out);
}

template<typename T>
struct vector1 : vector<T> {
private:
	T zero_element;
public:
	vector1(const T& zero) : zero_element(zero) { }
	T sum() const {
		// T tmp_sum(this->front().size());
		T tmp_sum = zero_element;
		BOOST_FOREACH(const auto& e , *this)
			for(unsigned i=0;i<tmp_sum.size();i++)
				tmp_sum[i] += e[i];
		return tmp_sum;
		// return std::accumulate(this->begin(), this->end(), T(this->front().size()));
	}
	T mean() const {
		auto tmp = sum();
		// BOOST_FOREACH(auto& e , tmp) e /= this->size();
		for(size_t idx = 0; idx < tmp.size(); idx++)
			tmp[idx] /= this->size();
		return tmp;
	}
};

template<typename T>
struct valarray1 : valarray<T> {
	// using valarray<T>::valarray;
	valarray1(const valarray<T>& other) : valarray<T>(other){}
	double mean() const {
		return this->sum()/static_cast<double>(this->size());
	}
	double var() const {
		double mu = mean();
		valarray<T> tmp(*this - mu);
		return (tmp*tmp).sum()/(this->size()-1.0);
	}
	void push_back(const T& a) {
		valarray1<T> tmp = *this;
		size_t s = this->size();
		this->resize(s+1);
		size_t i = 0;
		BOOST_FOREACH(const auto& e , tmp) (*this)[i++] = e;
		(*this)[s] = a;
	}
};
template<>
struct valarray1<bool> : valarray<bool> {
	// using valarray<bool>::valarray;
	double mean() const {
		unsigned short sum_ = 0;
		// BOOST_FOREACH(const auto e , *this) sum_ += e;
		for(size_t idx = 0; idx < this->size(); idx++)
			sum_ += (*this)[idx];
		return static_cast<double>(sum_)/this->size();
	}
	
};



#endif // FASTMAT_HPP_
