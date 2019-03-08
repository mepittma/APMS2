/*
Copyright (C) 2016 Guo Ci Teo <ci@u.nus.edu> and Hyungwon Choi
<hyung_won_choi@nuhs.edu.sg> National University of Singapore.

This file is part of SAINTq.

SAINTq is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

SAINTq is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with SAINTq.  If not, see <http://www.gnu.org/licenses/>.
*/


#pragma once

#include <vector>
#include <algorithm>
#include <iterator>

namespace bfdr{
namespace detail{

// pair of score, index
struct score_idx{
    double score;
    std::size_t idx;
    score_idx(){}// prevent zero initialzation
    score_idx(double score_,std::size_t idx_):
        score(score_),idx(idx_)
    {}
    bool operator>(const score_idx& o) const
    {return this->score > o.score;}
    bool operator<(const score_idx& o) const
    {return this->score < o.score;}
};


template<typename Out_RA_Range>
void bfdr_it_impl(std::vector<detail::score_idx>& v,Out_RA_Range& out_range){
    const auto n=v.size();

    std::sort(v.begin(),v.end(),std::greater<detail::score_idx>{});

    double last = out_range[v[0].idx] = 0;
    double psum=0;
    for(std::size_t i=0;i<n-1;++i){
        psum+=v[i].score;
        last=v[i].score==v[i+1].score ?
            last:
            1-psum/static_cast<double>(i+1);
        out_range[v[i+1].idx]=last;
    }
}

template<typename Out_RA_Range>
void bfdr_it_impl1(std::vector<detail::score_idx>& v,Out_RA_Range& out_range){
    const auto n=v.size();

    std::sort(v.begin(),v.end(),std::greater<detail::score_idx>{});


    std::vector<std::size_t> v_rank(n);
    /// v_rank is the number of scores bigger than self
    v_rank[0]=0;
    for(std::size_t i=0;i<n-1;++i)
        v_rank[i+1]=
            v[i].score==v[i+1].score ?
            v_rank[i] : i+1;

    std::vector<double> v_bfdr(n);
    v_bfdr[0]=0;
    {
        double psum=0;
        for(std::size_t i=0;i<n-1;++i){
            psum+=v[i].score;
            v_bfdr[i+1]=1-psum/static_cast<double>(i+1);
        }
    }

    for(std::size_t i=0;i<n;++i)
        out_range[v[i].idx] = v_bfdr[v_rank[i]];
}

}


class prob_acc{
public:
    prob_acc(const std::size_t n=0):v(),count(0)
    {v.reserve(n);}

    void operator()(const double p)
    {
        this->v.emplace_back(p,this->count);
        ++this->count;
    }

    std::vector<double> bfdr()
    {
        std::vector<double> out_range(this->count);
        detail::bfdr_it_impl(this->v,out_range);
        return out_range;
    }

private:
    std::vector<detail::score_idx> v;
    std::size_t count;
};

template<typename InputRange,typename Out_RA_Range>
void bfdr_it(const InputRange& scores,Out_RA_Range& out_range){
    ///// calculate BFDR
    const auto n = out_range.size();
    if(n==0) return;


    std::vector<detail::score_idx> v(n);
    {
        std::size_t i=0;
        for(auto it=scores.begin();it!=scores.end();++it){
            v[i].score=*it;
            v[i].idx=i;
            ++i;
        }
    }
    detail::bfdr_it_impl(v,out_range);
}


/*
template<typename InputRange,typename Out_RA_Range>
void bfdr_it(const InputRange& scores,Out_RA_Range& out_range){
    ///// calculate BFDR
    const auto n = out_range.size();
    if(n==0) return;

    typedef boost::counting_iterator<std::size_t> ci;
    static_assert(std::is_same<std::iterator_traits<ci>::iterator_category,
                  std::random_access_iterator_tag>::value,"");
    std::vector<std::size_t> idxs(ci(0), ci(n));
    std::sort(idxs.begin(),idxs.end(),
              [&scores](std::size_t i1,std::size_t i2)
              {return scores[i1]>scores[i2];});

    double last = out_range[idxs[0]] = 0;
    double psum=0;
    for(std::size_t i=0;i<n-1;++i){
        psum+=scores[idxs[i]];
        last=scores[idxs[i]]==scores[idxs[i+1]] ?
            last:
            1-psum/static_cast<double>(i+1);
        out_range[idxs[i+1]]=last;
    }
}
*/


template<typename ForwardRange>
std::vector<double> bfdr_vec(const ForwardRange& scores){
    const auto n=std::distance(scores.begin(),scores.end());
    std::vector<double> ret(n);
    bfdr_it(scores,ret);
    return ret;
}

}