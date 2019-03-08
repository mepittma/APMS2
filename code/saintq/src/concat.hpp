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
#include <cassert>

namespace concat{

template<typename T>
void
concat_impl(std::vector<T>&){}

template<typename T, typename Vec,typename... Tvecs>
void
concat_impl(std::vector<T>& ret
       ,const Vec& v1,const Tvecs&... args)
{
    ret.insert(ret.end(),v1.begin(),v1.end());
    concat_impl(ret,args...);
}

std::size_t get_size(){return 0;}

template<typename Vec,typename... Tvecs>
std::size_t
get_size(const Vec& v1,const Tvecs&... args)
{return v1.size()+get_size(args...);}

template<typename T,typename Vec,typename... Tvecs>
std::vector<T>
concat_ET(const Vec& v1,const Tvecs&... args)
///specify element type
{
    std::vector<T> ret;
    ret.reserve(get_size(v1,args...));
    concat_impl(ret,v1,args...);
    assert(ret.size()==get_size(v1,args...));
    return ret;
}


template<typename Vec,typename... Tvecs>
std::vector<typename std::remove_cv<typename Vec::value_type>::type>
concat(const Vec& v1,const Tvecs&... args)
{return concat_ET<typename std::remove_cv<typename Vec::value_type>::type>(v1,args...);}

}