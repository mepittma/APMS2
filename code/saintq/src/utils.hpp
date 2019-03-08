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


#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/io/ios_state.hpp>


#include <limits>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <string>

namespace util{

enum class newline_type : char {lf,crlf,cr,unknown};


template<typename CharT, typename Traits>
inline newline_type detect_line_ending_impl(::std::basic_istream<CharT,Traits>& is) {
    char tmp;
    while(is){
        is.get(tmp);
        if(tmp == is.widen('\r')) {	// old Mac or Windows
            is.get(tmp);
            if(tmp == is.widen('\n'))	// Windows
                return newline_type::crlf;
            return newline_type::cr;	// old Macs
        }
        if(tmp == is.widen('\n'))	// Unix and modern Macs
            return newline_type::lf;
    }
    return newline_type::unknown;
}


template<typename CharT, typename Traits>
class getline_fn{
public:
    getline_fn(::std::basic_istream<CharT,Traits>& i,
               const newline_type n):
        is(i), newline(n)
    {
        if(newline==newline_type::unknown)
            throw ::std::runtime_error("invalid delimiter");
    }
    template<typename Allocator>
    ::std::basic_istream<CharT,Traits>&
    operator()(::std::basic_string<CharT,Traits,Allocator>& line) const
    {
        switch(newline){
        case newline_type::lf:
            return ::std::getline(is, line, is.widen('\n'));
        case newline_type::crlf:
            ::std::getline(is, line, is.widen('\r'));
            assert(is.get()==is.widen('\n'));
            return is;
            return ::std::getline(is, line, is.widen('\r')).ignore(1);
        case newline_type::cr:
            return ::std::getline(is, line, is.widen('\r'));
        default:
            throw;
        }
    }
private:
    ::std::basic_istream<CharT,Traits>& is;
    const newline_type newline;
};

template<typename CharT, typename Traits>
inline
newline_type
detect_line_ending(::std::basic_istream<CharT,Traits>& infile)
{
    const auto state=infile.rdstate();
    const auto p= infile.tellg();
    const auto ret = detect_line_ending_impl(infile);
    infile.seekg(p);
    infile.setstate(state);
    return ret;
}


template<typename CharT, typename Traits>
inline
getline_fn<CharT, Traits>
make_getline(::std::basic_istream<CharT,Traits>& infile)
{
    return getline_fn<CharT, Traits>
        (infile,detect_line_ending(infile));
}

template<typename CharT, typename Traits,typename Allocator>
::std::basic_ifstream<CharT,Traits>&&
create_istream(const ::std::basic_string<CharT,Traits,Allocator>& fn)
{
    ::std::ifstream is(fn);
    is.exceptions(::std::ios_base::failbit);
    is.exceptions(::std::ios_base::badbit);
    return is;
}


template <class CharT, class Traits>
void check_init_istream(::std::basic_istream<CharT, Traits>& is)
{
    is.exceptions(::std::ios_base::failbit);
    is.exceptions(::std::ios_base::badbit);
}

template <class CharT, class Traits>
void check_init_ostream(::std::basic_ostream<CharT, Traits>& os)
{
    os.exceptions(::std::ios_base::failbit | ::std::ios_base::badbit);
}

void check_init_global_streams()
{
    check_init_ostream(::std::cout);
    check_init_ostream(::std::cerr);
    check_init_ostream(::std::clog);

    check_init_istream(::std::cin);

    check_init_ostream(::std::wcout);
    check_init_ostream(::std::wcerr);
    check_init_ostream(::std::wclog);

    check_init_istream(::std::wcin);
}


uint_least32_t count_lines(std::istream& is)
///// count lines from current position
{
    uint_least32_t linecount=0;
    const auto here=is.tellg();
    {
        boost::io::ios_iostate_saver ioss{is};
        boost::io::ios_exception_saver ioss2{is};

        is.exceptions(std::ios_base::badbit & is.exceptions());// only allow badbit exception
        auto prevpos=here;
        auto currpos=here;
        while(is.ignore(std::numeric_limits<std::streamsize>::max(),'\n')){
            prevpos=currpos;
            currpos=is.tellg();
            // const auto currpos=is.tellg();
            // if(currpos>=0)prevpos=currpos;
            ++linecount;
        }
        is.clear();
        // if last line is empty, will decrement linecount to be consistent with std::getline
        is.seekg(0,std::ios_base::end);// goto the end
        if(prevpos==is.tellg())
            --linecount;
    }
    is.seekg(here);
    return linecount;
}




template<typename C>
inline typename C::value_type sum(const C& v)
{
    static_assert(!std::is_same<typename C::value_type, bool>::value ,"container elements type should not be bool.");
    return std::accumulate(v.begin(),v.end(), typename C::value_type(0));
}



template<typename T>
inline T square(const T& x){return x*x;}

template<typename It,typename T=typename std::iterator_traits<It>::value_type>
std::pair<T,T> mean_sd(const It beg, const It end)
{
    // find mean and sd after omitting nans
    using namespace boost::accumulators;
    using ::std::sqrt; using ::std::isnan;
    accumulator_set<T,features<tag::mean, tag::variance, tag::count>> acc;
    for(auto it=beg;it!=end;++it)
        if(!isnan(*it)){
            acc((*it));
        }
    if(count(acc)==0)
        throw std::runtime_error("count==0");
    if(count(acc)==1)
        throw std::runtime_error("count==1,variance");
    return std::make_pair(mean(acc),sqrt(variance(acc)));

}



}