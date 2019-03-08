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


#include "concat.hpp"
#include "utils.hpp"
#include "bmarr.hpp"
#include "bfdr.hpp"

// http://dlib.net/optimization_ex.cpp.html
#include <dlib/optimization.h>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>

#include <boost/convert.hpp>
#include <boost/convert/stream.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
#include <boost/multi_array.hpp>
#include <boost/algorithm/string.hpp>

#include <boost/range/algorithm.hpp>
#include <boost/range.hpp>


#include <unordered_map>
#include <unordered_set>
#include <string>
#include <array>
#include <fstream>
#include <typeinfo>
#include <vector>
#include <iterator>
#include <iostream>
#include <cmath>
#include <cfenv>

namespace saint4{
using ::std::size_t;

enum class output_level : char {protein, fragment};
enum class input_level : char {protein=1, peptide=2, fragment=3};

const bool output_debug=false;
const char* const nanstr = "";

typedef std::string str_t;
// const double smallnum=sqrt(std::numeric_limits<double>::min());
const double small_sd=std::exp2(-10);



typedef double intensity_t;
typedef intensity_t quant_t;
typedef bmarrw::bmarr<quant_t, 2> QuantArr;
typedef bmarrw::bmarr<bool,2> BoolArr2D;
typedef bmarrw::bmarr<QuantArr, 2> PPQuantArr;
typedef bmarrw::bmarr<double,2> scores_table_t;
typedef dlib::matrix<double,0,1> column_vector;

struct state_pdf{double false_,true_;};


///// create the boost converter object

const boost::cnv::cstream& cnv=[]()
    ->decltype(cnv)
{
    static boost::cnv::cstream lcnv;
    lcnv(::std::boolalpha);
    return (lcnv);
}();

struct param_t{
    ::std::string input_filename;
    bool normalize_control;
    input_level input_lvl;
    str_t protein_colname,
        pep_colname
        ,frag_colname;
    std::uintmax_t compress_n_ctrl,
        compress_n_rep,
        min_n_pep,
        min_n_frag;
    double best_prop_pep,
        best_prop_frag;
};

double cnv_to_floating(const str_t& s)try
{return boost::convert<double>(s,cnv).value();}
catch(const boost::bad_optional_access& e)
{throw std::runtime_error{"cannot convert '"+s+"' to floating"};}

std::intmax_t cnv_to_integer(const str_t& s)try
{return boost::convert<std::intmax_t>(s,cnv).value();}
catch(const boost::bad_optional_access& e)
{throw std::runtime_error{"cannot convert '"+s+"' to integer"};}

bool cnv_to_bool(const str_t& s)try
{return boost::convert<bool>(s,cnv).value();}
catch(const boost::bad_optional_access& e)
{throw std::runtime_error{"cannot convert '"+s+"' to bool"};}

template<typename Tmap>
typename Tmap::mapped_type
get_and_erase(Tmap& m, const typename Tmap::key_type& key)
// get a key, return the mapped value and remove entry from the map
{
    const auto it=m.find(key);
    if(it==m.end())
        throw std::out_of_range{"key not found:	"+key};
    const auto val=it->second;
    m.erase(it);
    return val;
}

input_level str_to_input_level(const str_t& s){
    const std::map<str_t,input_level> to_input_level{
        {"protein",input_level::protein},
        {"peptide",input_level::peptide},
        {"fragment",input_level::fragment}
    };
    try{return to_input_level.at(s);}
    catch(const std::out_of_range& e)
    {throw std::out_of_range{"invalid input_level:	"+s};}
}


param_t read_param(std::istream& is,std::ostream& oserr)
{
    ///// read in file as list of lines
    std::vector<std::string> lines;
    for(std::string line;getline(is,line);)
        lines.push_back(line);

    ///// read as key-value pair, check for duplicate keys
    std::map<str_t,str_t> kv;
    namespace balg=boost::algorithm;
    BOOST_FOREACH(const auto& line, lines){
        const auto line2=balg::trim_copy(line);
        if(balg::starts_with(line2,"#") || line2=="")
            continue;
        const auto result = balg::find_first(line,"=");
        const auto key=balg::trim_copy(str_t(line.begin(),result.begin()));
        const bool inserted=kv.insert(
            std::make_pair(
                key,
                balg::trim_copy(str_t(result.end(),line.end())))
            ).second;
        if(!inserted)
            throw std::runtime_error("duplicate key: "+key);
    }

    ///// parse strings
    const auto input_lvl=str_to_input_level(get_and_erase(kv,"input_level"));
    const auto protein_colname=get_and_erase(kv,"protein_colname");
    const auto pep_colname=
        input_lvl==input_level::peptide||input_lvl==input_level::fragment?
        get_and_erase(kv,"pep_colname"):"";
    const auto frag_colname=
        input_lvl==input_level::fragment?
        get_and_erase(kv,"frag_colname"):str_t{};
    const auto compress_n_ctrl=cnv_to_integer(get_and_erase(kv,"compress_n_ctrl"));
    const auto compress_n_rep=cnv_to_integer(get_and_erase(kv,"compress_n_rep"));
    const auto min_n_pep=
        input_lvl==input_level::peptide||input_lvl==input_level::fragment?
        cnv_to_integer(get_and_erase(kv,"min_n_pep")):
        std::numeric_limits<std::intmax_t>::max();
    const auto best_prop_pep=
        input_lvl==input_level::peptide||input_lvl==input_level::fragment?
        cnv_to_floating(get_and_erase(kv,"best_prop_pep")):
        ::std::numeric_limits<double>::quiet_NaN();
    const auto min_n_frag=
        input_lvl==input_level::fragment?
        cnv_to_integer(get_and_erase(kv,"min_n_frag")):
        std::numeric_limits<std::intmax_t>::max();
    const auto best_prop_frag=
        input_lvl==input_level::fragment?
        cnv_to_floating(get_and_erase(kv,"best_prop_frag")):
        ::std::numeric_limits<double>::quiet_NaN();
    const auto input_filename=get_and_erase(kv,"input_filename");
    const auto normalize_control=cnv_to_bool(get_and_erase(kv,"normalize_control"));
    ///// check for unused values
    if(!kv.empty()){
        oserr<<"unused parameters:\n";
        BOOST_FOREACH(const auto& e, kv)
            oserr<<"\t"<<e.first<<"\t=\t"<<e.second<<"\n";
        throw std::runtime_error{"unused parameters"};
    }

    ///// check for invalid inputs
    if(compress_n_ctrl<2)
        throw std::runtime_error{"compress_n_ctrl<2"};
    if(compress_n_rep<1)
        throw std::runtime_error{"compress_n_rep must be positive"};

    switch(input_lvl){
    case input_level::fragment:{
        if(min_n_frag<1)
            throw std::runtime_error("min_n_frag must be positive");
        const bool is_prop=0<best_prop_frag && best_prop_frag<=1;
        if(!is_prop)
            throw std::runtime_error("check that 0<best_prop_frag<=1");
    }
        /// fallthrough
    case input_level::peptide:{
        if(min_n_pep<1)
            throw std::runtime_error("min_n_pep must be positive");
        const bool is_prop=0<best_prop_pep && best_prop_pep<=1;
        if(!is_prop)
            throw std::runtime_error("check that 0<best_prop_pep<=1");
    }
    default:;
    }
    return param_t{input_filename,
            normalize_control,
            input_lvl,
            protein_colname,
            pep_colname,
            frag_colname,
            static_cast<std::uintmax_t>(compress_n_ctrl),
            static_cast<std::uintmax_t>(compress_n_rep),
            static_cast<std::uintmax_t>(min_n_pep),
            static_cast<std::uintmax_t>(min_n_frag),
            best_prop_pep,
            best_prop_frag
            };
}


template<typename TArr>
std::ostream& print_arr(std::ostream& os, const TArr& a){
    if(a.empty())
        return os;
    os<<a[0];
    for(size_t j=1;j<a.size();++j) os<<"\n"<<a[j];
    return os;
}

template<typename TArr2D>
std::ostream& print_2d_arr(std::ostream& os, const TArr2D& a)
{
    for(size_t i=0;i<a.size();++i){
        os<<a[i][0];
        // for(size_t j=1;j<a.shape()[1];++j)
        for(size_t j=1;j<a[0].size();++j)
            os<<"\t"<<a[i][j];
        os<<"\n";
    }
    return os;
}


bool is_missing(const double x)
{return std::isnan(x);}
// {return x==0;}

//// don't use 0
const quant_t missing_v=::std::numeric_limits<intensity_t>::quiet_NaN();


template<typename Tstring>
inline intensity_t string_to_intensity(const Tstring& s)
{
    if(s==nanstr)
        return missing_v;
    const intensity_t ret =
        // boost::lexical_cast<intensity_t>(s);
        boost::convert<intensity_t>(s,cnv).value();
    if(ret==0)// catches zeros and return nan;
        return missing_v;
    if(!::std::isfinite(ret))// catches all nan and infs
        throw std::runtime_error("not finite");
        // throw boost::bad_lexical_cast();
    return ret;
}

template<typename Tstring>
inline quant_t string_to_quant(const Tstring& s)
try{
    static_assert(std::is_same<intensity_t,quant_t>::value,"");
    return string_to_intensity(s);
// }catch(boost::bad_optional_access const&){
// }catch(boost::bad_lexical_cast&){
}catch(const std::exception& e){
    std::cerr<<e.what()<<std::endl;
    std::cerr<<s<<std::endl;
    throw;
}

std::array<intmax_t,2> i__(intmax_t a,intmax_t b)
{return std::array<intmax_t,2>{{a,b}};}


template<typename Range>
bool all_missing(const Range& r)
{
    BOOST_FOREACH(const auto& e, r)
        if(!is_missing(e))
            return false;
    return true;
}

template<typename Arr2D>
bool all_missing_table(const Arr2D& r)
{
    BOOST_FOREACH(const auto& e1, r){
        BOOST_FOREACH(const auto& e, e1)
            if(!is_missing(e))
                return false;
    }
    return true;
}


#include "prey_to_prey_map.hpp"


struct read_s1{
    // all columns, the column name rows and anything above are excluded in each column
    std::vector<std::vector<str_t>> columns;
    // first column index of column containing intensities
    std::size_t startcol;
    // 3 × #IP table,
    // row 0, test,ctrl status,
    // row 1, bait name
    // row 2, IP name
    std::vector<std::vector<str_t>> ip_info;
    // row 2
    std::vector<str_t> colnames;
};

struct rows_peptide_level_less{
    const size_t prot_col_idx;
    template<typename T>
    bool operator()(const T&a, const T& b)const
    {return a[prot_col_idx]<b[prot_col_idx];}
};
struct rows_frag_level_less{
    const size_t prot_col_idx,pep_col_idx;
    template<typename T>
    bool operator()(const T&a, const T& b)const
    {
        if(a[prot_col_idx]<b[prot_col_idx])
            return true;
        if(a[prot_col_idx]==b[prot_col_idx])
            return a[pep_col_idx]<b[pep_col_idx];
        return false;
    }
};

template<typename Trows,typename Tvec>
void sort_table_rows(Trows& rows, const Tvec& colnames, const param_t& param)
{
    const auto prot_col_idx=boost::find(colnames,param.protein_colname)-colnames.begin();
    const auto pep_col_idx=boost::find(colnames,param.pep_colname)-colnames.begin();
    switch(param.input_lvl){
    case input_level::protein:
        return;
    case input_level::peptide:
        boost::sort(rows,
            rows_peptide_level_less{static_cast<size_t>(prot_col_idx)});
        return;
    case input_level::fragment:
        boost::sort(rows,
            rows_frag_level_less{static_cast<size_t>(prot_col_idx),
                            static_cast<size_t>(pep_col_idx)});
        return;
    default:throw;
    }
}

read_s1 read_table1(
    std::istream& is,
    const param_t& param,
    std::ostream& oserr)
///// check structure
// read_s1 read_table1(const char* const strptr,size_t strlen,std::ostream& oserr)
{
    using namespace std;
    namespace balg=boost::algorithm;
    const char* const delim="\t";
    vector<vector<str_t>> rows;// rows of cells
    {
        vector<str_t> ss;// temp storage for a row of cells
        str_t line;
        while(getline(is,line))
            rows.push_back(balg::split(ss, line, balg::is_any_of(delim)));
    }
    const auto ncol = rows.at(2).size();// number of columns in file
    ///// check for unequal columns
    for(size_t i=0; i<rows.size(); ++i)
        if(rows[i].size()!=ncol){
            std::ostringstream ostrstream;
            ostrstream<<"check row "<<i;
            oserr<<ostrstream.str()<<std::endl;
            throw std::runtime_error(ostrstream.str());
        }
    const auto& colnames=rows[2];

    ///// check colnames has colnames in the parameter file, must be done before sorting table
    switch(param.input_lvl){
    case input_level::fragment:
        if(boost::find(colnames,param.frag_colname)==colnames.end())
            throw std::runtime_error{"frag_colname:\t'"+param.frag_colname+"' not found"};
    case input_level::peptide:
        if(boost::find(colnames,param.pep_colname)==colnames.end())
            throw std::runtime_error{"pep_colname:\t'"+param.pep_colname+"' not found"};
    case input_level::protein:
        if(boost::find(colnames,param.protein_colname)==colnames.end())
            throw std::runtime_error{"protein_colname:\t'"+param.protein_colname+"' not found"};
        break;
    default:
        throw;
    }

    ///// sort table, colnames in the parameter file must be checked for existance first
    if(1){
        auto r=boost::make_iterator_range(rows.begin()+3,rows.end());
        sort_table_rows(r,colnames,param);
        if(output_debug){
            std::ofstream of{"DEBUG_sort.tsv"};
            print_2d_arr(of,rows);
        }
    }

    size_t startcol=0;// first IP column number
    BOOST_FOREACH(const auto& e, rows.at(0))
        if(!e.empty())break;
        else ++startcol;
    const auto nIP=ncol-startcol;
    vector<vector<str_t>> ip_info(3,vector<str_t>(nIP));
    for(size_t r=0;r<ip_info.size();++r)
        for(size_t c=0;c<nIP;++c){
            // const auto substr=rows[r][startcol+c];
            // ip_info[r][c]=str_t(substr.begin(),substr.end());
            ip_info[r][c]=rows[r][startcol+c];
        }
    vector<vector<str_t>> columns(ncol,vector<str_t>(rows.size()-ip_info.size()));
    for(size_t i=0;i<rows.size()-ip_info.size();++i)
        for(size_t j=0;j<ncol;++j)
            columns[j][i]=std::move(rows[i+ip_info.size()][j]);

    read_s1 ret;
    ret.columns=move(columns);
    ret.colnames=colnames;
    ret.startcol=startcol;
    ret.ip_info=ip_info;
    return ret;
}



// template <typename TVec>
template <template <typename T,typename=::std::allocator<T>> class TVec>
std::pair<
TVec<str_t>,
std::vector<std::vector<size_t>>>
calc_one_to_many_idxes(const TVec<str_t>& v)
/// input [a,a,a,b,b,c]
/// output ([a, b, c], [[0, 1, 2], [3, 4], [5]])
/// replicates must appear consecutively
    // const std::vector<std::string>  vvv{"a","a","a","b","b","c"};
{
    TVec<str_t> unique_list;
    boost::unique_copy(v,back_inserter(unique_list));
    if(std::set<str_t>(unique_list.cbegin(),unique_list.cend()).size()!=unique_list.size())
        throw std::runtime_error{"replicates must be in consecutive columns"};
    std::vector<std::vector<size_t>> unique_idx_to_idxes(unique_list.size());
    std::map<str_t,std::vector<size_t>> names_to_idxes;
    for(size_t i=0;i<v.size();++i)
        names_to_idxes[v[i]].push_back(i);
    for(size_t i=0; i!=unique_list.size();++i)
        unique_idx_to_idxes[i]=names_to_idxes.at(unique_list[i]);
    return std::make_pair(unique_list,unique_idx_to_idxes);
}

struct read_s2{
    std::vector<str_t>
    // list of prey proteins, no duplicate, same order as input file
    prey_prot_list,
    // list of test bait proteins, no duplicate, same order as input file
        bait_list,
    // list of fragments, same order as input file
        prey_prot_col,
        pep_column,
        frag_column,
    // list of test bait names, for the table of test quantification table, list size == #test IPs
        test_table_bait_names;
    QuantArr test_quant, ctrl_quant;
    std::vector<std::vector<size_t>> bait_idx_to_colset_idxes, prot_idx_to_rowset_idxes;
};

read_s2 read_table2(const read_s1 s1,std::ostream& oserr,
                    const param_t& param)
///// check data type
{
    const auto input_lvl=param.input_lvl;
    const auto& protein_colname=param.protein_colname;
    const auto& pep_colname=param.pep_colname;
    const auto& frag_colname=param.frag_colname;

    ///// get quant, for ctrls and tests
    const auto& columns=s1.columns;
    QuantArr quant;
    quant.resize(i__(columns.at(0).size(),columns.size()-s1.startcol));
    for(size_t i=0;i<quant.shape()[0];++i)
        for(size_t j=0;j<quant.shape()[1];++j)
            quant(i__(i,j))= string_to_quant(columns[j+s1.startcol][i]);

    ///// classify IPs
    const auto& bait_status=s1.ip_info[0];
    std::vector<unsigned> test_col_index, ctrl_col_index;
    for(unsigned i=0;i<bait_status.size();++i){
        const auto& s=bait_status[i];
        if(s=="T")
            test_col_index.push_back(i);
        else if(s=="C")
            ctrl_col_index.push_back(i);
        else{
            oserr<<"bait status row must be 'T' or 'C', but read in"<<s<<std::endl;
            throw std::runtime_error("");
        }
    }
    if(0){
        const auto m_s=util::mean_sd(quant.data(),quant.data()+quant.num_elements());
        const auto mean=m_s.first;
        const auto sd=m_s.second;
        exit(0);
        ///// normalize
        for(auto it=quant.data();it!=quant.data()+quant.num_elements();++it)
            *it=(*it-mean)/sd;
    }

    ///// separate quant into test and ctrls
    QuantArr test_quant{i__(quant.shape()[0],test_col_index.size())};
    for(size_t i=0;i<quant.shape()[0];++i)
        for(size_t j=0; j<test_col_index.size(); ++j)
            test_quant(i__(i,j))=quant(i__(i,test_col_index[j]));

    for(size_t i=0;i<test_quant.shape()[0];++i)
        for(size_t j=0; j<test_quant.shape()[1]; ++j)
            test_quant(i__(i,j));
    // for(size_t i=0;i<test_quant.shape()[0];++i)
    // 	test_quant[boost::indices[i][test_quant.shape()[1]]];

    ///// bait names for test quant table, with same ordering.
    const auto& bait_row=s1.ip_info[1];
    std::vector<str_t> test_table_bait_names(test_col_index.size());
    for(unsigned i=0;i<test_col_index.size();++i)
        test_table_bait_names[i]=bait_row[test_col_index[i]];

    auto bbb=calc_one_to_many_idxes(test_table_bait_names);
    auto& bait_list=bbb.first;
    auto& bait_idx_to_colset_idxes=bbb.second;

    /// ctrl quant
    QuantArr ctrl_quant{i__(quant.shape()[0],ctrl_col_index.size())};
    for(size_t i=0;i<quant.shape()[0];++i)
        for(size_t j=0; j<ctrl_col_index.size(); ++j)
            ctrl_quant(i__(i,j))=quant(i__(i,ctrl_col_index[j]));

    // zero_self_interactions


    ///// get prey list
    std::map<str_t,unsigned> colname_to_col_idx;
    for(unsigned i=0;i<s1.startcol;++i){
        const auto it=colname_to_col_idx.find(s1.colnames[i]);
        if(it!=colname_to_col_idx.end()){
            oserr<<"duplicated column names: "<<it->first<<std::endl;
            throw std::runtime_error("");
        }
        colname_to_col_idx[s1.colnames[i]]=i;
    }

    const auto& prey_prot_col=columns[colname_to_col_idx.at(protein_colname)];
    const auto& pep_column=
        input_lvl==input_level::peptide || input_lvl==input_level::fragment
        ?
        columns[colname_to_col_idx.at(pep_colname)]:
        std::vector<str_t>{};
    const auto& frag_column=
        input_lvl==input_level::fragment?
        columns[colname_to_col_idx.at(frag_colname)]:
        std::vector<str_t>{};


    auto ppp=calc_one_to_many_idxes(prey_prot_col);
    auto& prey_prot_list=ppp.first;
    auto& prot_idx_to_rowset_idxes=ppp.second;

    return {prey_prot_list,
            bait_list,
            prey_prot_col, pep_column, frag_column,
            test_table_bait_names,
            test_quant, ctrl_quant,
            bait_idx_to_colset_idxes,
            prot_idx_to_rowset_idxes};
    // read_s2 ret;
    // return ret;
}



template<typename T>
inline bool isnan_template(const T d)
{return ::std::isnan(d);}

template<typename Tvec,typename T=typename Tvec::value_type>
T TPP_mtd(const Tvec& v){
    // TPP method for rolling up probabilities
    std::vector<T> v2(v.begin(),v.end());
    std::replace_if(v2.begin(),v2.end(),
                    // static_cast<bool(*)(T)>(&std::isnan),
                    &isnan_template<T>,
                    T(0.0));
    using ::std::sqrt;
    static const auto smallnum=sqrt(std::numeric_limits<T>::min());
    T tmp =1;
    BOOST_FOREACH(const auto& e, v2){
        if(tmp<smallnum)
            tmp=0;//
        tmp*=1-e;
    }
    return 1-tmp;
}

// template<typename T=double,typename Tvec>
template<typename Tvec,typename T=typename Tvec::value_type>
T mean_of_top_n(const Tvec& v,const std::uintmax_t n){
    // mean of top n elements of v, with nan replaced with zero
    std::vector<T> v2(v.begin(),v.end());
    const auto n2=std::min<std::uintmax_t>(n,v2.size());
    std::replace_if(v2.begin(),v2.end(),
                    // static_cast<bool(*)(T)>(&std::isnan),
                    &isnan_template<T>,
                    static_cast<T>(0.0));
    std::sort(v2.begin(),v2.end(),std::greater<T>{});
    return std::accumulate(
        v2.cbegin(),
        v2.cbegin()+n2,
        static_cast<T>(0.0))/static_cast<T>(n2);
}


template<typename T2DArr>
void impute_ctrl_prot_or_pep(T2DArr& t,const quant_t minval){
    using namespace boost::accumulators;
    typedef accumulator_set<typename T2DArr::value_type::value_type
            ,features<tag::min, tag::count>> acc_type;
    //// get prey protein and rep min
    std::vector<quant_t> rep_min(t[0].size());
    acc_type acc;
    for(unsigned i=0;i<rep_min.size();++i){
        acc_type acc_seq;
        for(unsigned j=0;j<t.size();++j){
            const auto x=t[j][i];
            if(!is_missing(x)){
                acc(x);
                acc_seq(x);
            }
        }
        rep_min[i]= count(acc_seq)==0 ? missing_v : min(acc_seq);
    }
    const auto prot_min=count(acc)==0 ? missing_v: min(acc);
    //// impute
    for(unsigned i=0;i<t.size();++i)
        for(unsigned j=0;j<t[i].size();++j){
            auto& x=t[i][j];
            if(is_missing(x)){
                if(!is_missing(rep_min[j]))
                    x=rep_min[j]*0.90;
                else if(!is_missing(prot_min))
                    x=prot_min*0.90;
                else
                    x=minval;
            }
        }
}

template<typename T2DArr>
void impute_single_prot0(T2DArr& t,const quant_t minval){
    ///// not used
    using namespace boost::accumulators;
    typedef accumulator_set<typename T2DArr::value_type::value_type
            ,features<tag::min, tag::count>> acc_type;
    ///// get protein and fragment min
    std::vector<quant_t> seq_min(t.size());
    acc_type acc;
    for(unsigned i=0;i<t.size();++i){
        acc_type acc_seq;
        for(unsigned j=0;j<t[i].size();++j){
            const auto x=t[i][j];
            if(!is_missing(x)){
                acc(x);
                acc_seq(x);
            }
        }
        seq_min[i]= count(acc_seq)==0 ? missing_v : min(acc_seq);
    }
//	if(count(acc)==0)
//	const auto prot_min=min(acc);
    const auto prot_min=count(acc)==0 ? missing_v: min(acc);
    ///// get fragment min, if possible
    for(unsigned i=0;i<t.size();++i)
        for(unsigned j=0;j<t[i].size();++j){
            auto& x=t[i][j];
            if(is_missing(x)){
                if(!is_missing(seq_min[i]))
                    x=seq_min[i];
                else if(!is_missing(prot_min))
                    x=prot_min/2;
                else
                    x=minval;
            }
        }
}

quant_t quant_min(const QuantArr& t)
///// minimum of a 2D array
{
    using namespace boost::accumulators;
    typedef accumulator_set<quant_t,features<tag::min, tag::count>> acc_type;
    acc_type acc;
    for(auto it=t.data();it!=t.data()+t.num_elements();++it)
        if(!is_missing(*it))
            acc(*it);
    if(count(acc)==0)
        throw std::runtime_error{"all missing in table"};
    return min(acc);
}

struct prey_frag_idx{
    /// prey prot index to row index
private:
    std::vector<size_t> iii;
public:
    prey_frag_idx(const std::vector<std::vector<size_t>>& x):
        iii(x.size())
    {
        for(size_t i=0; i<x.size();++i)
            iii[i]=x[i].front();
    }

    std::size_t operator()(std::size_t prey_idx,std::size_t frag_idx) const
    {return iii[prey_idx]+frag_idx;}
};

template<typename Range, typename T=typename Range::value_type>
std::pair<bool,T>
is_all_consecutive(const Range& r)
// test if all identical strings in a range are placed consecutively
{
    if(r.empty())
        throw std::runtime_error("empty range");
    auto prev=r.front();
    std::unordered_set<T> s;
    s.insert(r.front());
    BOOST_FOREACH(const auto& e, r){
        if(prev!=e){
            const bool inserted=s.insert(e).second;
            if(!inserted)
                return std::make_pair(false,e);
            prev=e;
        }
    }
    return std::make_pair(true,T{});
}

struct indexing_t{
private:
    typedef boost::multi_array_types::index_range range_t;
public:
    std::vector<size_t>
    bait_prot_row_idxs,
        prey_prot_row_idxs,
        pep_row_idxs,
        prot_pep_range;
    size_t num_bait_prots()const{return bait_prot_row_idxs.size()-1;}
    size_t num_prey_prots()const{return prey_prot_row_idxs.size()-1;}
    std::pair<size_t,size_t> prey_prot_row_range(size_t i) const
    {
        const auto it=prey_prot_row_idxs.begin()+i;
        return std::make_pair(it[0],it[1]);
    }
    range_t prey_prot_row_index_range(size_t i) const
    {
        const auto p=prey_prot_row_range(i);
        typedef boost::multi_array_types::index index;
        return index(p.first) <= range_t() < index(p.second);
        return range_t(p.first, p.second);
    }
    size_t prey_prot_num_peps(const size_t prot) const
    {
        const auto it=prot_pep_range.begin()+prot;
        return it[1]-it[0];
    }
    std::vector<size_t> prey_prot__frag_lengths(const size_t prot) const
    {
        const auto start_pep=prot_pep_range[prot];// first peptide index
        std::vector<size_t> v(prey_prot_num_peps(prot));
        for(size_t i=0;i<v.size();++i)
            v[i]=pep_row_idxs[start_pep+i+1]-pep_row_idxs[start_pep+i];
        // std::adjacent_difference();
        return v;
    }
};

indexing_t get_indexing(
    const read_s2& s2,
    const input_level input_lvl, std::ostream& oserr)
//// get indexing information
// start row/column number of all proteins + end index
{
    const auto& prey_prot_col=s2.prey_prot_col;
    const auto& pep_column=s2.pep_column;
    const auto& test_table_bait_names=s2.test_table_bait_names;

    const auto num_rows=prey_prot_col.size();

    ///// protein row boundaries
    std::vector<size_t> prey_prot_row_idxs;
    {
        prey_prot_row_idxs.push_back(0);
        std::unordered_set<str_t> s;
        for(size_t i=1;i<num_rows;++i){
            const auto& current=prey_prot_col[i];
            if(current!=prey_prot_col[i-1]){
                prey_prot_row_idxs.push_back(i);
                const bool inserted=s.insert(current).second;
                if(!inserted)
                    throw std::runtime_error("proteins not placed consecutively: "+current);
            }
        }
        prey_prot_row_idxs.push_back(num_rows);
    }
    //// for protein level data, each protein should only have 1 row
    if(input_lvl==input_level::protein)
        if(prey_prot_row_idxs.size()-1!=num_rows){
            throw std::runtime_error("protein level data should not have duplicate proteins");
        }

    ///// for fragment level input, check if all identical peptides within a protein are placed consecutively
    if(input_lvl==input_level::fragment){
        const auto beg=pep_column.cbegin();
        for(size_t i=0;i<prey_prot_row_idxs.size()-1;++i){
            const auto res=is_all_consecutive(
                boost::make_iterator_range(beg+prey_prot_row_idxs[i],beg+prey_prot_row_idxs[i+1]));
            if(!res.first){
                oserr<<res.second<<prey_prot_row_idxs[i]<<"\n"<<prey_prot_row_idxs[i+1]<<"\n";
                throw std::runtime_error{"identical peptides within a protein not placed consecutively"};
            }
        }
    }


    ///// peptide row boundaries
    std::vector<size_t> pep_row_idxs;
    if(input_lvl==input_level::peptide
       ||input_lvl==input_level::fragment){
        if(pep_column.size()!=num_rows)
            throw;
        pep_row_idxs.push_back(0);
        for(size_t i=1;i<num_rows;++i){
            const auto& current=pep_column[i];
            if(current!=pep_column[i-1] ||
               prey_prot_col[i]!=prey_prot_col[i-1])
                pep_row_idxs.push_back(i);
        }
        pep_row_idxs.push_back(num_rows);
    }
    ///// for peptide level data, in each protein, all peptides should only have one row
    if(input_lvl==input_level::peptide){
        if(pep_row_idxs.size()-1!=num_rows)
            throw std::runtime_error("peptide level data should not have more than one row for a protein");
    }

    ///// peptide number each protein owns, stored as a range of peptides
    std::vector<size_t> prot_pep_range;
    if(input_lvl==input_level::peptide
       ||input_lvl==input_level::fragment){
        prot_pep_range=std::vector<size_t>(prey_prot_row_idxs.size());
        for(size_t j=0,i=0; j<pep_row_idxs.size(); ++j)
            if(pep_row_idxs[j]==prey_prot_row_idxs[i])
                prot_pep_range[i++]=j;
    }

    ///// check if bait names are placed consecutively
    {
        const auto res=is_all_consecutive(test_table_bait_names);
        if(!res.first)
            throw std::runtime_error{"bait names must be places consecutively:\t"+res.second};
    }

    std::vector<size_t> bait_prot_row_idxs;
    bait_prot_row_idxs.push_back(0);
    for(size_t i=1;i<test_table_bait_names.size();++i)
        if(test_table_bait_names[i]!=test_table_bait_names[i-1])
            bait_prot_row_idxs.push_back(i);
    bait_prot_row_idxs.push_back(test_table_bait_names.size());


    return indexing_t{
        bait_prot_row_idxs, prey_prot_row_idxs,
            pep_row_idxs, prot_pep_range};
}

struct model_t{
    std::vector<double> mu_F, sd_F, mu_d, sd_T;
    /// does a protein-protein has any measurement
    BoolArr2D test_pp_score_pred;
    QuantArr test_quant2;
    // PPQuantArr test_pp_quant;
};

model_t calc(
    const read_s2 s2,
    const indexing_t& idx_info,
    const input_level input_lvl,
    std::ostream& /*oserr*/,
    const decltype(param_t::compress_n_ctrl) compress_n_ctrl/*number of ctrls to use, must be at least 2*/,
    const bool normalize_control
    )
{

    ///// numerator for normalization constant
    static_assert(std::numeric_limits<quant_t>::has_quiet_NaN,"");
    quant_t numer=std::numeric_limits<quant_t>::quiet_NaN();
    {
        namespace bacc=boost::accumulators;
        namespace tag=boost::accumulators::tag;
        typedef bacc::accumulator_set<quant_t
                                      ,bacc::features<tag::mean, tag::count>> acc_t;
        acc_t acc;
        const auto& test_quant=s2.test_quant;
        for(auto it=test_quant.data();it!=test_quant.data()+test_quant.num_elements();++it)
            if(!is_missing(*it))
                acc(*it);
        if(bacc::count(acc)>0)
            numer=bacc::mean(acc);
    }


    ///// denominator for normalization constant
    static_assert(std::numeric_limits<quant_t>::has_quiet_NaN,"");
    quant_t denom=std::numeric_limits<quant_t>::quiet_NaN();
    {
        namespace bacc=boost::accumulators;
        namespace tag=boost::accumulators::tag;
        typedef bacc::accumulator_set<quant_t
                                      ,bacc::features<tag::mean, tag::count>> acc_t;
        acc_t acc;
        const auto& ctrl_quant = s2.ctrl_quant;
        for(auto it=ctrl_quant.data();it!=ctrl_quant.data()+ctrl_quant.num_elements();++it)
            if(!is_missing(*it))
                acc(*it);
        if(bacc::count(acc)>0)
            denom=bacc::mean(acc);
    }

    const auto norm_const= numer/denom;
    auto ctrl_quant2 = s2.ctrl_quant;
    ///// normalization of controls
    if(normalize_control)
        for(auto it=ctrl_quant2.data();it!=ctrl_quant2.data()+ctrl_quant2.num_elements();++it)
            if(!is_missing(*it))
                *it *= norm_const;

    /////impute controls
    switch(input_lvl){
    case input_level::protein:
    case input_level::peptide:{
        /// a vector of prey protein of quants for controls
        // note: array_view cannot be default ctored.
        std::vector<QuantArr::array_view<2>::type> ctrls;
        ctrls.reserve(s2.prey_prot_list.size());
        const auto& prey_prot_row_idxs=idx_info.prey_prot_row_idxs;
        for(auto it=prey_prot_row_idxs.cbegin();it<prey_prot_row_idxs.cend()-1;++it){
            typedef boost::multi_array_types::index_range range_t;
            ctrls.emplace_back(ctrl_quant2[
                                   boost::indices
                                   [range_t(it[0],it[1])]
                                   [range_t(0,ctrl_quant2.shape()[1])]
                                   ]);
        }
        if(ctrls.size()!=idx_info.num_prey_prots())
            throw;
        const auto ctrl_min=quant_min(ctrl_quant2);
        BOOST_FOREACH(auto& e, ctrls)
            impute_ctrl_prot_or_pep(e, ctrl_min);
    }
        break;
    case input_level::fragment:{
        const auto& pep_row_idxs=idx_info.pep_row_idxs;
        std::vector<QuantArr::array_view<2>::type> ctrls;
        ctrls.reserve(pep_row_idxs.size()-1);
        for(auto it=pep_row_idxs.cbegin();it<pep_row_idxs.cend()-1;++it){
            typedef boost::multi_array_types::index_range range_t;
            ctrls.emplace_back(ctrl_quant2[
                                   boost::indices
                                   [range_t(it[0],it[1])]
                                   [range_t(0,ctrl_quant2.shape()[1])]
                                   ]);
        }
        if(ctrls.size()!=pep_row_idxs.size()-1)
            throw;
        const auto ctrl_min=quant_min(ctrl_quant2);
        BOOST_FOREACH(auto& e, ctrls)
            impute_ctrl_prot_or_pep(e, ctrl_min);
    }
        break;
    default:throw;
    }


    ///// log2 transform controls
    for(auto it=ctrl_quant2.data();it!=ctrl_quant2.data()+ctrl_quant2.num_elements();++it)
        if(!is_missing(*it))
            *it=log2(*it);

    ///// calculate means, sd of controls
    namespace bacc=boost::accumulators;
    namespace tag=boost::accumulators::tag;
    typedef bacc::accumulator_set<quant_t
            ,bacc::features<tag::mean, tag::variance, tag::count>> acc_t;
    std::vector<double> ctrl_means, ctrl_sds, ctrl_sds0;
    BOOST_FOREACH(const auto& r, ctrl_quant2){
        acc_t acc;
        auto r2=r;
        // boost::sort(r2,std::greater<quant_t>{});
        boost::partial_sort(r2,r2.begin()+std::min<size_t>(compress_n_ctrl,r2.size())
                            ,std::greater<quant_t>{});
        for(std::size_t i=0;i < std::min<size_t>(compress_n_ctrl,r2.size());++i)
            acc(r2[i]);
        ctrl_means.push_back(bacc::mean(acc));
        const auto sd=sqrt(bacc::variance(acc));
        ctrl_sds.push_back(sd);
        if(sd>small_sd)
            ctrl_sds0.push_back(sd);
    }


    ///// increase small sds for ctrls
    {
        ctrl_sds0.at(0);
        boost::sort(ctrl_sds0);
        const auto ctrl_sd_median=ctrl_sds0[ctrl_sds0.size()/2];//median
        BOOST_FOREACH(auto& sd, ctrl_sds)
            sd=std::max(sd,ctrl_sd_median);
    }

    //////////////////////////////////////////////////
    ///// test baits
    ///// impute test
    auto test_quant2=s2.test_quant;
    // if(input_lvl==input_level::fragment)
    // 	impute_test_frag(test_quant2);

    ///// log2 transform test quants
    for(auto it=test_quant2.data();it!=test_quant2.data()+test_quant2.num_elements();++it)
        if(!is_missing(*it))
            *it=log2(*it);


    ///// calculate means, sd of tests
    std::vector<double> test_means, test_sds0, test_sds;
    std::vector<size_t> test_counts;
    BOOST_FOREACH(const auto& r, test_quant2){
        acc_t acc;
        BOOST_FOREACH(const auto& e, r)
            if(!is_missing(e))
                acc(e);
        test_means.push_back(bacc::count(acc)>=1 ? bacc::mean(acc) : missing_v);
        const auto sd=bacc::count(acc)>=2 ? sqrt(bacc::variance(acc)) : missing_v;
        test_sds.push_back(sd);
        test_counts.push_back(bacc::count(acc));
    }

    BOOST_FOREACH(const auto& sd, test_sds)
        if(!is_missing(sd))
            test_sds0.push_back(sd);


    ///// increase small sd for tests
    test_sds0.at(0);
    {
        const auto nth_it=test_sds0.begin()+test_sds0.size()/2;
        boost::nth_element(test_sds0,nth_it);
        const auto test_sd_median=*nth_it;
        BOOST_FOREACH(auto& sd, test_sds)
            sd=std::max(is_missing(sd) ? 0 : sd, test_sd_median);
    }

    ///// get difference between means
    std::vector<double> mu_d;
    {
        //// get difference between means
        // auto& mu_d=ret.mu_d;
        mu_d.resize(s2.ctrl_quant.size());
        for(size_t i=0;i<mu_d.size();++i)
            mu_d[i]=is_missing(test_means[i])?
                std::log2(4):
                std::max(test_means[i]-ctrl_means[i],std::log2(4));
    }


    ///// bait protein - prey protein table
    std::vector<QuantArr::array_view<2>::type> aaa;
    aaa.reserve(s2.prey_prot_list.size() * s2.bait_list.size());// array_view has no default ctor, so do this
    BOOST_FOREACH(const auto& prot_idxes, s2.prot_idx_to_rowset_idxes){
        typedef boost::multi_array_types::index_range range_t;
        const range_t prey_prot_range(prot_idxes.front(),prot_idxes.back()+1);
        BOOST_FOREACH(const auto& bait_idxes, s2.bait_idx_to_colset_idxes){
            const auto bait_range=range_t(bait_idxes.front(),bait_idxes.back()+1);
            aaa.emplace_back(test_quant2[boost::indices[prey_prot_range][bait_range]]);
        }
    }


    if(0){
        BOOST_FOREACH(const auto& e, aaa){
            print_2d_arr(std::cout,e);
            std::cout<<"\n\n";
        }
    }
    // const_array_view has no copy assignment
    // bmarrw::bmarr<QuantArr::array_view<2>::type, 2> test_pp_quant_view(
    // 	i__(s2.prey_prot_list.size(),s2.bait_list.size()),aaa);
    boost::multi_array_ref<QuantArr::array_view<2>::type, 2>
        test_pp_quant_view{aaa.data(),i__(s2.prey_prot_list.size(),s2.bait_list.size())};

    ///// in test baits, check if a protein-protein has any measurement
    BoolArr2D test_pp_score_pred(i__(s2.prey_prot_list.size(),s2.bait_list.size()));
    for(size_t i=0;i<test_pp_score_pred.shape()[0];++i)
        for(size_t j=0; j<test_pp_score_pred.shape()[1];++j){
            // print_2d_arr(std::cout,test_pp_quant_view[i][j]);
            test_pp_score_pred[i][j] = !all_missing_table(test_pp_quant_view[i][j]);
        }
    ///// table of (prey protein × bait proteins) -> table of (replicates × prey protein group members)
    // test_pp_quant
    // {
    // 	///// table of (prey protein × bait proteins) -> table of (replicates × prey protein group members)
    // 	// auto& test_pp_quant=ret.test_pp_quant;
    // 	test_pp_quant.resize(i__(s2.prey_prot_list.size(),s2.bait_list.size()));
    // 	for(size_t i=0;i<test_pp_quant.shape()[0];++i)
    // 		for(size_t j=0; j<test_pp_quant.shape()[1];++j){
    // 			const auto& src=test_pp_quant_view[i][j];
    // 			auto& dest=test_pp_quant[i][j];
    // 			dest.resize(i__(src.shape()[1],src.shape()[0]));
    // 			///// copy, change indexing order to rep->(pep/frags)
    // 			for(size_t ii=0;ii<src.shape()[1];++ii)
    // 				for(size_t jj=0;jj<src.shape()[0];++jj)
    // 					dest[ii][jj]=src[jj][ii];
    // 		}
    // }


    return model_t{
        ctrl_means,ctrl_sds,mu_d,test_sds,
            test_pp_score_pred,
            test_quant2
            };
}




struct Model{
    Model(const std::vector<std::vector<size_t>>& x,
          BoolArr2D& test_pp_score_pred__,
          PPQuantArr& test_pp_quant__):
        row_idx_func(x),
        test_pp_score_pred(test_pp_score_pred__),
        test_pp_quant(test_pp_quant__)
        {}
    prey_frag_idx row_idx_func;
    std::vector<double> mu_F, sd_F, mu_d, sd_T;
    /// does a protein-protein has any measurement
    BoolArr2D& test_pp_score_pred;
    PPQuantArr& test_pp_quant;
};

Model calculations(
    /*not const modify test_quant*/ read_s2 s2,
    // const input_level input_lvl,
    std::ostream& /*oserr*/,
    BoolArr2D& test_pp_score_pred,
    PPQuantArr& test_pp_quant,
    const unsigned L/*number of ctrls to use*/)
{
    using std::vector;
    Model ret(
        ///// make indexing function
        s2.prot_idx_to_rowset_idxes,
        test_pp_score_pred,
        test_pp_quant);
    //// print quantities
    if(output_debug){
        {std::ofstream os{"DEBUG_ctrl_quant.tsv"};print_2d_arr(os,s2.ctrl_quant);}
        {std::ofstream os{"DEBUG_test_quant.tsv"};print_arr(os, s2.test_table_bait_names);print_2d_arr(os,s2.test_quant);}
    }

    ///// make indexing function
    // ret.row_idx_func=prey_frag_idx(s2.prot_idx_to_rowset_idxes);


    /////impute controls
    auto ctrl_quant2 = s2.ctrl_quant;
    {
        /// a vector of protein of quants for controls
        std::vector<QuantArr::array_view<2>::type> ctrls; ctrls.reserve(s2.prey_prot_list.size());
        BOOST_FOREACH(const auto& prot_idxes, s2.prot_idx_to_rowset_idxes){
            typedef boost::multi_array_types::index_range range;
            const range prey_prot_range(prot_idxes.front(),prot_idxes.back()+1);
            ctrls.emplace_back(ctrl_quant2[
                                   boost::indices
                                   [prey_prot_range]
                                   [range(0,ctrl_quant2.shape()[1])]
                                   ]);
        }
        const auto ctrl_min=quant_min(ctrl_quant2);
        BOOST_FOREACH(auto& e, ctrls)
            impute_ctrl_prot_or_pep(e, ctrl_min);
    }

    if(output_debug){
        std::ofstream os{"DEBUG_ctrl_quant_imputed.tsv"};print_2d_arr(os,ctrl_quant2);
    }

    //// log2 transform controls
    for(auto it=ctrl_quant2.data();it!=ctrl_quant2.data()+ctrl_quant2.num_elements();++it)
        if(!is_missing(*it))
            *it=log2(*it);

    //// calculate means, sd of controls
    namespace bacc=boost::accumulators;
    namespace tag=boost::accumulators::tag;
    typedef bacc::accumulator_set<quant_t
            ,bacc::features<tag::mean, tag::variance, tag::count>> acc_type;
    vector<double> ctrl_means, ctrl_sds, ctrl_sds0;
    // vector<size_t> ctrl_counts;
    BOOST_FOREACH(const auto& r, ctrl_quant2){
        acc_type acc;
        // BOOST_FOREACH(const auto& e, r)
        // 	acc(e);
        auto r2=r;
        boost::sort(r2,std::greater<quant_t>{});
        for(std::size_t i=0;i< std::min<size_t>(L,r.size());++i)
            acc(r2[i]);
        ctrl_means.push_back(bacc::mean(acc));
        const auto sd=sqrt(bacc::variance(acc));
        ctrl_sds.push_back(sd);
        if(sd>small_sd)
            ctrl_sds0.push_back(sd);
        // ctrl_counts.push_back(count(acc));
    }

    /// increase small sds for ctrls
    {
        ctrl_sds0.at(0);
        // const auto ctrl_sd_max=*boost::max_element(ctrl_sds);
        // BOOST_FOREACH(auto& e, ctrl_sds)
        // 	if(e<ctrl_sd_max/2) e=ctrl_sd_max;
        // for(size_t i=0;i<ctrl_sds.size();++i)
        // 	if(ctrl_sds[i]<ctrl_sd_max/2)
        // 		ctrl_sds[i]=ctrl_sd_max;
        boost::sort(ctrl_sds0);
        const auto ctrl_sd_median=ctrl_sds0[ctrl_sds0.size()/2];//median
        for(size_t i=0;i<ctrl_sds.size();++i)
            if(ctrl_sds[i]<ctrl_sd_median)
                ctrl_sds[i]=ctrl_sd_median;
    }
    if(output_debug){
        {std::ofstream os{"DEBUG_ctrl_means.tsv"}; print_arr(os,ctrl_means);}
        {std::ofstream os{"DEBUG_ctrl_sd.tsv"}; print_arr(os,ctrl_sds);}
    }

    ret.mu_F=ctrl_means;
    ret.sd_F=ctrl_sds;

    ///// log2 transform test quants
    for(auto it=s2.test_quant.data();it!=s2.test_quant.data()+s2.test_quant.num_elements();++it)
        if(!is_missing(*it))
            *it=log2(*it);

    //// calculate means, sd of tests
    vector<double> test_means, test_sds0;
    auto& test_sds=ret.sd_T;
    vector<size_t> test_counts;
    BOOST_FOREACH(const auto& r, s2.test_quant){
        acc_type acc;
        BOOST_FOREACH(const auto& e, r)
            if(!is_missing(e)) acc(e);
        test_means.push_back(bacc::count(acc)>=1 ? bacc::mean(acc) : missing_v);
        const auto sd=sqrt(bacc::variance(acc));
        test_sds.push_back(bacc::count(acc)>=2 ? sd : 0);
        if(bacc::count(acc)>=2)
            test_sds0.push_back(sd);
        test_counts.push_back(bacc::count(acc));
    }

    ///// increase small sd for tests
    test_sds0.at(0);
    boost::sort(test_sds0);
    {
        const auto test_sd_median=test_sds0[test_sds0.size()/2];
        BOOST_FOREACH(auto& e, test_sds)
            e=std::max(e,test_sd_median);
    }

    {
        //// get difference between means
        auto& mu_d=ret.mu_d;
        mu_d.resize(s2.ctrl_quant.size());
        for(size_t i=0;i<mu_d.size();++i)
            mu_d[i]=is_missing(test_means[i])?
                std::log2(4):
                std::max(test_means[i]-ctrl_means[i],std::log2(4));
    }

    ////// bait protein - prey protein table
    std::vector<QuantArr::array_view<2>::type> aaa;
    aaa.reserve(s2.prey_prot_list.size() * s2.bait_list.size());// array_view has no default ctor, so do this
    BOOST_FOREACH(const auto& prot_idxes, s2.prot_idx_to_rowset_idxes){
        typedef boost::multi_array_types::index_range range;
        const range prey_prot_range(prot_idxes.front(),prot_idxes.back()+1);
        BOOST_FOREACH(const auto& bait_idxes, s2.bait_idx_to_colset_idxes){
            const auto bait_range=range(bait_idxes.front(),bait_idxes.back()+1);
            aaa.emplace_back(s2.test_quant[boost::indices[prey_prot_range][bait_range]]);
        }
    }


    if(0){
        BOOST_FOREACH(const auto& e, aaa){
            print_2d_arr(std::cout,e);
            std::cout<<"\n\n";
        }
    }

    boost::multi_array_ref<QuantArr::array_view<2>::type, 2>
        test_pp_quant_view{aaa.data(),i__(s2.prey_prot_list.size(),s2.bait_list.size())};
    {
        ///// table of (prey protein × bait proteins) -> table of (replicates × prey protein group members)
        // auto& test_pp_quant=ret.test_pp_quant;
        test_pp_quant.resize(i__(s2.prey_prot_list.size(),s2.bait_list.size()));
        for(size_t i=0;i<test_pp_quant.shape()[0];++i)
            for(size_t j=0; j<test_pp_quant.shape()[1];++j){
                const auto& src=test_pp_quant_view[i][j];
                auto& dest=test_pp_quant[i][j];
                dest.resize(i__(src.shape()[1],src.shape()[0]));
                ///// copy, change indexing order to rep->(pep/frags)
                for(size_t ii=0;ii<src.shape()[1];++ii)
                    for(size_t jj=0;jj<src.shape()[0];++jj)
                        dest[ii][jj]=src[jj][ii];
            }
    }

    {
        ///// in test baits, check if a protein-protein has any measurement
        // auto& test_pp_score_pred=ret.test_pp_score_pred;
        test_pp_score_pred.resize(i__(s2.prey_prot_list.size(),s2.bait_list.size()));
        for(size_t i=0;i<test_pp_score_pred.shape()[0];++i)
            for(size_t j=0; j<test_pp_score_pred.shape()[1];++j)
                test_pp_score_pred[i][j]= !all_missing_table(test_pp_quant_view[i][j]);
    }

    //print_2d_arr(std::cerr,test_pp_score_pred);
    // Z.resize(i__(s2.prey_prot_list.size(),s2.bait_list.size()));

    return ret;
}



double quant_pdf(const double x,const double mu,const double sd)
{
    const double tmp=util::square((x-mu)/sd)/2;
    // return -tmp-std::log2(sd);//loglog

    const double tmp2=std::min<double>(tmp,32);
    return std::exp2(-tmp2)/sd;
}

double quant_pdf_true(double x,double mu,double sd)
{return quant_pdf(std::min(x,mu),mu,sd);}
double quant_pdf_false(double x,double mu,double sd)
{return quant_pdf(std::max(x,mu),mu,sd);}


enum class MRF_state : char{zero=0,one=1};

struct parameters{
    bmarrw::bmarr<MRF_state,2> Z;
    double /*beta0,*/ beta1, gamma;
    std::ostream& print_MRF_parameters(std::ostream& os/*,const parameters& p*/) const
    {
        return os
            //<<"beta0:"<<this->beta0
            <<"\tbeta1:"<<this->beta1
            <<"\tgamma:"<<this->gamma<<std::endl;
    }
};




typedef bmarrw::bmarr<state_pdf,2> likelihood_t;

likelihood_t get_quant_pdfs(const model_t& m){
    const auto dimarr=m.test_quant2.shape();
    likelihood_t ret{i__(dimarr[0],dimarr[1])};

    for(size_t i=0;i<dimarr[0];++i)
        for(size_t j=0;j<dimarr[1];++j){
            const auto e=m.test_quant2[i][j];
            const auto mu_F=m.mu_F[i];
            const auto mu_T=mu_F+m.mu_d[i];
            const auto sd_T=m.sd_T[i];
            const auto sd_F=m.sd_F[i];
            ret[i][j]=is_missing(e)?
                state_pdf{
                std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN()}:
            state_pdf{
                quant_pdf_false(e,mu_F,sd_F),
                    quant_pdf_true(e,mu_T,sd_T)};
        }
 	return ret;
}

double row_col_loglikelihood(
    // const bmarrw::bmarr<MRF_state,2>& Z,
    const likelihood_t& lik,
    const MRF_state z,
    const double MRF_true,const double MRF_false,
    const size_t i, const size_t j)
{
    const auto ll=lik(i__(i,j));
    return log(
        (z==::saint4::MRF_state::zero?
         MRF_false*ll.false_:MRF_true*ll.true_)
        /(MRF_true+MRF_false));
}

template<typename Tparameters>
double loglikelihood(const Tparameters& p,const model_t& m,const likelihood_t& lik)
{
    double loglik=0;
    const auto MRF_true=exp2(p.beta1/* + p.gamma * gsum*/);
    const auto MRF_false = exp2(0/*p.beta0*/);
    const auto dimarr=m.test_quant2.shape();
    assert(
        p.Z.shape()[0]==dimarr[0] &&
        p.Z.shape()[1]==dimarr[1]
        );
    for(size_t i=0;i<dimarr[0];++i)
        for(size_t j=0;j<dimarr[1];++j){
            if(is_missing(m.test_quant2(i__(i,j))))
                continue;
            const auto& z=p.Z(i__(i,j));
            // const auto gsum=static_cast<unsigned>(z);
            loglik += row_col_loglikelihood(lik,z,MRF_true,MRF_false,i,j);
        }
    return loglik;
}

double loglikelihood_Z(
    const likelihood_t& lik,
    const MRF_state z,
    const parameters& p,const model_t& /*m*/,
    const size_t i, const size_t j)
{
    const auto gsum=static_cast<unsigned>(z);
    // if(p.gamma!=0){}

    const auto MRF_true=exp2(p.beta1 + p.gamma * gsum);
    const auto MRF_false = exp2(0/*p.beta0*/);
    return row_col_loglikelihood(lik,z,MRF_true,MRF_false,i,j);
}

void icm_Z(parameters& p,const model_t& m,const likelihood_t& lik)
{
    const auto dimarr=m.test_quant2.shape();
    for(size_t i=0;i<dimarr[0];++i)
        for(size_t j=0;j<dimarr[1];++j){
            if(!m.test_quant2[i][j])
                continue;
            std::array<double,2> ll;
            for(unsigned st=0;st<ll.size();++st)
                ll[st]=loglikelihood_Z(lik,static_cast<MRF_state>(st), p, m, i, j);
            p.Z(i__(i,j))=static_cast<MRF_state>(boost::max_element(ll)-ll.begin());
        }
}


struct llik_MRF_gamma_0__{
    const parameters& p0;
    const model_t& m;
    const likelihood_t& lik;
    double operator()(const column_vector& beta1) const
    {
        param p=p0;
        p.beta1 = beta1(0);
        return loglikelihood(p,m,lik);
    }
    struct param{
        const decltype(parameters::Z)& Z;
        // const double beta0;
        double beta1;
        const double gamma;
        param(const parameters& p):
            Z(p.Z),/*beta0(p.beta0),*/beta1(p.beta1),gamma(0){}
    };
};


double wrt_MRF_gamma_0(parameters& p,const model_t& m,const likelihood_t& lik)
{
    column_vector starting_point(1);
    starting_point=p.beta1;
    const auto func=llik_MRF_gamma_0__{p,m,lik};
    const auto ret=dlib::find_max_box_constrained(
        dlib::bfgs_search_strategy{},
        dlib::objective_delta_stop_strategy(1e-7),
        func,
        dlib::derivative(func),
        starting_point,
        -2,2
        );
    p.beta1=starting_point(0);
    return ret;
}


struct results_t{
    scores_table_t avg_scores;
    BoolArr2D test_pp_score_pred;
    // QuantArr test_quant;
};

results_t scoring(
    const model_t m,
    const indexing_t& idx_info,
    const param_t& param,
    std::ostream& /*oserr*/)
{
    const auto input_lvl=param.input_lvl;
    const auto
        compress_n_rep=param.compress_n_rep,
        min_n_pep=param.min_n_pep,
        min_n_frag=param.min_n_frag;
    const auto lik=get_quant_pdfs(m);
    {
        std::vector<double> ttt;
        for(auto it=lik.data();it!=lik.data()+lik.num_elements();++it)
            if(!std::isnan(it->true_))
                ttt.push_back(it->true_);
        boost::sort(ttt);
        // {std::ofstream of("debug"); for(auto&e:ttt)of<<e<<"\n";}
    }

    const auto dimarr = m.test_quant2.shape();
    auto p=parameters();
    p.Z.resize(i__(dimarr[0],dimarr[1]));


    ///// ICM
    double oldllik, newllik = loglikelihood(p,m,lik);
    for (unsigned iter = 0;iter < 10000 ; ++iter) {
        oldllik = newllik;
        icm_Z(p,m,lik);
        newllik=wrt_MRF_gamma_0(p,m,lik);
        // newllik = loglikelihood(p,m,lik);

        if( (newllik >= oldllik) &&
            std::abs(newllik - oldllik) < 0.01 ) break;
    }

    ///// scoring
    const double MRF_true = exp(p.beta1/* + p.gamma * gsum*/);
    const double MRF_false = exp(0/*beta0*/);
    bmarrw::bmarr<double,2> scores{
        i__(idx_info.num_prey_prots(),idx_info.num_bait_prots())};

    bmarrw::bmarr<double,2> max_scores{i__(scores.shape()[0],scores.shape()[1])};

    for(size_t i=0;i<scores.shape()[0];++i)//global prey prot index
        for(size_t j=0;j<scores.shape()[1];++j){//global bait prot index
            // for ith prey and jth protein
            const auto prey_idx_range=idx_info.prey_prot_row_idxs.begin()+i;
            const auto prey_size=prey_idx_range[1]-prey_idx_range[0];
            const auto bait_idx_range=idx_info.bait_prot_row_idxs.begin()+j;
            const auto rep_size=bait_idx_range[1]-bait_idx_range[0];

            const auto prey_prot_row_offset=prey_idx_range[0];
            const auto col_offset=bait_idx_range[0];
            /// all scores associated with a prey prot to bait prot
            /// can be prot, pep, or frag level data
            bmarrw::bmarr<double,2> raw_score{i__(prey_size,rep_size)};
            for(size_t i2=0;i2<prey_size;++i2)//prey prot level row idx
                for(size_t j2=0;j2<rep_size;++j2){// bait prot level rep col idx
                    const auto ll=lik(i__(prey_prot_row_offset+i2,col_offset+j2));
                    const double unnorm_score_false = MRF_false * ll.false_;
                    const double unnorm_score_true = MRF_true * ll.true_;
                    raw_score(i__(i2,j2))= unnorm_score_true / (unnorm_score_true + unnorm_score_false);
                }
            switch(input_lvl){
                typedef boost::multi_array_types::index_range range_t;
            case input_level::fragment:{
                const auto ppr=idx_info.prot_pep_range.begin()+i;
                bmarrw::bmarr<double,2> pep_score{i__(ppr[1]-ppr[0],rep_size)};
                for(size_t i2=0;i2<pep_score.shape()[0];++i2)// prey prot level pep idx
                    for(size_t j2=0;j2<rep_size;++j2){// bait prot level rep col idx
                        /// global row index of this pep group
                        const auto pep_r=idx_info.pep_row_idxs.begin()+ppr[0]+i2;
                        const auto pep_grp_view=raw_score[
                            boost::indices
                            [range_t(pep_r[0],pep_r[1])-prey_prot_row_offset]
                            [j2]];
                        const auto best_n_frag2=std::lround(
                            std::max(
                                static_cast<double>(min_n_frag),
                                param.best_prop_frag
                                *
                                static_cast<double>(pep_grp_view.size())));
                    pep_score(i__(i2,j2))= mean_of_top_n(pep_grp_view,best_n_frag2);
                        // pep_score(i__(i2,j2))= TPP_mtd(pep_grp_view);
                    }
                raw_score=pep_score;
            }
                /// fallthrough
            case input_level::peptide:{
                const auto dims=raw_score.shape();
                const auto pep_size=dims[0];
                bmarrw::bmarr<double,2> preyprot_baitrep_score(i__(1,rep_size));
                for(size_t i2=0;i2<rep_size;++i2){
                    const auto column_view=raw_score[boost::indices[range_t(0,pep_size)][i2]];
                    const auto best_n_pep2=std::lround(
                        std::max(
                            static_cast<double>(min_n_pep),
                            param.best_prop_pep
                            *
                            static_cast<double>(column_view.size())));
                    preyprot_baitrep_score(i__(0,i2))=mean_of_top_n(column_view, best_n_pep2);
                    // preyprot_baitrep_score(i__(0,i2))=TPP_mtd(column_view);
                }
                raw_score=preyprot_baitrep_score;
            }
                /// fallthrough
            case input_level::protein:{
                assert(raw_score.shape()[0]==1);
                // scores(i__(i,j))=mean_of_top_n(raw_score[0],compress_n_rep);
                const auto score=mean_of_top_n(raw_score[0],compress_n_rep);
                scores(i__(i,j))=::std::round(score*10000.0)/10000.0;//rounding to 4th digit
            }
                break;
            default:
                throw;
            }
        }
    return results_t{scores, m.test_pp_score_pred};
}





template<typename Tcontainer,typename Tstring>
std::ostream& print_column_impl(std::ostream& os,const Tcontainer& v,const Tstring& delim, std::true_type)
{
    const auto nan_output_string=".";
    if(v.empty())
        return os;
    auto it=v.begin();
    if(is_missing(*it))
        os<<nan_output_string;
    else
        os<<*it;
    for(++it;it!=v.end();++it)
        if(is_missing(*it))
            os<<delim<<nan_output_string;
        else
            os<<delim<<*it;
    return os;
}
template<typename Tcontainer,typename Tstring>
std::ostream& print_column_impl(std::ostream& os,const Tcontainer& v,const Tstring& delim, std::false_type)
{
    if(v.empty())
        return os;
    auto it=v.begin();
    os<<*it;
    for(++it;it!=v.end();++it)
        os<<delim<<*it;
    return os;
}

template<typename Tcontainer,typename Tstring>
std::ostream& print_column(std::ostream& os,const Tcontainer& v,const Tstring& delim)
{return print_column_impl(os,v,delim, std::is_floating_point<typename Tcontainer::value_type>{});}


template<typename Tcontainer,typename Tstring>
struct print_col_type{
    const Tcontainer& v;
    const Tstring& delim;
    print_col_type(const Tcontainer& v_,const Tstring& delim_):
        v(v_),delim(delim_){}
    std::ostream& operator()(std::ostream& os) const
    {return print_column(os,v,delim);};
};

template<typename Tcontainer,typename Tstring>
print_col_type<Tcontainer,Tstring>
make_print_column(const Tcontainer& v,const Tstring& delim)
{return print_col_type<Tcontainer,Tstring>(v,delim);}

template<typename Tcontainer,typename Tstring>
std::ostream& operator<<(std::ostream& os,const print_col_type<Tcontainer,Tstring>& obj)
{return obj(os);}

template<typename T>
std::vector<intensity_t>
get_test_intensity(const T& t)
//// flatten a rectangular array into a column by summing
// t is a rep × pep table of a protein quants
{
    std::vector<intensity_t> ret(t.size());
    QuantArr t2(i__(t.size(),t[0].size()));
    for(size_t i=0;i<t.size();++i)
        for(size_t j=0;j<t[i].size();++j){
            const auto& e=t[i][j];
            /// sub missing with 0, undo log2 transform.
            t2[i][j]=is_missing(e)?0:exp2(e);
        }
    for(size_t i=0;i<t2.size();++i)
        ret[i]=util::sum(t2[i]);
    return ret;
}

template<typename T>
std::vector<intensity_t>
get_ctrl_intensity(const T& t)
//// flatten a rectangular array into a column by summing
// t is a rep × pep table of a protein quants
{
    std::vector<intensity_t> ret(t[0].size());
    QuantArr t2(i__(t[0].size(),t.size()));
    for(size_t i=0;i<t.size();++i)
        for(size_t j=0;j<t[i].size();++j){
            const auto& e=t[i][j];
            //transpose
            t2[j][i]=is_missing(e)?0:e;
        }
    for(size_t i=0;i<t2.size();++i)
        ret[i]=util::sum(t2[i]);
    return ret;
}

template<typename T>
QuantArr
exp2_with_nan(const T& t)
// for each element, if not missing, exp2, else do nothing
{
    QuantArr t2(i__(t.size(),t[0].size()));
    for(size_t i=0;i<t.size();++i)
        for(size_t j=0;j<t[i].size();++j){
            const auto& e=t[i][j];
            t2[i][j]=is_missing(e)?e:exp2(e);
        }
    return t2;
}

struct Exp2_with_nan{
    QuantArr t2;
    template<typename T>
    Exp2_with_nan(const T& t):
        t2(i__(t.size(),t[0].size()))
    {}

    template<typename T>
    QuantArr& operator()()
    // for each element, if not missing exp2, else do nothing
    {
        for(size_t i=0;i<t2.size();++i)
            for(size_t j=0;j<t2[i].size();++j){
                const auto& e=t2[i][j];
                t2[i][j]=is_missing(e)?e:exp2(e);
            }
        return t2;
    }
};

template<typename T>
QuantArr
transpose_sub0(const T& t)
// transpose and for each element, if missing set to 0
{
    QuantArr t2(i__(t[0].size(),t.size()));
    for(size_t i=0;i<t.size();++i)
        for(size_t j=0;j<t[i].size();++j){
            const auto& e=t[i][j];
            t2[j][i]=is_missing(e)?0:e;
        }
    return t2;
}

template<typename T>
QuantArr
transpose(const T& t)
{
    QuantArr t2(i__(t[0].size(),t.size()));
    for(size_t i=0;i<t.size();++i)
        for(size_t j=0;j<t[i].size();++j)
            t2[j][i]=t[i][j];
    return t2;
}

template<typename T>
QuantArr
missing_to_0(const T& t)
/// sub missing with 0
{
    QuantArr t2(i__(t.size(),t[0].size()));
    for(size_t i=0;i<t.size();++i)
        for(size_t j=0;j<t[i].size();++j){
            const auto& e=t[i][j];
            t2[i][j]=is_missing(e)?0:e;
        }
    return t2;
}

template<typename T>
std::vector<intensity_t>
row_sums(const T& t)
{
    std::vector<intensity_t> ret(t.size());
    for(size_t i=0;i<t.size();++i)
        ret[i]=util::sum(t[i]);
    return ret;
}


double odds(double p)
{
    p=std::min(p,0.9999);
    return p/(1-p);
}

scores_table_t get_bfdr(
    const BoolArr2D& test_pp_score_pred,
    const scores_table_t& scores)
{
    const auto dimarr=test_pp_score_pred.shape();

    //// get bfdr
    scores_table_t bfdr_2D_arr(i__(dimarr[0],dimarr[1]));

    std::vector<double> scores_without_nans;
    scores_without_nans.reserve(bfdr_2D_arr.num_elements());
    for(size_t i=0;i<dimarr[0];++i)
        for(size_t j=0;j<dimarr[1];++j)
            if(test_pp_score_pred(i__(i,j)))
                scores_without_nans.push_back(scores[i][j]);
    const auto bfdrs=bfdr::bfdr_vec(scores_without_nans);
    auto bfdr_it=bfdrs.begin();
    for(size_t i=0;i<dimarr[0];++i)
        for(size_t j=0;j<dimarr[1];++j)
            if(test_pp_score_pred(i__(i,j)))
                bfdr_2D_arr[i][j] = *bfdr_it++;

    return bfdr_2D_arr;
}



void output(
    std::ostream& os,
    const results_t& res,
    const read_s2& s2,
    const indexing_t& idx_info,
    const input_level input_lvl,
    const output_level lvl,
    const bool more)
{
    const auto colnames_p0_more=
        {"Bait", "Prey",/*"PreyFrag",*/ "Intensity", "IntensitySum", "AvgIntensity", "ctrlIntensity"};
    const auto colnames_p0_less = {"Bait", "Prey"};
    const auto& colnames_p0=more?colnames_p0_more:colnames_p0_less;
    const auto colnames_nums={"#Rep", "#Pep", "#Frag"};
    const auto colnames_scores={"AvgP", "BFDR"};

    static_assert(static_cast<int>(input_level::protein)==1,"");
    static_assert(static_cast<int>(input_level::peptide)==2,"");
    static_assert(static_cast<int>(input_level::fragment)==3,"");

    const auto column_names = concat::concat(
        std::vector<const char*>(colnames_p0.begin(),colnames_p0.end()),
        // more?colnames_p0_more:colnames_p0_less,
        boost::make_iterator_range(
            colnames_nums.begin(),
            colnames_nums.begin()+static_cast<int>(input_lvl)
            ),
        colnames_scores
        );

    const auto& test_pp_score_pred=res.test_pp_score_pred;
    const auto& test_quant=s2.test_quant;
    const auto& avg_scores=res.avg_scores;
    // const auto& max_scores=res.max_scores;
    const auto delim = "\t";
    const auto dimarr=test_pp_score_pred.shape();
    // const auto& pep_column=s2.pep_column;
    const auto& prey_prot_list=s2.prey_prot_list;
    const auto& bait_list=s2.bait_list;
    const auto& ctrl_quant=s2.ctrl_quant;

    const scores_table_t bfdr_2D_arr=get_bfdr(test_pp_score_pred,avg_scores);

    os<<std::setprecision(4)<<std::fixed;
    os<<make_print_column(column_names, delim)<<"\n";


    for(size_t i=0;i<dimarr[0];++i)
        for(size_t j=0;j<dimarr[1];++j){
            if(!test_pp_score_pred(i__(i,j)))
                continue;
            typedef boost::multi_array_types::index_range range_t;
            const auto prey_prot_rr_ir=idx_info.prey_prot_row_index_range(i);
            const auto bait_prot_r=idx_info.bait_prot_row_idxs.begin()+j;
            const auto q=transpose(test_quant[boost::indices
                                              [prey_prot_rr_ir]
                                              [range_t(bait_prot_r[0],bait_prot_r[1])]]);
            const auto test_intensities=/*sum frags*/row_sums(missing_to_0(q));
            const auto test_intensities_table_frag_rep=transpose(q);
            const auto sum_int=util::sum(test_intensities);
            typedef boost::multi_array_types::index_range range_t;
            const auto ctrl_q=ctrl_quant[boost::indices
                                         [prey_prot_rr_ir]
                                         [range_t(0,ctrl_quant.shape()[1])]];
            const auto ctrl_intensities_table=transpose_sub0(ctrl_q);
            const auto ctrl_intensities=/*sum frags*/row_sums(ctrl_intensities_table);
            /// protein level
            const auto avg_score=avg_scores[i][j];
            const auto bfdr__ = bfdr_2D_arr[i][j];
/*
            std::vector<str_t> prey_frags(q.shape()[1]);
            for(size_t frag=0; frag<q.shape()[1]; ++frag)
                prey_frags[frag]=pep_column[prey_prot_rr_ir.start()+frag];
*/
            // if(lvl==output_level::protein)
                os<<bait_list[j]
                  <<delim<<prey_prot_list[i];
            if(more){
                // os<<delim<<make_print_column(prey_frags,"|");
                os<<delim<<make_print_column(test_intensities,"|");
                os<<delim<<sum_int;
                os<<delim<<sum_int/static_cast<double>(test_intensities.size());
                os<<delim<<make_print_column(ctrl_intensities,"|");
            }
            os<<delim<<q.size();
            if(input_lvl==input_level::peptide||input_lvl==input_level::fragment)
                os<<delim<<idx_info.prey_prot_num_peps(i);
            if(input_lvl==input_level::fragment)
                os<<delim<<make_print_column(idx_info.prey_prot__frag_lengths(i),"|");
            os<<delim<<avg_score
                  // <<delim<<max_scores[i][j]
                  // <<delim<<odds(avg_score)
              <<delim<<bfdr__
                  <<"\n";
            // os<<i<<","<<j<<""<<"\n";
            /// peptide/fragment level
            if(lvl==output_level::fragment)
                for(size_t frag=0; frag<q.shape()[1]; ++frag)
                    os<<bait_list[j]
                      <<delim<<prey_prot_list[i]
                      // <<delim<<prey_frags[frag]
                      <<delim<<make_print_column(test_intensities_table_frag_rep[frag],"|")
                      <<delim<<sum_int//rm
                      <<delim<<sum_int/static_cast<double>(test_intensities.size())//rm
                      <<delim<<q.size()
                      <<delim<<make_print_column(ctrl_q[frag],"|")//rm
                      <<delim<<avg_score
                      // <<delim<<max_scores[i][j]
                      // <<delim<<odds(avg_score)
                      <<delim<<bfdr__
                      <<"\n";
        }
}


}