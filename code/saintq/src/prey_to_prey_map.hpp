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


// #include <string>
// #include <vector>
// namespace saint4{


std::vector<std::vector<std::size_t> > get_prey_to_prey_map(
    const std::vector<str_t>& prot_list,
    std::istream& is)
{
    using namespace ::std;
    namespace str=boost::algorithm;
    const char* const delim="\t";
    const auto nprot=prot_list.size();

    {
        ///// check for duplicates in protein list
        std::unordered_map<str_t,unsigned> tmp;
        for(size_t i=0;i<nprot;++i)
            ++tmp[prot_list[i]];
        BOOST_FOREACH(const auto& e , tmp){
            if(e.second!=1)
                // throw std::runtime_error("duplicated preys"+ boost::convert<string>(e,cnv).value());
                throw std::runtime_error("duplicated preys"+ e.first);
        }
    }

    std::unordered_map<str_t,size_t> prot_index;
    for(size_t i=0;i<nprot;++i)
        prot_index[prot_list[i]]=i;

    vector<vector<str_t>> rows;// rows of cells
    {
        vector<str_t> ss;// temp storage for a row of cells
        str_t line;
        while(getline(is,line))
            rows.push_back(str::split(ss, line, str::is_any_of(delim)));
    }

    bmarrw::bmarr<unsigned,2> p2p_2d_arr{i__(nprot,nprot)};
    // iRefWeb
    for(size_t r=0;r<rows.size();++r)/*try*/{
        const auto& row=rows[r];
        if(row.size()!=2)throw std::runtime_error("number of columns!=2");
        if(r==0) continue;// skip first row
        const auto i=prot_index.at(row[0]);
        const auto j=prot_index.at(row[1]);
        ++p2p_2d_arr(i__(i,j));
    }//catch (const std::out_of_range& e){cerr<<e.what()<<"\n";}

    std::vector<std::vector<size_t> > ret(nprot);
    for(size_t i=0;i<nprot;++i)
        for(size_t j=0;j<nprot;++j)
            if(p2p_2d_arr(i__(i,j))>0)
                ret[i].push_back(j);

    return ret;
}


// }