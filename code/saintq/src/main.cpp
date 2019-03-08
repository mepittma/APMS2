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


#include "main.hpp"


int main(const int argc,const char*const*const argv) try {

    util::check_init_global_streams();
    if(argc!=2){
        std::cerr<<"Usage:	saintq <parameter file>\n";
        return EXIT_FAILURE;
    }
    ///// read param file
    std::ifstream param_istream{argv[1]};
    if(!param_istream.good()){
        std::cerr<<"Cannot open file: “"<<argv[1]<<"”"<<std::endl;
        return EXIT_FAILURE;
    }

    util::check_init_istream(param_istream);
    const auto param=::saint4::read_param(param_istream, std::cerr);
    param_istream.close();


    const auto input_lvl = param.input_lvl;
    const ::std::string& input_filename=param.input_filename;

    std::ifstream ifs(input_filename);
    util::check_init_istream(ifs);
    auto s1=::saint4::read_table1(ifs,param,std::cerr);
    ifs.close();
    auto s2=::saint4::read_table2(std::move(s1),std::cerr,param);


    if(0){
        std::ifstream ifs2("PPP2CA_newinput2_prot_level.tsv");
        util::check_init_istream(ifs2);
        const auto prey_to_prey_map=::saint4::get_prey_to_prey_map(s2.prey_prot_list,ifs2);
    }


    const auto idx_info=::saint4::get_indexing(s2,input_lvl,std::cerr);
    auto model=calc(s2, idx_info, input_lvl,::std::cerr, param.compress_n_ctrl,
                    param.normalize_control);
    const auto prot_prot_scores=scoring(
        std::move(model),
        idx_info,
        param,
        std::cerr);

    if(false){
        std::ofstream ofs{"scores_list_DEBUG__"+input_filename+"__.tsv"};
        ::saint4::output(ofs,
                         prot_prot_scores,
                         s2,
                         idx_info,
                         input_lvl,
                         ::saint4::output_level::protein,
                         true/*more output*/);
    }
    {
        std::ofstream ofs{"scores_list__"+input_filename+"__.tsv"};
        ::saint4::output(ofs,
                         prot_prot_scores,
                         s2,
                         idx_info,
                         input_lvl,
                         ::saint4::output_level::protein,
                         false/*more output*/);
    }


}catch(const std::exception& e){
    const char* ename = typeid(e).name();
    std::cout<<"\n"<<ename<<"\t"<<e.what()<<"\n";
    return EXIT_FAILURE;
}