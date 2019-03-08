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
/*
drop in replacement for boost::multi_array.
rationale:
- boost::multi_array's copy assignment operator requires both objects to have the same dimensions
- and you have to manually resize the destination object if that is not satisfied
- but sometimes you may want to assign to a multi_array of different dimension
- for example, when initializing a vector of multi_array s.
  std::vector<boost::multi_array<int,2>> v(100);
  for(size_t i=0; i<100; ++i) v[i]=f(i);
*/
#include <boost/multi_array.hpp>

#include <boost/container/vector.hpp>

#include <array>

namespace bmarrw{

template<typename T,std::size_t NumDims,typename Allocator = ::std::allocator<T>>
class bmarr;


template<typename T,std::size_t NumDims,typename Allocator>
class bmarr{
public:
    typedef typename boost:: template multi_array_ref<T,NumDims> ref_type;

// from multi_array_ref.hpp
    typedef typename ref_type::value_type value_type;
    typedef typename ref_type::reference reference;
    typedef typename ref_type::iterator iterator;
    typedef typename ref_type::reverse_iterator reverse_iterator;
    typedef typename ref_type::const_reference const_reference;
    typedef typename ref_type::const_iterator const_iterator;
    typedef typename ref_type::const_reverse_iterator const_reverse_iterator;
    typedef typename ref_type::element element;
    typedef typename ref_type::size_type size_type;
    typedef typename ref_type::difference_type difference_type;
    typedef typename ref_type::index index;
    typedef typename ref_type::extent_range extent_range;

    typedef typename ref_type::storage_order_type storage_order_type;
    typedef typename ref_type::index_list index_list;
    typedef typename ref_type::size_list size_list;

    template <std::size_t NDims>
    struct const_array_view {
        typedef typename ref_type::template const_array_view<NDims>::type type;
    };
    template <std::size_t NDims>
    struct array_view {
        typedef typename ref_type::template array_view<NDims>::type type;
    };
private:
    boost::container::vector<T,Allocator> vec;
    // multi_array_ref does not have external resources
    ref_type vec_ref;

    template<typename Tint>
    std::array<std::uintmax_t,NumDims> make_array(const Tint* arr) const
    {
        std::array<std::uintmax_t,NumDims> ret;
        for(std::size_t i=0;i<ret.size();++i)
            ret[i]=arr[i];
        return ret;
    }
    std::array<std::uintmax_t,NumDims> make_array() const
    {
        std::array<std::uintmax_t,NumDims> ret;
        ret.fill(0);
        return ret;
    }

public:
    bmarr():vec{},vec_ref(NULL,make_array()) {}//default ctor

    // template <typename Tarr>
    // bmarr(Tarr const& other)://ctor
    bmarr(boost::const_multi_array_ref<T,NumDims> const& other):
        vec(other.num_elements()),
        vec_ref(vec.data(),make_array(other.vec_ref.shape()))
    {}


    // template <class ExtentList>
    // bmarr(ExtentList const& extents):
    ///// ctor for objects with default ctor
    template <class Tint, template <typename,std::size_t> class Tarr>
    bmarr(Tarr<Tint,NumDims> const& extents):
        vec(std::accumulate(extents.begin(),extents.end(),1,
                            std::multiplies<std::uintmax_t>{})
            ,boost::container::default_init
            ),
        vec_ref(vec.data(),extents)
    {static_assert(std::is_integral<Tint>::value,"indexing array must have integral elements");}


    ///// ctor for objects without default ctor
    template <class Tint, template <typename,std::size_t> class Tarr, class Trange>
    bmarr(Tarr<Tint,NumDims> const& extents, const Trange& vec2):
        vec(vec2.begin(),vec2.end()),
        vec_ref(vec.data(),extents)
    {
        static_assert(std::is_integral<Tint>::value,"indexing array must have integral elements");
        const auto n=std::accumulate(
            extents.begin(),extents.end(),1,
            std::multiplies<std::uintmax_t>{});
        if(static_cast<size_t>(n)!=vec.size())
            throw std::runtime_error{"array dims!=vector size"};
    }

    bmarr(bmarr const& other)://copy ctor
        vec(other.vec),
        vec_ref(vec.data(),make_array(other.vec_ref.shape()))
    {}

    bmarr(bmarr&& other)://move ctor
        vec(std::move(other.vec)),vec_ref(other.vec_ref)
    {new (&vec_ref) ref_type(other.vec_ref);}

    bmarr& operator=(bmarr const& other)
    {
        vec=other.vec;
        new (&vec_ref) ref_type(vec.data(),make_array(other.vec_ref.shape()));
        return *this;
    }

    bmarr& operator=(bmarr&& other)
    {
        vec=std::move(other.vec);
        new (&vec_ref) ref_type(other.vec_ref);
        return *this;
    }

    boost::multi_array_types::size_type
    size() const {return vec_ref.size();}
    const boost::multi_array_types::size_type*
    shape() const {return vec_ref.shape();}

    iterator begin(){return vec_ref.begin();}
    iterator end(){return vec_ref.end();}
    const_iterator begin()const{return vec_ref.begin();}
    const_iterator end()const{return vec_ref.end();}

    reference operator[](index idx)
    {return vec_ref[idx];}
    const_reference operator[](index idx) const
    {return vec_ref[idx];}

    template <int NDims>
    typename array_view<NDims>::type
    operator[](const boost::detail::multi_array::
               index_gen<NumDims,NDims>& indices)
    {return vec_ref[indices];}
    template <int NDims>
    typename const_array_view<NDims>::type
    operator[](const boost::detail::multi_array::
               index_gen<NumDims,NDims>& indices) const
    {return vec_ref[indices];}


    template <class IndexList>
    element& operator()(const IndexList& indices)
    {return vec_ref(indices);}
    template <class IndexList>
    const element& operator()(const IndexList& indices) const
    {return vec_ref(indices);}


    // template <class ExtentList>
    // void resize(ExtentList const& extents)
    template <class Tint, template <typename,std::size_t> class Tarr>
    void resize(Tarr<Tint,NumDims> const& extents)
    {
        static_assert(std::is_integral<Tint>::value,"indexing array must have integral elements");
        vec.resize(std::accumulate(extents.begin(),extents.end(),1,
                                   std::multiplies<std::uintmax_t>{})
                   ,boost::container::default_init);
        ::new(static_cast<void*>(&vec_ref)) ref_type(vec.data(),extents);
    }

    std::size_t num_elements() const {return vec.size();}

    const T* data() const {return vec.data();}

    T* data() {return vec.data();}

};

}