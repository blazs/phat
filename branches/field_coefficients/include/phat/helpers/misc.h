/*  Copyright 2013 IST Austria
    Contributed by: Ulrich Bauer, Michael Kerber, Jan Reininghaus

    This file is part of PHAT.

    PHAT is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    PHAT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with PHAT.  If not, see <http://www.gnu.org/licenses/>. */

#pragma once

// STL includes
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <list>
#include <map>
#include <algorithm>
#include <queue>
#include <cassert>
#include <sstream>
#include <algorithm>
#include <iomanip>
#include <cmath>
#include <cstdlib>

#ifndef ZP_PRIME
#define ZP_PRIME 11
#endif

// VS2008 and below unfortunately do not support stdint.h
#if defined(_MSC_VER)&& _MSC_VER < 1600
    typedef __int8 int8_t;
    typedef unsigned __int8 uint8_t;
    typedef __int16 int16_t;
    typedef unsigned __int16 uint16_t;
    typedef __int32 int32_t;
    typedef unsigned __int32 uint32_t;
    typedef __int64 int64_t;
    typedef unsigned __int64 uint64_t;
#else
    #include <stdint.h>
#endif

 
// basic types. index can be changed to int32_t to save memory on small instances
namespace phat {
    typedef int8_t zp;
    typedef int64_t index;
    typedef int8_t dimension;
    typedef std::vector< index > _column_t;

    class _f_column_t{
    protected:
        std::vector< std::pair<index,zp> > data;
    public:
        typedef std::vector< std::pair<index,zp> >::iterator iterator;
        typedef std::vector< std::pair<index,zp> >::reverse_iterator reverse_iterator;
        
        index& operator[] (const index x) {
            return data[x].first;
        }
        void push_back(const index x){
            data.push_back(std::pair<index,zp>(x,1));
        }
        void resize(const index num_rows ){
            data.resize(num_rows);
        }
        index size(){
            return data.size();
        }
        void clear(){
            data.clear();
        }
        iterator  begin(){
            return data.begin();
        }
        iterator  end(){
            return data.end();
        }
        
        reverse_iterator  rbegin(){
            return data.rbegin();
        }
        reverse_iterator  rend(){
            return data.rend();
        }
        

    };
}

// OpenMP (proxy) functions
#if defined _OPENMP
    #include <omp.h>
#else
    #define omp_get_thread_num() 0
    #define omp_get_max_threads() 1
    #define omp_get_num_threads() 1
	void omp_set_num_threads( int ) {};
    #include <time.h>
    #define omp_get_wtime() (float)clock() / (float)CLOCKS_PER_SEC
#endif

#include <phat/helpers/thread_local_storage.h>



