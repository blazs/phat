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
    typedef std::pair<index, zp> _entry_t;
    typedef std::vector< _entry_t > _f_column_t;

  
}



