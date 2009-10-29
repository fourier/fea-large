#ifndef __STD_H__
#define __STD_H__

// Precompiled headers file

//
// STL includes sections
// 
#include <cmath>
#include <memory>
#include <iostream>
#include <fstream>
#include <vector>
#include <cassert>
#include <stdexcept>
#include <cmath>
#include <valarray>
#include <map>
#include <list>
#include <algorithm>
#include <string>
#include <sstream>
#include <stdexcept>

//
// Boost includes section
//
// Boost defines
//
// Boost documentation stands for:
// "...by default, the library constructs a global extent_gen object boost::extents. 
// In case of concern about memory used by these objects, defining BOOST_MULTI_ARRAY_NO_GENERATORS 
// before including the library header inhibits its construction."
// So let's define it:
#define BOOST_MULTI_ARRAY_NO_GENERATORS 1
// to avoid panics with fatal error LNK1104: cannot open file 'libboost_filesystem-vc71-sgd-1_35.lib'
// while trying to link against boost::filesystem
#define BOOST_ALL_NO_LIB 1 

// boost headers
#include <boost/static_assert.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/utility.hpp>
#include <boost/multi_array.hpp>
//#include <boost/filesystem.hpp>

// Lua headers
#include <lua.hpp>

// other headers
#include "global.h"
#include "jobtype.h"

#endif // __STD_H__
