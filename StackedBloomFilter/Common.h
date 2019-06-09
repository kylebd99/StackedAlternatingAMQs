#pragma once
#include <cmath>
#include <iostream>
#include <sstream>
#include <cstdint>
#include <vector>
#include <functional>
#include <algorithm>
#include <numeric>
#include "windows.h"
#include "CityHash.h"

#define DBOUT( s )            \
{                             \
   std::wostringstream os_;    \
   os_ << s;                   \
   OutputDebugStringW( os_.str().c_str() );  \
}