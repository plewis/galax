//
//  cladeinfo.hpp
//  galax
//
//  Created by Paul O. Lewis on 1/4/15.
//  Copyright (c) 2015 Paul O. Lewis. All rights reserved.
//

#ifndef GALAX_CLADEINFO_HPP
#define GALAX_CLADEINFO_HPP

#include "split.hpp"

namespace galax
{

class CladeInfo
	{
	public:
							CladeInfo() {}
							CladeInfo(Split & s, double v) : _split(s), _value(v) {}
							~CladeInfo() {}

        bool                operator<(const CladeInfo & other) const {return (bool)(_value < other._value);}
        bool                operator>(const CladeInfo & other) const {return (bool)(_value > other._value);}

		Split               _split;
        double              _value;
	};

}

#endif
