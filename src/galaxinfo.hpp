//
//  galaxinfo.hpp
//  galax
//
//  Created by Paul O. Lewis on 1/4/15.
//  Copyright (c) 2015 Paul O. Lewis. All rights reserved.
//

#ifndef GALAX_GALAXINFO_HPP
#define GALAX_GALAXINFO_HPP

namespace galax
{

class GalaxInfo
	{
	public:
							GalaxInfo() {}
							GalaxInfo(const std::string & s, std::vector<double> & v) : _name(s), _value(v) {}
							~GalaxInfo() {}

        bool                operator<(const GalaxInfo & other) const {return (bool)(_value[_sortby_index] < other._value[_sortby_index]);}
        bool                operator>(const GalaxInfo & other) const {return (bool)(_value[_sortby_index] > other._value[_sortby_index]);}

		std::string         _name;

        // _value for clades:  0=Ipct, 1=D, 2=w, 3=I
        // _value for subsets: 0=Ipct, 1=unique, 2=coverage
        std::vector<double> _value;


        static unsigned     _sortby_index;
	};

}

#endif

