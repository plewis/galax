//
//  xgalax.hpp
//  galax
//
//  Created by Paul O. Lewis on 10/25/14.
//  Copyright (c) 2014 Paul O. Lewis. All rights reserved.
//

#ifndef GALAX_XGALAX_HPP
#define GALAX_XGALAX_HPP

namespace galax
{

class XGalax : public std::exception
	{
	public:
							XGalax() throw() {}
							XGalax(const std::string s) throw() : _msg() {_msg = s;}
		virtual				~XGalax() throw() {}
		const char *        what() const throw() {return _msg.c_str();}
        
    private:
    
		std::string         _msg;
	}; 

}

#endif
