//
//  node.cpp
//  galax
//
//  Created by Paul O. Lewis on 10/25/14.
//  Copyright (c) 2014 Paul O. Lewis. All rights reserved.
//

#ifndef GALAX_NODE_HPP
#define GALAX_NODE_HPP

#include <string>
#include "split.hpp"

namespace galax
{
    template <class T> class Tree;
    template <class T> class TreeManip;

    class Node
	{

    friend class Tree<Node>;
    friend class TreeManip<Node>;

	public:
		Node();
		virtual             ~Node();

	protected:

		virtual void        clear();

		Node *              _left_child;
		Node *              _right_sib;
		Node *              _parent;
		int                 _number;
		std::string         _name;
		double              _edge_length;
		double              _edge_support;
        Split               _split;
	};

	inline Node::Node()
	{
		//std::cout << "Creating Node object" << std::endl;
		clear();
	}

	inline Node::~Node()
	{
		//std::cout << "Destroying Node object" << std::endl;
	}

	inline void Node::clear()
	{
		_left_child = 0;
		_right_sib = 0;
		_parent = 0;
		_number = 0;
		_edge_length = 0.0;
		_edge_support = 0.0;
	}

}

#endif /* defined(GALAX_NODE_HPP) */

