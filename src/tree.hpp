//
//  tree.hpp
//  galax
//
//  Created by Paul O. Lewis on 10/25/14.
//  Copyright (c) 2014 Paul O. Lewis. All rights reserved.
//

#ifndef GALAX_TREE_HPP
#define GALAX_TREE_HPP

#include <boost/shared_ptr.hpp>
#include <iostream>
#include <vector>

namespace galax
{

    template <class T> class TreeManip;

	template <class T>
	class Tree
	{
    friend class TreeManip<T>;

	public:
		Tree();
		virtual             ~Tree();

	private:

		T *                 getRoot();
		void                clear();

		bool                _is_rooted;
		unsigned            _nleaves;
		std::vector<T *>    _preorder;
		std::vector<T>      _nodes; 
  
	public:
		typedef boost::shared_ptr< Tree<T> > TreeShPtr;
	};

	template <class T>
	inline Tree<T>::Tree()
	{
		//std::cerr << "Constructing a Tree" << std::endl;
		clear();
	}

	template <class T>
	inline Tree<T>::~Tree()
	{
		//std::cerr << "Destroying a Tree" << std::endl;
	}

	template <class T>
	inline T * Tree<T>::getRoot()
	{
		if (_preorder.size() > 0)
			return _preorder[0]->_parent;
		else
			return 0;
	}

	template <class T>
	inline void Tree<T>::clear()
	{
		_is_rooted = false;
		_nodes.clear();
		_preorder.clear();
	}

}

#endif //defined(GALAX_TREE_HPP)

