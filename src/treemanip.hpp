//
//  treemanip.hpp
//  galax
//
//  Created by Paul O. Lewis on 10/25/14.
//  Copyright (c) 2014 Paul O. Lewis. All rights reserved.
//

#ifndef GALAX_TREEMANIP_HPP
#define GALAX_TREEMANIP_HPP

#include <boost/format.hpp>
#include <boost/foreach.hpp>
#include <boost/regex.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/foreach.hpp>
#include <boost/range/iterator_range_core.hpp>
#include <set>
#include <numeric>      // accumulate
#include <stdlib.h>     // atoi
#include "split.hpp"
#include "node.hpp"
#include "tree.hpp"
#include "galaxinfo.hpp"
#include "galaxutil.hpp"
#include "xgalax.hpp"

namespace galax
{

typedef std::vector< Split >                      SplitVector;
typedef std::vector< double >                     CountVector;
typedef std::map< SplitVector, CountVector >      CCDMapType;
typedef std::vector< double >                     TimeVector;
typedef std::vector< double >                     WeightVector;
typedef std::map< std::string, double >           CladeInfoMap;

typedef std::vector< SplitVector >                TreeIDType;
typedef std::set< TreeIDType >                    TreeIDSetType;
typedef std::vector< TreeIDSetType >              SubsetTreeSetType;

template <class T>
class TreeManip
    {
	public:
        typedef boost::shared_ptr< TreeManip<T> >       TreeManipShPtr;

                                                        TreeManip() {clear();}
		virtual                                         ~TreeManip() {clear();}
        
        typename Tree<T>::TreeShPtr                     getTree() {return _tree;}
        void                                            setTree(typename Tree<T>::TreeShPtr t) {_tree = t;}
                
        T *                                             rerootAt(int node_index);
        void                                            buildFromNewick(const std::string newick, unsigned root_at);
        void                                            buildFromSplitVector(const std::vector<Split> & splits, unsigned root_at);
        void                                            buildStarTree(unsigned nleaves, unsigned root_at);
		std::string                                     makeNewick(unsigned ndecimals, bool edge_support) const;

        void                                            addToCCDMap(CCDMapType & ccdmap, SubsetTreeSetType & treeCCD, unsigned subset_index, unsigned num_subsets);
        void                                            addToProfile(TimeVector & profile_times, WeightVector & profile_weights, const CladeInfoMap & clade_map);

        std::string                                     debugDescribeNode(T * node) const;
        std::string                                     debugDescribeTree() const;

        void                                            clear() {_tree.reset();}

    protected:

        unsigned                                        countNewickLeaves(const std::string newick);
        void                                            extractNodeNumberFromName(T * nd, std::set<unsigned> & used);
        void                                            refreshPreorder(T * root_node = 0);
        void                                            rerootHelper(T * m, T * t);

    private:

        typename Tree<T>::TreeShPtr                     _tree;
    };


template <class T>
inline std::string TreeManip<T>::debugDescribeNode(T * node) const
    {
    if (!node)
        throw XGalax("Can not describe node because node does not exist (node is NULL)");
        
    std::string s;
    const T & nd = *node;
    std::string node_type = (nd._parent ? (nd._left_child ? "Internal node" : "Tip node") : "Root node");
    
    // Show node name
    if (nd._name.empty())
        s += boost::str(boost::format("\n%s %d:") % node_type % nd._number);
    else
        s += boost::str(boost::format("\n%s \"%s\":") % node_type % nd._name);

    // Show node right sibling
    if (nd._right_sib)
        {
        if (nd._right_sib->_name.empty())
            s += boost::str(boost::format("\n  _right_sib: %d") % nd._right_sib->_number);
        else
            s += boost::str(boost::format("\n  _right_sib: \"%s\"") % nd._right_sib->_name);
        }
    else
        s += "\n  _right_sib: NULL";
        
    // Show node left child
    if (nd._left_child)
        {
        if (nd._left_child->_name.empty())
            s += boost::str(boost::format("\n  _left_child: %d") % nd._left_child->_number);
        else
            s += boost::str(boost::format("\n  _left_child: \"%s\"") % nd._left_child->_name);
        }
    else
        s += "\n  _left_child: NULL";

    // Show node parent
    if (nd._parent)
        {
        if (nd._parent->_name.empty())
            s += boost::str(boost::format("\n  _parent: %d") % nd._parent->_number);
        else
            s += boost::str(boost::format("\n  _parent: \"%s\"") % nd._parent->_name);
        }
    else
        s += "\n  _parent: NULL";

    // Show other info about node
    s += boost::str(boost::format("\n  _number: %d") % nd._number);
    s += boost::str(boost::format("\n  _edge_length: %g") % nd._edge_length);
    s += boost::str(boost::format("\n  _height: %g") % nd._height);
    s += boost::str(boost::format("\n  _edge_support: %s") % nd._edge_support);
    return s;
    }

template <class T>
inline std::string TreeManip<T>::debugDescribeTree() const
    {
    if (!_tree)
        throw XGalax("Can not describe tree because tree does not exist (TreeManip<T>::_tree is NULL)");
        
    std::string s = boost::str(boost::format("Description of the %d nodes in %s tree:\n") % _tree->_nodes.size() % (_tree->_is_rooted ? "rooted" : "unrooted"));;
    s += debugDescribeNode(_tree->getRoot());
    for (typename std::vector<T *>::const_iterator it = _tree->_preorder.begin(); it != _tree->_preorder.end(); ++it)
        {
        s += debugDescribeNode(*it);
        }
        
    return s;
    }

template <class T>
inline void TreeManip<T>::refreshPreorder(T * root_node)
    {
    // Create vector of node pointers in preorder sequence
    _tree->_preorder.clear();   //TODO: why is it necessary to clear?
    _tree->_preorder.reserve(_tree->_nodes.size() - 1); // _preorder does not include root node

	if (!root_node)
		return;

    int min_tip_number = 0;
        
    T * first_preorder = root_node->_left_child;
        
    // sanity check: first preorder node should be the only child of the root node
    assert(first_preorder->_right_sib == 0);
        
    T * nd = first_preorder;
    _tree->_preorder.push_back(nd);
    
    while (true)
        {
        if (!nd->_left_child && !nd->_right_sib)
            {
            // nd has no children (it is a tip) and no siblings, so next preorder is the right sibling of
            // the first ancestral node that has a right sibling.
            if (nd->_number < min_tip_number)
                min_tip_number = nd->_number;

            T * anc = nd->_parent;
            while (anc && !anc->_right_sib)
                anc = anc->_parent;
            if (anc)
                {
                // We found an ancestor with a right sibling
                _tree->_preorder.push_back(anc->_right_sib);
                nd = anc->_right_sib;
                }
            else
                {
                // nd is last preorder node in the tree
                break;
                }
            }
        else if (nd->_right_sib && !nd->_left_child)
            {
            // nd has no children (it is a tip), but does have a sibling on its right
            if (nd->_number < min_tip_number)
                min_tip_number = nd->_number;
            _tree->_preorder.push_back(nd->_right_sib);
            nd = nd->_right_sib;
            }
        else if (nd->_left_child && !nd->_right_sib)
            {
            // nd has children (it is an internal node) but no siblings on its right
            _tree->_preorder.push_back(nd->_left_child);
            nd = nd->_left_child;
            }
        else
            {
            // nd has both children and siblings on its right
            _tree->_preorder.push_back(nd->_left_child);
            nd = nd->_left_child;
            }
            
        }   // end while loop

    // If min_tip_number is less than 0, add -min_tip_number to each tip node number
    // so that the smallest tip number is 0
    assert(min_tip_number == 0 || min_tip_number == -1);
    int offset = -min_tip_number;
        
    // renumber internal nodes in postorder sequence and update _split data members along the way
    int curr_internal = _tree->_nleaves;
    if (_tree->_is_rooted)
        ++curr_internal;

    for (typename std::vector<T *>::reverse_iterator rit = _tree->_preorder.rbegin(); rit != _tree->_preorder.rend(); ++rit)
        {
        T * nd = *rit;

        // reset parent's split if this node is right-most child
        T * parent = nd->_parent;
        assert(parent);
        if (!nd->_right_sib)
            {
            parent->_height = 0.0;
            parent->_split.resetNUnits(_tree->_nleaves);
            }
        
        if (nd->_left_child)
            {
            // nd is an internal node
            nd->_height /= 2.0; // assumes no polytomies
            parent->_height += (nd->_height + nd->_edge_length);
            
            // node numbers for internal nodes begin at _tree->_nleaves (if rooted tree, root node has number _tree->_nleaves)
            nd->_number = curr_internal++;

            // set split weight (edge length)
            nd->_split.setWeight(nd->_edge_length);

            // update parent's split
            parent->_split |= nd->_split;
            }
        else
            {
            // nd is a leaf node
            nd->_number += offset;
            assert(nd->_number >= 0);
            
            parent->_height += (nd->_height + nd->_edge_length);

            // set split
            nd->_split.resetNUnits(_tree->_nleaves);
            nd->_split.setBit(nd->_number);
            
            // update parent's split
            parent->_split.setBit(nd->_number);
            }
        }

    //if (_tree->_is_rooted)
    //    root_node->_number = curr_internal;
    //else
    //    {
    //    // If the first (right-most) bit is "on" in the split for the first preorder node (only
    //    // child of root node), will need to invert all splits
    //    //
    //    //  2     3     4     5      |   1     5     4     3
    //    //   \   /     /     /       |    \   /     /     /
    //    //    \ /     /     /        |     \ /     /     /
    //    //   00110   /     /         |    10001   /     /
    //    //      \   /     /          |       \   /     /
    //    //       \ /     /           |        \ /     /
    //    //      01110   /            |       11001   /
    //    //         \   /             |          \   /
    //    //          \ /              |           \ /
    //    //         11110 <- ok       |          11101 <- needs to be inverted
    //    //           |               |            |
    //    //           |               |            |
    //    //           1               |            2
    //    //
    //    if (_tree->_preorder[0]->_split.isBitSet(0))
    //        {
    //        BOOST_FOREACH(T * nd, _tree->_preorder)
    //            {
    //            nd->_split.invertSplit();
    //            }
    //        }
    //    }
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Makes node `m' (the mover) the rightmost child of node `t' (the target). Assumes both `m' and `t' are non-NULL. 
|   Because the purpose of this function is to move nodes below the new root node to the appropriate point above 
|   the new root node, this method also assumes that `t' is always a descendant (perhaps remote) of `m'.
*/
template <class T>
inline void TreeManip<T>::rerootHelper(T * m, T * t)
	{
	assert(m);
	assert(t);

	// Save nodes to which m attaches
	T * m_lChild    = m->_left_child;
	T * m_rSib      = m->_right_sib;
	T * m_par		= m->_parent;

	// Starting with t, walk down tree to identify x, the child of m that is on the path from m to t
	T * x = t;
	while (x->_parent != m)
		{
		x = x->_parent;
		assert(x != NULL);
		}
	T * x_rSib = x->_right_sib;

	// identify x_lSib, the immediate left sibling of x (will be NULL if x is _left_child of m)
	T * x_lSib = NULL;
	if (x != m_lChild)
		{
		x_lSib = m_lChild;
		while (x_lSib->_right_sib != x)
			{
			x_lSib = x_lSib->_right_sib;
			assert(x_lSib != NULL);
			}
		}

	// identify m_lSib, the immediate left sibling of m (will be NULL if m is root node or is _left_child of its parent)
	T * m_lSib = 0;
	if (m_par && m != m_par->_left_child)
		{
		m_lSib = m_par->_left_child;
		while (m_lSib->_right_sib != m)
			{
			m_lSib = m_lSib->_right_sib;
			assert(m_lSib != NULL);
			}
		}

	// Put x where m is now
	if (m_par == NULL)
		{
		// m is the root node
		assert(m_rSib == NULL);
		assert(m_lSib == NULL);
		x->_right_sib = NULL;
		x->_parent = NULL;
		if (x == m_lChild)
			m->_left_child = x_rSib;
		else
			x_lSib->_right_sib = x_rSib;
		}
	else if (m == m_par->_left_child)
		{
		// m is leftmost child of its parent
		x->_right_sib = m_rSib;
		x->_parent = m_par;
		m->_right_sib = NULL;
		m->_parent = NULL;
		m_par->_left_child = x;
		if (x == m_lChild)
			m->_left_child = x_rSib;
		else
			x_lSib->_right_sib = x_rSib;
		}
	else
		{
		// m is not leftmost child of its parent
		m_lSib->_right_sib = x;
		x->_right_sib = m_rSib;
		x->_parent = m_par;
		m->_right_sib = NULL;
		m->_parent = NULL;
		if (x == m_lChild)
			m->_left_child = x_rSib;
		else
			x_lSib->_right_sib = x_rSib;
		}

	// Make m the new rightmost child of t
	m->_parent = t;
	if (t->_left_child == NULL)
		t->_left_child = m;
	else
		{
		// Find rightmost child of t
		m_lSib = t->_left_child;
		while (m_lSib->_right_sib != NULL)
			m_lSib = m_lSib->_right_sib;

		// Make rightmost child of t the left sib of m
		m_lSib->_right_sib = m;
		}
	}

template <class T>
inline T * TreeManip<T>::rerootAt(int node_index)
	{
	typename std::vector<T>::iterator it = std::find_if(_tree->_nodes.begin(), _tree->_nodes.end(), boost::lambda::bind(&T::_number, boost::lambda::_1) == node_index);
    if (it == _tree->_nodes.end())
        throw XGalax(boost::str(boost::format("no node found with index equal to %d") % node_index));
        
    T & nd = *it;
    if (nd._left_child)
        throw XGalax(boost::str(boost::format("cannot currently root trees at internal nodes (e.g. node %d)") % (nd._number + 1))); //TODO remove this restriction
    
    T * t = &nd;
    T * m = nd._parent;
	while (nd._parent)
		{
		// Begin by swapping the mover's edge length with nd's edge length
        double tmp_edgelen = m->_edge_length;
        std::string tmp_edgesup = m->_edge_support;
        m->_edge_length = nd._edge_length;
        nd._edge_length = tmp_edgelen;
        nd._edge_support = tmp_edgesup;

		// Make the "mover" node m (along with all of its children except nd) the rightmost child of the "target" node t
		rerootHelper(m, t);

		// The next target node is always the previous mover, and the next mover node is always nd's parent
		t = m;
		m = nd._parent;
		}
    return &nd;
    }

template <class T>
inline void TreeManip<T>::extractNodeNumberFromName(T * nd, std::set<unsigned> & used)
	{
    assert(nd);
    std::pair<bool,unsigned> p = convertStringToUnsignedInteger(nd->_name);
    bool success = p.first;
    if (success)
        {
        // conversion succeeded
        if (numberAlreadyUsed(p.second, used))
            throw XGalax(boost::str(boost::format("leaf number %d used more than once") % p.second));
        else
            {
            nd->_number = p.second - 1;
            }
        }
    else
        throw XGalax(boost::str(boost::format("node name (%s) not interpretable as a positive integer") % nd->_name));
    }

template <class T>
inline void TreeManip<T>::addToProfile(TimeVector & profile_times, WeightVector & profile_weights, const CladeInfoMap & clade_map)
	{
    BOOST_FOREACH(T * nd, _tree->_preorder)
        {
        if (nd->_left_child)
            {
            T * a = nd->_left_child;
            T * b = a->_right_sib;

            // Assuming trees are binary (no polytomies), so _left_child should have only one right sibling
            if (b->_right_sib)
                throw XGalax("Expecting all input trees to be binary, but found a polytomy");

            std::string s = nd->_split.createPatternRepresentation();
            CladeInfoMap::const_iterator it = clade_map.find(s);
            if (it != clade_map.end())
                {
                double h = nd->_height;
                profile_times.push_back(h);

                double w = it->second;
                profile_weights.push_back(w);
                }
            //else
            //    {
            //    //temporary!
            //    double h = nd->_height;
            //    profile_times.push_back(h);
            //    profile_weights.push_back(-1.0);
            //    }
            }   // if (nd->_left_child)
        }   // BOOST_FOREACH
    }

template <class T>
inline void TreeManip<T>::addToCCDMap(CCDMapType & ccdmap, SubsetTreeSetType & treeCCD, unsigned subset_index, unsigned num_subsets)
	{
    if (treeCCD.size() == 0)
        {
        //std::cerr << "~~~ Initializing treeCCD ~~~" << std::endl;
        for (unsigned i = 0; i < num_subsets; ++i)
            {
            TreeIDSetType v;
            treeCCD.push_back(v);
            }
        }

    // Add a vector to hold conditional clade definitions (each of which is a SplitVector)
    TreeIDType tree_vector;

    BOOST_FOREACH(T * nd, _tree->_preorder)
        {
        if (nd->_left_child)
            {
            T * a = nd->_left_child;
            T * b = nd->_left_child->_right_sib;

            // Assuming trees are binary (no polytomies), so _left_child should have only one right sibling
            if (b->_right_sib)
                throw XGalax("Expecting all input trees to be binary, but found a polytomy");

            SplitVector v;
            v.push_back(nd->_split);

            // increment the unconditional clade count
            efficientIncrement(ccdmap, v, subset_index, num_subsets);

            if (a->_split < b->_split)
                {
                v.push_back(a->_split);
                v.push_back(b->_split);
                }
            else
                {
                v.push_back(b->_split);
                v.push_back(a->_split);
                }

            // increment the conditional clade count
            efficientIncrement(ccdmap, v, subset_index, num_subsets);

            // add iter to treeCCD for the current tree
            if (nd->_split.countOnBits() > 2)
                tree_vector.push_back(v);

            }   // if (nd->_left_child)
        }   // BOOST_FOREACH


    std::sort(tree_vector.begin(), tree_vector.end());
    treeCCD[subset_index].insert(tree_vector);
        
    }

#define QUICK_AND_DIRTY_NODE_NUMBERS

template <class T>
inline unsigned TreeManip<T>::countNewickLeaves(const std::string newick)
    {
    boost::regex taxonexpr("[(,]\\s*(\\d+|\\S+?|['].+?['])\\s*(?=[,):])");
    boost::sregex_iterator m1(newick.begin(), newick.end(), taxonexpr);
    boost::sregex_iterator m2;
    return (unsigned)std::distance(m1, m2);
    }

template <class T>
inline void TreeManip<T>::buildFromNewick(const std::string newick, unsigned root_at)
	{
    unsigned nleaves_in_newick = countNewickLeaves(newick);
    bool rooted = (root_at == 0);
    unsigned num_nodes = 2*nleaves_in_newick - (rooted ? 0 : 2);
    if (_tree)
        {
        // Just need to perform sanity checks
        assert(root_at > 0 || _tree->_is_rooted);       // if _tree exists, must be already rooted if root_at implies rooted
        assert(_tree->_nleaves == nleaves_in_newick);   // expecting existing tree to have same number of leaves
        assert(_tree->_nodes.size() == num_nodes);      // ensure the number of nodes allocated is correct
        }
    else
        {
        _tree.reset(new Tree<T>());
        _tree->_nleaves = countNewickLeaves(newick);
        if (_tree->_nleaves == 0)
            throw XGalax("Expecting newick tree description to have at least 1 leaf");
        _tree->_nodes.resize(num_nodes);
        }
    _tree->_is_rooted = (root_at == 0);
    _tree->_preorder.assign(_tree->_preorder.size(), 0);

#ifndef QUICK_AND_DIRTY_NODE_NUMBERS
    std::set<unsigned> used; // used to ensure that two tips do not have the same number //TODO eliminate
#endif

    // Assign all nodes a default node number (leaves will replace this number with the number equivalent of their name)
    // and set their
    //TODO is this loop even necessary?
    BOOST_FOREACH(T & nd, _tree->_nodes)
        {
        nd.clear();
        nd._number = _tree->_nleaves;
        }
        
    unsigned curr_leaf = 0;
    unsigned num_edge_lengths = 0;
    unsigned curr_node_index = 0;

	try
		{
		// Root node
		T * nd = &_tree->_nodes[curr_node_index];

		T * root_nd = nd;
        if (_tree->_is_rooted)
            {
            nd = &_tree->_nodes[++curr_node_index];
            nd->_parent = &_tree->_nodes[curr_node_index - 1];
            nd->_parent->_left_child = nd;
            }

		// Some flags to keep track of what we did last
		enum {
			Prev_Tok_LParen		= 0x01,	// previous token was a left parenthesis ('(') 
			Prev_Tok_RParen		= 0x02,	// previous token was a right parenthesis (')') 
			Prev_Tok_Colon		= 0x04,	// previous token was a colon (':') 
			Prev_Tok_Comma		= 0x08,	// previous token was a comma (',') 
			Prev_Tok_Name		= 0x10,	// previous token was a node name (e.g. '2', 'P._articulata') 
			Prev_Tok_EdgeLen	= 0x20	// previous token was an edge length (e.g. '0.1', '1.7e-3') 
			};
		unsigned previous = Prev_Tok_LParen;

		// Some useful flag combinations 
		unsigned LParen_Valid = (Prev_Tok_LParen | Prev_Tok_Comma);
		unsigned RParen_Valid = (Prev_Tok_RParen | Prev_Tok_Name | Prev_Tok_EdgeLen);
		unsigned Comma_Valid  = (Prev_Tok_RParen | Prev_Tok_Name | Prev_Tok_EdgeLen);
		unsigned Colon_Valid  = (Prev_Tok_RParen | Prev_Tok_Name);
		unsigned Name_Valid   = (Prev_Tok_RParen | Prev_Tok_LParen | Prev_Tok_Comma);

		// loop through the characters in newick, building up tree as we go
		std::string::const_iterator newick_start = newick.begin();
		std::string::const_iterator i = newick.begin();
		for (; i != newick.end(); ++i)
			{
			char ch = *i;
			if (iswspace(ch))
				continue;
			switch(ch)
				{
				case ';':
					break;  

				case ')':
					// If nd is bottommost node, expecting left paren or semicolon, but not right paren
					if (!nd->_parent)
						throw XGalax(boost::str(boost::format("Too many right parentheses at position %d in tree description") % (unsigned)std::distance(newick_start, i)));
                        
					// Expect right paren only after an edge length, a node name, or another right paren
					if (!(previous & RParen_Valid))
						throw XGalax(boost::str(boost::format("Unexpected right parenthesisat position %d in tree description") % (unsigned)std::distance(newick_start, i)));

					// Go down a level
					nd = nd->_parent;
					if (!nd->_left_child->_right_sib)
						throw XGalax(boost::str(boost::format("Internal node has only one child at position %d in tree description") % (unsigned)std::distance(newick_start, i)));
					previous = Prev_Tok_RParen;
					break;

				case ':':
					// Expect colon only after a node name or another right paren
					if (!(previous & Colon_Valid))
						throw XGalax(boost::str(boost::format("Unexpected colon at position %d in tree description") % (unsigned)std::distance(newick_start, i)));
					previous = Prev_Tok_Colon;
					break;

				case ',':
					// Expect comma only after an edge length, a node name, or a right paren
					if (!nd->_parent || !(previous & Comma_Valid))
						throw XGalax(boost::str(boost::format("Unexpected comma at position %d in tree description") % (unsigned)std::distance(newick_start, i)));
                        
					// Create the sibling
                    curr_node_index++;
                    if (curr_node_index == _tree->_nodes.size())
						throw XGalax(boost::str(boost::format("Too many nodes specified by tree description (%d nodes allocated for %d leaves)") % _tree->_nodes.size() % _tree->_nleaves));
					nd->_right_sib = &_tree->_nodes[curr_node_index];
					nd->_right_sib->_parent = nd->_parent;
					nd = nd->_right_sib;
					previous = Prev_Tok_Comma;
					break;

				case '(':
					// Expect left paren only after a comma or another left paren
					if (!(previous & LParen_Valid))
						throw XGalax(boost::str(boost::format("Not expecting left parenthesis at position %d in tree description") % (unsigned)std::distance(newick_start, i)));
                        
					// Create new node above and to the left of the current node
					//assert(!nd->_left_child);
                    curr_node_index++;
                    if (curr_node_index == _tree->_nodes.size())
                        throw XGalax(boost::str(boost::format("malformed tree description (more than %d nodes specified)") % _tree->_nodes.size()));
					nd->_left_child = &_tree->_nodes[curr_node_index];
					nd->_left_child->_parent = nd;
					nd = nd->_left_child;
					previous = Prev_Tok_LParen;
					break;

				case '\'':
					// Encountered an apostrophe, which always indicates the start of a 
					// node name (but note that node names do not have to be quoted)

					// Expect node name only after a left paren (child's name), a comma (sib's name)
					// or a right paren (parent's name)
					if (!(previous & Name_Valid))
						throw XGalax(boost::str(boost::format("Not expecting node name at position %d in tree description") % (unsigned)std::distance(newick_start, i)));

					// Get the rest of the name
					nd->_name.clear();
					for (++i; i != newick.end(); ++i)
						{
						ch = *i;
						if (ch == '\'')
							break;
						else if (iswspace(ch))
							nd->_name += ' ';
						else
							nd->_name += ch;
						}
					if (ch != '\'')
						throw XGalax(boost::str(boost::format("Expecting single quote to mark the end of node name at position %d in tree description") % (unsigned)std::distance(newick_start, i)));
                    
					if (!nd->_left_child)
						{
#ifdef QUICK_AND_DIRTY_NODE_NUMBERS
                        int x = atoi(nd->_name.c_str());
                        assert(x > 0);
                        nd->_number = x - 1;
#else
                        extractNodeNumberFromName(nd, used);
#endif
                        curr_leaf++;
						}

					previous = Prev_Tok_Name;
					break;

				default:
					// Expecting either an edge length or an unquoted node name
					if (previous == Prev_Tok_Colon)
						{
						// Edge length expected (e.g. "235", "0.12345", "1.7e-3")
                        // May be multiple edge lengths if created by BayesPhylogenies (e.g. "0.0262473566|0.0262473566")
                        std::string::const_iterator j = i;
                        for (; i != newick.end(); ++i)
							{
							ch = *i;
							if (ch == ',' || ch == ')' || iswspace(ch))
								{
								--i;
								break;
								}
                            else if (ch == '|')
                                {
                                // this version saves only the last edge length when multiple edge lengths are supplied,
                                j = ++i;
                                ch = *i;
                                }
							bool valid = (ch =='e' || ch == 'E' || ch =='.' || ch == '-' || ch == '+' || isdigit(ch));
							if (!valid)
                                throw XGalax(boost::str(boost::format("Invalid branch length character (%c) at position %d in tree description") % ch % (unsigned)std::distance(newick_start, i)));
							}
                        std::string edge_length_str = std::string(j,i+1);
						nd->_edge_length = atof(edge_length_str.c_str());

						++num_edge_lengths;
						previous = Prev_Tok_EdgeLen;
						}
					else
						{
						// Get the node name
						nd->_name.clear();
						for (; i != newick.end(); ++i)
							{
							ch = *i;
							if (ch == '(')
                                throw XGalax(boost::str(boost::format("Unexpected left parenthesis inside node name at position %d in tree description") % (unsigned)std::distance(newick_start, i)));
							if (iswspace(ch) || ch == ':' || ch == ',' || ch == ')')
								{
								--i;
								break;
								}
							nd->_name += ch;
							}

						// Expect node name only after a left paren (child's name), a comma (sib's name) or a right paren (parent's name)
						if (!(previous & Name_Valid))
                            throw XGalax(boost::str(boost::format("Unexpected node name (%s) at position %d in tree description") % nd->_name % (unsigned)std::distance(newick_start, i)));
                            
						if (!nd->_left_child)
							{
#ifdef QUICK_AND_DIRTY_NODE_NUMBERS
                            int x = atoi(nd->_name.c_str());
                            //if (x == 0)
                            //    {
                            //    std::cerr << "problematic tree description: " << nd->_name << std::endl;
                            //    }
                            //assert(x > 0);
                            nd->_number = x - 1;
#else
                            extractNodeNumberFromName(nd, used);
#endif
                            curr_leaf++;
							}

						previous = Prev_Tok_Name;
						}
				}
				if (i == newick.end())
					{
					break;
					}
			}   // loop over characters in newick string

        if (_tree->_is_rooted)
            {
            refreshPreorder(root_nd);
            }
        else
            {
            // Root at leaf whose _number = root_at - 1
            root_nd = rerootAt(root_at - 1);
            refreshPreorder(root_nd);
            }
            
		}
	catch(XGalax x)
		{
		clear();
		throw x;
		}

	}

template <class T>
inline void TreeManip<T>::buildStarTree(unsigned nleaves, unsigned root_at)
	{
    bool rooted = (root_at == 0);
    unsigned num_nodes = 2*nleaves - (rooted ? 0 : 2);
	try
		{
        if (nleaves < 3)
            throw XGalax(boost::str(boost::format("TreeManip<T>::buildStarTree expected splits to contain at least 3 taxa, but found %d instead") % nleaves));

        _tree->clear();

        _tree->_nodes.resize(num_nodes);
        _tree->_is_rooted = rooted;
        _tree->_preorder.assign(_tree->_preorder.size(), 0);

        unsigned curr_node_index = 0;
        unsigned curr_tip_number = 0;
        T * root = &_tree->_nodes[curr_node_index];
        root->_edge_length = 1.0;
        root->_edge_support = "";
        if (!rooted)
            root->_number = root_at - 1;

        T * hub = &_tree->_nodes[++curr_node_index];
        root->_left_child = hub;
        hub->_parent = root;
        hub->_edge_length = 1.0;
        hub->_edge_support = "";
        hub->_number = nleaves;

        for (unsigned i = (rooted ? 0 : 1); i < nleaves; ++i)
            {
            T * nd = &_tree->_nodes[++curr_node_index];
            nd->_right_sib = hub->_left_child;
            nd->_parent = hub;
            nd->_edge_length = 1.0;
            nd->_edge_support = "";
            if (!rooted && curr_tip_number == root_at - 1)
                curr_tip_number++;
            nd->_number = curr_tip_number++;
            hub->_left_child = nd;
            }

        assert(curr_tip_number == nleaves || root_at == nleaves);

        refreshPreorder(root);
        }
	catch(XGalax x)
		{
		clear();
		throw x;
		}
	}

template <class T>
inline void TreeManip<T>::buildFromSplitVector(const std::vector<Split> & split_vect, unsigned root_at)
    {
    bool rooted = (root_at == 0);
    if (split_vect.size() == 0)
        throw XGalax("Tried to build tree from zero-length split vector");

	try
		{
        // Start with clean slate in case _tree already exists
        _tree.reset(new Tree<T>());
        _tree->_is_rooted = rooted;
        _tree->_nleaves = split_vect[0].getNTaxa();
        if (_tree->_nleaves == 0)
            throw XGalax("Expecting splits in supplied split_vect to have at least 1 taxon");
        unsigned max_nodes = 2*_tree->_nleaves - (rooted ? 0 : 2);
        _tree->_nodes.resize(max_nodes);

        // Make a copy of split_vect because we will need to sort splits
        std::vector<Split> splits(split_vect.begin(), split_vect.end());

        // Assume splits are already correctly polarized (i.e. if unrooted, tip that will serve as root
        // is unset in all splits)
        std::vector<Split>::iterator sit = splits.begin();

        // Assume splits have at least 3 taxa
        unsigned ntips = sit->getNTaxa();
        if (ntips < 3)
            throw XGalax(boost::str(boost::format("TreeManip<T>::buildFromSplitVector expected splits to contain at least 3 taxa, but found %d instead") % ntips));

        // Sort the split objects from smallest to largest so that more inclusive splits will
        // follow the less inclusive splits. For example,
        // ---**-*- less inclusive split
        // ---***** more inclusive split
        std::sort(splits.begin(), splits.end());

        // Build a star tree to begin with
        buildStarTree(ntips, root_at);

        // Star tree used ntips + 1 nodes if unrooted, ntips + 2 nodes if rooted
        unsigned curr_node_index = ntips + (rooted ? 1 : 0);

        T * subroot = _tree->_preorder[0];
        T * root = subroot->_parent;

        // Loop over all splits, pulling out taxa specified under a new ancestral node for each
        BOOST_FOREACH(Split & s, splits)
            {
            //std::cerr << boost::str(boost::format("%12.5f %s") % s.getWeight() % s.createPatternRepresentation()) << std::endl;

            // The first split should correspond to the subroot, but we continue here because it makes
            // no sense to detach and then reattach all leaf nodes
            if (subroot->_split.subsumedIn(s))
                continue;

            // Create a new node to hold the taxa in the current split
            assert(curr_node_index + 1 < _tree->_nodes.size());
            T * anc = &_tree->_nodes[++curr_node_index];
            assert(anc);
            anc->_edge_length = 1.0;
            double w = s.getWeight();
            double i = s.getInfo();
            double d = s.getDisparity();
            double c = s.getCertainty();
            std::string str = boost::str(boost::format("P=%.5f I=%.5f D=%.5f IC=%.5f") % w % i % d % c);
            anc->_edge_support = str;

            // Detach nodes in s from tree and add each to anc
            for (typename std::vector<T *>::iterator nit = _tree->_preorder.begin(); nit != _tree->_preorder.end();)
                {
                T * nd = *nit;
                Split & ss = nd->_split;
                T * next = nd->_right_sib;
                if (ss.subsumedIn(s))
                    {
                    // prepare tree for nd removal
                    if (nd == nd->_parent->_left_child)
                        nd->_parent->_left_child = nd->_right_sib;
                    else
                        {
                        T * prev_sib = nd->_parent->_left_child;
                        while (prev_sib->_right_sib != nd)
                            prev_sib = prev_sib->_right_sib;
                        prev_sib->_right_sib = nd->_right_sib;
                        }

                    // move nd to anc
                    nd->_parent = anc;
                    nd->_right_sib = anc->_left_child;
                    anc->_left_child = nd;

                    // skip ahead and continue with nd's right sibling
                    while (nit != _tree->_preorder.end() && (*nit) != next)
                        ++nit;
                    }
                else
                    ++nit;
                }

            // Add anc to subroot (only descendant of tip serving as the root)
            anc->_right_sib = subroot->_left_child;
            anc->_parent = subroot;
            subroot->_left_child = anc;
            refreshPreorder(root);
            }
        }
	catch(XGalax x)
		{
		clear();
		throw x;
		}
	}

// a ************-*
//
//   12345678901234
// b ***-*-*---**-*
// c *---*------*-*
// d *----------*-*
//
//   12345678901234
// e ---*-*-***----
// f ---*---*-*----
//
//   12345678901234
// g -**---*---*---
// h -**-------*---
//
// (((7,(2,3,11)),(5,(1,12,14)))),(6,9,(4,8,10)),13)

template <class T>
inline std::string TreeManip<T>::makeNewick(unsigned ndecimals, bool edge_support) const
    {
	bool rooted = _tree->_is_rooted;
    boost::format edgelen_format(boost::str(boost::format(":%%.%df") % ndecimals));
    boost::format edgelen_support_format(boost::str(boost::format("\"%%s\":%%.%df") % ndecimals));

    std::string s = "(";
    //std::cerr << s << std::endl;
	unsigned open_parens = 1;

	// Start with root node, which may actually represent an extant tip (unrooted trees are
	// rooted at one of the tips)
    typename std::vector<T *>::const_iterator nit = _tree->_preorder.begin();
	const T * nd = *nit;    // nd points to subroot node (node just above root)

	if (!rooted)
		{
        // In this case (unrooted tree), root node is actually a tip
        // output node number (plus 1) as the node name
        s += boost::str(boost::format("%d") % (nd->_parent->_number + 1));

        // In an unrooted tree, where the root node is actually an upside-down extant tip,
        // the root node's edge length is actually held by its only child
        if (ndecimals > 0)
            s += boost::str(edgelen_format % nd->_edge_length);
        //std::cerr << s << std::endl;
		}

    if (!rooted)
        {
	    s += ",";
        //std::cerr << s << std::endl;
        }

	// Now visit all other nodes
	for (++nit; nit != _tree->_preorder.end(); ++nit)
		{
        nd = *nit;

		if (nd->_left_child)
			{
			// nd is internal
			s += "(";
            //std::cerr << s << std::endl;
			++open_parens;
			}
        else
            {
			// nd is a leaf
			// Output tip node number plus 1, then the edge length
            s += boost::str(boost::format("%d") % (nd->_number + 1));
            if (ndecimals > 0)
                s += boost::str(edgelen_format % nd->_edge_length);
            //std::cerr << s << std::endl;

			if (nd->_right_sib)
                {
				s += ",";
                //std::cerr << s << std::endl;
                }
			else
				{
				// Descend toward root until we find an ancestor with a right sibling,
				// outputting edge lengths as we go
				const T * anc = nd;
				while (anc->_parent->_parent && !anc->_right_sib)
					{
					anc = anc->_parent;
					assert(anc);

                    s += ")";
                    if (edge_support && anc->_parent->_parent)
                        s += boost::str(edgelen_support_format % anc->_edge_support % anc->_edge_length);
                    else if (ndecimals > 0)
                        s += boost::str(edgelen_format % anc->_edge_length);
                    //std::cerr << s << std::endl;
                    assert(open_parens > 0);
                    --open_parens;
					}
                if (anc->_parent->_parent)
                    {
                    s += ",";
                    //std::cerr << s << std::endl;
                    }
                else
                    break;
				}
			}
		}

    //std::cerr << "done: " << s << std::endl;
	return s;
	}
}

#endif //defined(GALAX_TREEMANIP_HPP)

