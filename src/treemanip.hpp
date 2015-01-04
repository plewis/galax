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

template <class T>
class TreeManip
    {
	public:
                                        TreeManip() {clear();}
		virtual                         ~TreeManip() {clear();}
        
        typename Tree<T>::TreeShPtr     getTree() {return _tree;}
        void                            setTree(typename Tree<T>::TreeShPtr t) {_tree = t;}
                
        T *                             rerootAt(int node_index);
        void                            buildFromNewick(const std::string newick, bool rooted = false);

        void                            addToCCDMap(unsigned subset_index, unsigned num_subsets);
        void                            showCCDMap(unsigned subset_index);
        void                            estimateMergedInfo(std::string & infostr, std::vector<unsigned> & tree_counts, std::vector< std::string > & treefile_names);

        std::string                     debugDescribeNode(T * node) const;
        std::string                     debugDescribeTree() const;

        void                            clear() {_tree.reset();}

    protected:

        void                            extractNodeNumberFromName(T * nd, std::set<unsigned> & used);
        void                            refreshPreorder(T * root_node = 0);
        void                            rerootHelper(T * m, T * t);

    public:
        typedef boost::shared_ptr< TreeManip<T> >       TreeManipShPtr;
        typedef std::vector< Split >                    SplitVector;
        typedef std::vector< double >                   CountVector;
        typedef std::map< SplitVector, CountVector >    CCDMapType;

        static const unsigned                           ALLSUBSETS;

    private:

        typename Tree<T>::TreeShPtr _tree;
        CCDMapType                  _ccdmap;
    };


template<class T>
const unsigned TreeManip<T>::ALLSUBSETS = (unsigned)-1;

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
    _tree->_preorder.clear();
    _tree->_preorder.reserve(_tree->_nodes.size() - 1); // _preorder does not include root node

	if (!root_node)
		return;
        
    T * first_preorder = root_node->_left_child;
        
    // sanity check: first preorder node should be the only child of the root node
    assert(first_preorder->_right_sib == 0);
        
    T * nd = first_preorder;
    _tree->_preorder.push_back(nd);
    
    while (true)
        {
        if (!nd->_left_child && !nd->_right_sib)
            {
            // nd has no children and no siblings, so next preorder is the right sibling of 
            // the first ancestral node that has a right sibling.
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
        
    // renumber internal nodes in postorder sequence and update _split data members along the way
    int curr_internal = _tree->_nleaves;
    for (typename std::vector<T *>::reverse_iterator rit = _tree->_preorder.rbegin(); rit != _tree->_preorder.rend(); ++rit)
        {
        T * nd = *rit;

        // reset parent's split if this node is right-most child
        T * parent = nd->_parent;
        assert(parent);
        if (!nd->_right_sib)
            {
            parent->_split.resetNUnits(_tree->_nleaves);
            }
        
        if (nd->_left_child)
            {
            // nd is an internal node
            
            // node numbers for internal nodes are negative integers
            nd->_number = curr_internal++;

            // update parent's split
            parent->_split |= nd->_split;
            }
        else
            {
            // nd is a leaf node
            
            // set split
            nd->_split.resetNUnits(_tree->_nleaves);
            nd->_split.setBit(nd->_number);
            
            // update parent's split
            parent->_split.setBit(nd->_number);
            }
        }

    if (_tree->_is_rooted)
        root_node->_number = curr_internal;
    else
        {
        // If the first (right-most) bit is "on" in the split for the first preorder node (only
        // child of root node), will need to invert all splits before adding them to treestats:
        //
        //  2     3     4     5      |   1     5     4     3
        //   \   /     /     /       |    \   /     /     /
        //    \ /     /     /        |     \ /     /     /
        //   00110   /     /         |    10001   /     /
        //      \   /     /          |       \   /     /
        //       \ /     /           |        \ /     /
        //      01110   /            |       11001   /
        //         \   /             |          \   /
        //          \ /              |           \ /
        //         11110 <- ok       |          11101 <- needs to be inverted
        //           |               |            |
        //           |               |            |
        //           1               |            2
        //
        if (_tree->_preorder[0]->_split.isBitSet(0))
            {
            BOOST_FOREACH(T * nd, _tree->_preorder)
                {
                nd->_split.invertSplit();
                }
            }
        }
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
        throw XGalax(boost::str(boost::format("cannot currently root trees at internal nodes (e.g. node %d)") % nd._number)); //TODO remove this restriction
    
    T * t = &nd;
    T * m = nd._parent;
	while (nd._parent)
		{
		// Begin by swapping the mover's edge length with nd's edge length
        double tmp_edgelen = m->_edge_length;
        m->_edge_length = nd._edge_length;
        nd._edge_length = tmp_edgelen;

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
inline void TreeManip<T>::estimateMergedInfo(std::string & s, std::vector<unsigned> & tree_counts, std::vector< std::string > & treefile_names)
	{
    assert(_ccdmap.size() > 0);

    unsigned num_subsets = (unsigned)tree_counts.size();
    unsigned total_trees = std::accumulate(tree_counts.begin(), tree_counts.end(), 0);

    // Determine maximum possible entropy
    Split first_clade = (_ccdmap.begin()->first)[0];
    unsigned ntaxa = first_clade.countOnBits() + first_clade.countOffBits();
    if (!_tree->_is_rooted)
        ntaxa -= 1;
    double total_entropy = lognrooted(ntaxa);

    s = "";
    s += "\nLindley Information\n";
    s += boost::str(boost::format("  number of taxa: %d\n") % ntaxa);
    s += boost::str(boost::format("  total prior entropy: %.5f\n") % total_entropy);

    Split clade;

    std::vector<GalaxInfo> clade_info;

    std::vector<double> clade_H(num_subsets+1, 0.0);
    std::vector<double> clade_Hp(num_subsets+1, 0.0);
    std::vector<double> clade_denom(num_subsets+1, 0.0);
    std::vector<double> w(num_subsets+1, 0.0);
    std::vector<double> I(num_subsets+1, 0.0);
    std::vector<double> Ipct(num_subsets+1, 0.0);
    std::vector<double> total_I(num_subsets+1, 0.0);

    double total_D = 0.0;
    bool first = true;
    const char indiv_template[]  = "%20s %12d %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12s\n";
    const char merged_template[] = "%20s %12d %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f\n";
    for (CCDMapType::iterator it = _ccdmap.begin(); it != _ccdmap.end(); ++it)
        {
        // v contains either 1 split (unconditional) or 3 splits (conditional)
        const SplitVector & v = it->first;

        std::vector<double> & count = it->second;
        double sum_counts = std::accumulate(count.begin(), count.end(), 0.0);

        // Depends on _ccdmap keys (split vectors) being sorted such that an entry for an unconditional clade (e.g. ABC)
        // precedes any conditional clade entries (e.g. A|BC) that have that clade as parent. Thus, when next unconditional
        // clade entry is encountered, we know that there are no more conditional clade entries for the previous parent
        // clade and we can at that point compute the information for the previous parent.
        if (first)
            {
            first = false;

            // Initialize data for first clade
            assert(v.size() == 1);  // first entry should be an unconditional clade entry
            clade = v[0];

            //std::string clade_pattern_representation = clade.createPatternRepresentation();
            //s += boost::str(boost::format("\nClade %s\n") % clade_pattern_representation);

            unsigned subset_index = 0;
            BOOST_FOREACH(double subset_count, count)
                {
                clade_Hp[subset_index] = lognrooted(clade.countOnBits());
                clade_H[subset_index] = 0.0;
                clade_denom[subset_index] = subset_count;
                ++subset_index;
                }

            // Merged case
            clade_Hp[num_subsets] = lognrooted(clade.countOnBits());
            clade_H[num_subsets] = 0.0;
            clade_denom[num_subsets] = sum_counts;
            }
        else
            {
            if (v[0] == clade)
                {
                // Conditional clade entry whose parent is clade
                assert(v.size() == 3); // parent clade, left clade, right clade

                // Individual subsets
                unsigned subset_index = 0;
                BOOST_FOREACH(double subset_count, count)
                    {
                    if (subset_count > 0)
                        {
                        double p = subset_count/clade_denom[subset_index];
                        clade_H[subset_index] -= p*log(p);
                        clade_Hp[subset_index] -= p*(lognrooted((v[1].countOnBits())) + lognrooted(v[2].countOnBits()));
                        }
                    ++subset_index;
                    }

                // Merged case
                double p = sum_counts/clade_denom[num_subsets];
                clade_H[num_subsets] -= p*log(p);
                clade_Hp[num_subsets] -= p*(lognrooted((v[1].countOnBits())) + lognrooted(v[2].countOnBits()));
                }
            else
                {
                // Just found next unconditional clade entry, so now is the time to finish computing
                // information content for the previous clade

                assert(v.size() == 1);
                unsigned subset_index = 0;
                double sum_Ipct = 0.0;
                BOOST_FOREACH(double subset_count, count)
                    {
                    // Compute I for previous clade
                    w[subset_index] = clade_denom[subset_index]/tree_counts[subset_index]; // marginal posterior clade probability

                    I[subset_index] = w[subset_index]*(clade_Hp[subset_index] - clade_H[subset_index]);
                    total_I[subset_index] += I[subset_index];

                    Ipct[subset_index] = 100.0*I[subset_index]/total_entropy;
                    sum_Ipct += Ipct[subset_index];

                    // Initialize data for next clade
                    clade_Hp[subset_index] = lognrooted(clade.countOnBits());
                    clade_H[subset_index] = 0.0;
                    clade_denom[subset_index] = subset_count;

                    ++subset_index;
                    }

                // handle merged case
                // Compute I for previous clade
                w[num_subsets] = clade_denom[num_subsets]/total_trees; // marginal posterior clade probability

                I[num_subsets] = w[num_subsets]*(clade_Hp[num_subsets] - clade_H[num_subsets]);
                total_I[num_subsets] += I[num_subsets];

                Ipct[num_subsets] = 100.0*I[num_subsets]/total_entropy;

                // Compute clade-specific D
                double D = (sum_Ipct/num_subsets) - Ipct[num_subsets];
                total_D += D;

                if (I[num_subsets] > 0.0)
                    {
                    s += boost::str(boost::format("\nClade %s\n") % clade.createPatternRepresentation());
                    s += boost::str(boost::format("%20s %12s %12s %12s %12s %12s %12s %12s %12s\n")
                        % "treefile"
                        % "trees"
                        % "Hp"
                        % "H"
                        % "Hp - H"
                        % "w"
                        % "I"
                        % "Ipct"
                        % "D");
                    for (subset_index = 0; subset_index < num_subsets; ++subset_index)
                        {
                        s += boost::str(boost::format(indiv_template)
                            % treefile_names[subset_index]
                            % tree_counts[subset_index]
                            % clade_Hp[subset_index]
                            % clade_H[subset_index]
                            % (clade_Hp[subset_index] - clade_H[subset_index])
                            % w[subset_index]
                            % I[subset_index]
                            % Ipct[subset_index]
                            % "---");
                        }
                    if (num_subsets > 1)
                        {
                        s += boost::str(boost::format(merged_template)
                            % "merged"
                            % total_trees
                            % clade_Hp[num_subsets]
                            % clade_H[num_subsets]
                            % (clade_Hp[num_subsets] - clade_H[num_subsets])
                            % w[num_subsets]
                            % I[num_subsets]
                            % Ipct[num_subsets]
                            % D);
                        std::vector<double> tmp;
                        tmp.push_back(Ipct[num_subsets]);
                        tmp.push_back(D);
                        clade_info.push_back(GalaxInfo(clade.createPatternRepresentation(), tmp));
                        }
                    }

                // Initialize data for next clade
                clade = v[0];

                subset_index = 0;
                BOOST_FOREACH(double subset_count, count)
                    {
                    clade_Hp[subset_index] = lognrooted(clade.countOnBits());
                    clade_H[subset_index] = 0.0;
                    clade_denom[subset_index] = subset_count;
                    ++subset_index;
                    }

                // Merged case
                clade_Hp[num_subsets] = lognrooted(clade.countOnBits());
                clade_H[num_subsets] = 0.0;
                clade_denom[num_subsets] = sum_counts;
                }
            }
        }

    // Compute I for final clade
    double sum_Ipct = 0.0;
    for (unsigned subset_index = 0; subset_index < num_subsets; ++subset_index)
        {
        if (clade_denom[subset_index] > 0)
            {
            w[subset_index] = clade_denom[subset_index]/tree_counts[subset_index];

            I[subset_index] = w[subset_index]*(clade_Hp[subset_index] - clade_H[subset_index]);
            total_I[subset_index] += I[subset_index];

            Ipct[subset_index] = 100.0*I[subset_index]/total_entropy;
            sum_Ipct += Ipct[subset_index];
            }
        }

    w[num_subsets] = clade_denom[num_subsets]/total_trees;

    I[num_subsets] = w[num_subsets]*(clade_Hp[num_subsets] - clade_H[num_subsets]);
    total_I[num_subsets] += I[num_subsets];

    Ipct[num_subsets] = 100.0*I[num_subsets]/total_entropy;

    // Compute clade-specific D
    double D = (sum_Ipct/num_subsets) - Ipct[num_subsets];
    total_D += D;

    if (I[num_subsets] > 0.0)
        {
        s += boost::str(boost::format("\nClade %s\n") % clade.createPatternRepresentation());
        s += boost::str(boost::format("%20s %12s %12s %12s %12s %12s %12s %12s %12s\n")
            % "treefile"
            % "trees"
            % "Hp"
            % "H"
            % "Hp - H"
            % "w"
            % "I"
            % "Ipct"
            % "D");
        for (unsigned subset_index = 0; subset_index < num_subsets; ++subset_index)
            {
            s += boost::str(boost::format(indiv_template)
                % treefile_names[subset_index]
                % tree_counts[subset_index]
                % clade_Hp[subset_index]
                % clade_H[subset_index]
                % (clade_Hp[subset_index] - clade_H[subset_index])
                % w[subset_index]
                % I[subset_index]
                % Ipct[subset_index]
                % "---");
            }
        if (num_subsets > 1)
            {
            s += boost::str(boost::format(merged_template)
                % "merged"
                % total_trees
                % clade_Hp[num_subsets]
                % clade_H[num_subsets]
                % (clade_Hp[num_subsets] - clade_H[num_subsets])
                % w[num_subsets]
                % I[num_subsets]
                % Ipct[num_subsets]
                % D);
                std::vector<double> tmp;
                tmp.push_back(Ipct[num_subsets]);
                tmp.push_back(D);
                clade_info.push_back(GalaxInfo(clade.createPatternRepresentation(), tmp));
            }
        }

    // Report totals for each subset
    std::vector<GalaxInfo> subset_info;
    for (unsigned subset_index = 0; subset_index < num_subsets; ++subset_index)
        {
        std::vector<double> tmp;
        tmp.push_back(total_I[subset_index]);
        subset_info.push_back(GalaxInfo(treefile_names[subset_index], tmp));
        }
    GalaxInfo::_sortby_index = 0;
    std::sort(subset_info.begin(), subset_info.end(), std::greater<GalaxInfo>());

    const char indiv_summary[]  = "%20s %12.5f %12.5f %12s\n";
    const char merged_summary[] = "%20s %12.5f %12.5f %12.5f\n";
    s += std::string("\nTotals\n");
    s += boost::str(boost::format("%20s %12s %12s %12s\n")
        % "treefile"
        % "I"
        % "Ipct"
        % "D");
    BOOST_FOREACH(GalaxInfo & c, subset_info)
        {
        double info = c._value[0];
        double pct = 100.0*info/total_entropy;
        s += boost::str(boost::format(indiv_summary)
            % c._name
            % info
            % pct
            % "---");
        }
    double totalIpct = 0.0;
    if (num_subsets > 1)
        {
        totalIpct = (100.0*total_I[num_subsets]/total_entropy);
        s += boost::str(boost::format(merged_summary)
            % "merged"
            % total_I[num_subsets]
            % (100.0*total_I[num_subsets]/total_entropy)
            % total_D);
        }

    // Report merged info for each clade sorted from highest to lowest
    GalaxInfo::_sortby_index = 0;
    std::sort(clade_info.begin(), clade_info.end(), std::greater<GalaxInfo>());
    const char clade_summary[]  = "%12.5f %12.5f %12.5f %12.5f %s\n";
    s += std::string("\nClades sorted by merged info (top 95% shown):\n");
    s += boost::str(boost::format("%12s %12s %12s %12s %s\n")
        % "I"
        % "%"
        % "cum. %"
        % "D"
        % "clade");
    double cumpct = 0.0;
    BOOST_FOREACH(GalaxInfo & c, clade_info)
        {
        double info = c._value[0];
        double diff = c._value[1];
        double pct = 100.0*info/totalIpct;
        cumpct += pct;
        if (cumpct > 95)
            break;
        s += boost::str(boost::format(clade_summary)
            % info
            % pct
            % cumpct
            % diff
            % c._name);
        }

    // Report diff for each clade sorted from highest to lowest
    GalaxInfo::_sortby_index = 1;
    std::sort(clade_info.begin(), clade_info.end(), std::greater<GalaxInfo>());
    s += std::string("\nClades sorted by D (top 95% shown):\n");
    s += boost::str(boost::format("%12s %12s %12s %12s %s\n")
        % "D"
        % "%"
        % "cum. %"
        % "Ipct"
        % "clade");
    cumpct = 0.0;
    BOOST_FOREACH(GalaxInfo & c, clade_info)
        {
        double info = c._value[0];
        double diff = c._value[1];
        double pct = 100.0*diff/total_D;
        cumpct += pct;
        if (cumpct > 95)
            break;
        s += boost::str(boost::format(clade_summary)
            % diff
            % pct
            % cumpct
            % info
            % c._name);
        }

    }

template <class T>
inline void TreeManip<T>::showCCDMap(unsigned subset_index)
	{
    std::cerr << "CCD map has " << _ccdmap.size() << " elements." << std::endl;
    for (CCDMapType::iterator it = _ccdmap.begin(); it != _ccdmap.end(); ++it)
        {
        const SplitVector & v = it->first;
        double count = 0.0;
        if (subset_index < ALLSUBSETS)
            count = it->second[subset_index];
        else
            count = std::accumulate(it->second.begin(), it->second.end(), 0.0);
        std::cerr << "\n+--------------------------" << std::endl;
        for (SplitVector::const_iterator vit = v.begin(); vit != v.end(); ++vit)
            {
            std::cerr << "| " << vit->createPatternRepresentation() << std::endl;
            }
        std::cerr << "| " << count << std::endl;
        }
    }

template <class T>
inline void TreeManip<T>::addToCCDMap(unsigned subset_index, unsigned num_subsets)
	{
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
            efficientIncrement(_ccdmap, v, subset_index, num_subsets);

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
            efficientIncrement(_ccdmap, v, subset_index, num_subsets);

            }   // if (nd->_left_child)
        }   // BOOST_FOREACH
    }

#define QUICK_AND_DIRTY_NODE_NUMBERS

template <class T>
inline void TreeManip<T>::buildFromNewick(const std::string newick, bool rooted)
	{
    // Assume that if _tree already exists, it has the correct number of nodes
    if (!_tree)
        {
        _tree.reset(new Tree<T>());
        _tree->_nleaves = countNewickLeaves(newick);
        if (_tree->_nleaves == 0)
            throw XGalax("Expecting newick tree description to have at least 1 leaf");
        unsigned max_nodes = 2*_tree->_nleaves - (rooted ? 0 : 2);
        _tree->_nodes.resize(max_nodes);
        }
    _tree->_is_rooted = rooted;
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
        
	// This will point to the first tip node encountered so that we can reroot at this node before returning
	T * first_tip = 0;  //TODO: pass in root tip number, this needs to be consistent
    
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
						if (!first_tip)
							first_tip = nd;
						}

					previous = Prev_Tok_Name;
					break;

				default:
					// Expecting either an edge length or an unquoted node name
					if (previous == Prev_Tok_Colon)
						{
						// Edge length expected (e.g. "235", "0.12345", "1.7e-3")
                        std::string::const_iterator j = i;
                        for (; i != newick.end(); ++i)
							{
							ch = *i;
							if (ch == ',' || ch == ')' || iswspace(ch))
								{
								--i;
								break;
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
                            assert(x > 0);
                            nd->_number = x - 1;
#else
                            extractNodeNumberFromName(nd, used);
#endif
                            curr_leaf++;
							if (!first_tip)
								first_tip = nd;
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
            // Root at leaf whose _number = 0
            root_nd = rerootAt(0);
            refreshPreorder(root_nd);
            }
            
		}
	catch(XGalax x)
		{
		clear();
		throw x;
		}

	}
}

#endif //defined(GALAX_TREEMANIP_HPP)
