//
//  galaxdata.hpp
//  galax
//
//  Created by Paul O. Lewis on 7/11/15.
//  Copyright (c) 2015 Paul O. Lewis. All rights reserved.
//

#ifndef GALAX_GALAXDATA_HPP
#define GALAX_GALAXDATA_HPP

#define DEBUG_SHOW_ESTIMATED_TREE_POSTERIORS 0

#include <vector>
#include "galaxinfo.hpp"
#include "galaxutil.hpp"
#include "split.hpp"

namespace galax
{

typedef std::vector< Split >                    SplitVector;
typedef std::set< Split >                       SplitSet;
typedef std::vector< double >                   CountVector;
typedef std::vector< unsigned >                 TreeCountsVector;
typedef std::vector< std::string >              NameVector;
typedef std::map< SplitVector, CountVector >    CCDMapType;
typedef std::vector< SplitVector >              TreeIDType;
typedef std::set< TreeIDType >                  TreeIDSetType;
typedef std::vector< TreeIDSetType >            SubsetTreeSetType;

// This class serves to make Galax::estimateInfo less complex by consolidating vectors that store intermediate results
// into a struct that can handle initialization and updating
class GalaxData
	{
	public:
                                    GalaxData(const TreeCountsVector & tree_counts, const NameVector & tree_file_names, unsigned rooted_num_taxa, unsigned outgroup);
                                    ~GalaxData() {}

        void                        newClade(const SplitVector & v, CountVector & c);
        void                        conditionalClade(const SplitVector & v, CountVector & c);
        std::vector<double>         finalizeClade(std::string & detailedinfostr);
        std::pair<unsigned,double>  estimateCoverageForSubset(CCDMapType & ccdmap, SubsetTreeSetType & treeCCD, SubsetTreeMapType & treeMap, std::string & details, unsigned subset_index);
        std::pair<unsigned,double>  estimateMergedCoverage(CCDMapType & ccdmap, SubsetTreeSetType & treeCCD, SubsetTreeMapType & treeMap, std::string & details);
        void                        estimateCoverage(CCDMapType & ccdmap, SubsetTreeSetType & treeCCD, SubsetTreeMapType & treeMap, std::string & details);

        double                      calcNaiveEntropyForSubset(SubsetTreeMapType & treeMap, unsigned subset_index);
        double                      calcNaiveEntropyForMerged(SubsetTreeMapType & treeMap, TreeMapType & m);

        void                        saveDetailedInfoForClade(std::string & detailedinfostr, double D);
        void                        reportTotals(std::string & infostr, std::vector<GalaxInfo> & clade_info);

        double                      getTotalEntropy() const {return _total_entropy;}
        void                        setShowDetails(bool show) {_show_details = show;}

        // Assigned or calculated once in the constructor
        const TreeCountsVector &    _tree_counts;
        const NameVector &          _tree_file_names;
        unsigned                    _num_subsets;
        unsigned                    _total_trees;
        double                      _total_entropy;
        unsigned                    _outgroup;

        // Set using a setter function
        bool                        _show_details;

        // Reassigned or recalculated for each new clade
        Split                       _clade;

        // These vectors store information for one particular clade for each subset plus the merged case
        // All of these vectors are reinitialized for every new clade
        std::vector<double>         _clade_H;
        std::vector<double>         _clade_Hp;
        std::vector<double>         _clade_denom;
        std::vector<double>         _clade_prob;
        std::vector<double>         _I;
        std::vector<double>         _Ipct;

        // These vectors are cumulative across clades and are not reset for each clade.
        std::vector<double>         _total_I;
        double                      _total_D;
        std::vector<double>         _unique;
        std::vector<double>         _coverage;
	};

GalaxData::GalaxData(const TreeCountsVector & tree_counts, const NameVector & tree_file_names, unsigned rooted_num_taxa, unsigned outgroup)
  : _tree_file_names(tree_file_names),
    _tree_counts(tree_counts),
    _total_trees(0),
    _total_entropy(0.0),
    _outgroup(outgroup),
    _total_D(0.0),
    _num_subsets((unsigned)_tree_counts.size()),
    _clade_H(_num_subsets+1, 0.0),
    _clade_Hp(_num_subsets+1, 0.0),
    _clade_denom(_num_subsets+1, 0.0),
    _clade_prob(_num_subsets+1, 0.0),
    _I(_num_subsets+1, 0.0),
    _Ipct(_num_subsets+1, 0.0),
    _total_I(_num_subsets+1, 0.0),
    _unique(_num_subsets+1, 0.0),
    _coverage(_num_subsets+1, 0.0)
    {
    // Determine total number of trees sampled across all subsets
    _total_trees = std::accumulate(_tree_counts.begin(), _tree_counts.end(), 0);

    // Determine maximum possible entropy
    //Split first_clade = (ccdmap.begin()->first)[0];
    //unsigned ntaxa = first_clade.countOnBits() + first_clade.countOffBits();
    //if (!_rooted)
    //    ntaxa -= 1;
    _total_entropy = lognrooted(rooted_num_taxa);
    }

void GalaxData::newClade(const SplitVector & v, std::vector<double> & c)
    {
    assert(v.size() == 1);  // this should be an unconditional clade entry
    _clade = v[0];   // keep track of the current parent clade

    // Individual subsets
    unsigned subset_index = 0;
    BOOST_FOREACH(double subset_count, c)
        {
        _clade_Hp[subset_index] = lognrooted(_clade.countOnBits());
        _clade_H[subset_index] = 0.0;
        _clade_denom[subset_index] = subset_count;
        ++subset_index;
        }

    // Merged case
    double sum_counts = std::accumulate(c.begin(), c.end(), 0.0);
    _clade_Hp[_num_subsets] = lognrooted(_clade.countOnBits());
    _clade_H[_num_subsets] = 0.0;
    _clade_denom[_num_subsets] = sum_counts;
    }

void GalaxData::conditionalClade(const SplitVector & v, std::vector<double> & c)
    {
    assert(v.size() == 3); // parent clade, left clade, right clade
    assert(v[0] == _clade); // parent clade, left clade, right clade

    // Individual subsets
    unsigned subset_index = 0;
    BOOST_FOREACH(double subset_count, c)
        {
        if (subset_count > 0)
            {
            double p = subset_count/_clade_denom[subset_index];
            _clade_H[subset_index]  -= p*log(p);
            _clade_Hp[subset_index] -= p*(lognrooted((v[1].countOnBits())) + lognrooted(v[2].countOnBits()));
            }
        ++subset_index;
        }

    // Merged case
    double sum_counts = std::accumulate(c.begin(), c.end(), 0.0);
    double p = sum_counts/_clade_denom[_num_subsets];
    _clade_H[_num_subsets] -= p*log(p);

    _clade_Hp[_num_subsets] -= p*(lognrooted((v[1].countOnBits())) + lognrooted(v[2].countOnBits()));
    }

void GalaxData::saveDetailedInfoForClade(std::string & detailedinfostr, double D) //, std::string pattern, unsigned num_subsets, unsigned total_trees, std::vector<double> & clade_Hp, std::vector<double> & clade_H, std::vector<double> & w, std::vector<double> & I, std::vector<double> & Ipct, double D)
    {
    const char indiv_template[]  = "%20s %12d %12.5f %12.5f %12.5f %12.5f %12s\n";
    const char merged_template[] = "%20s %12d %12.5f %12.5f %12.5f %12.5f %12.5f\n";
    std::string pattern = _clade.createPatternRepresentation();

    if (_I[_num_subsets] > 0.0)   // don't bother outputting numbers for clades with no information
        {
        detailedinfostr += boost::str(boost::format("\nClade %s\n") % pattern);
        detailedinfostr += boost::str(boost::format("%20s %12s %12s %12s %12s %12s %12s\n")
            % "treefile"
            % "trees"
            % "KL"
            % "p*"
            % "I"
            % "Ipct"
            % "D");
        }
    for (unsigned subset_index = 0; subset_index < _num_subsets; ++subset_index)
        {
        if (_I[subset_index] > 0.0)   // don't bother outputting numbers for clades with no information
            {
            assert(_tree_file_names.size() > 0);
            assert(_tree_counts.size() > 0);
            assert(_clade_Hp.size() > 0);
            assert(_clade_H.size() > 0);
            assert(_clade_prob.size() > 0);
            assert(_I.size() > 0);
            assert(_Ipct.size() > 0);
            double KL_C = _clade_Hp[subset_index] - _clade_H[subset_index];
            detailedinfostr += boost::str(boost::format(indiv_template)
                % _tree_file_names[subset_index]
                % _tree_counts[subset_index]
                % KL_C
                % _clade_prob[subset_index]
                % _I[subset_index]
                % _Ipct[subset_index]
                % "---");
            }
        }
    //if (_num_subsets > 1)
    if (_I[_num_subsets] > 0.0)   // don't bother outputting numbers for clades with no information
        {
        double KL_C = _clade_Hp[_num_subsets] - _clade_H[_num_subsets];
        detailedinfostr += boost::str(boost::format(merged_template)
            % "merged"
            % _total_trees
            % KL_C
            % _clade_prob[_num_subsets]
            % _I[_num_subsets]
            % _Ipct[_num_subsets]
            % D);
        }
    }

std::vector<double> GalaxData::finalizeClade(std::string & detailedinfostr)
    {
    double sum_I = 0.0;
    for (unsigned subset_index = 0; subset_index < _num_subsets; ++subset_index)
        {
        // Compute I for previous clade
        _clade_prob[subset_index] = _clade_denom[subset_index]/_tree_counts[subset_index]; // marginal posterior clade probability

        _I[subset_index] = _clade_prob[subset_index]*(_clade_Hp[subset_index] - _clade_H[subset_index]);
        _total_I[subset_index] += _I[subset_index];

        _Ipct[subset_index] = 100.0*_I[subset_index]/_total_entropy;
        sum_I += _I[subset_index];
        }

    // handle merged case
    // Compute I for previous clade
    _clade_prob[_num_subsets] = _clade_denom[_num_subsets]/_total_trees; // marginal posterior clade probability

    _I[_num_subsets] = _clade_prob[_num_subsets]*(_clade_Hp[_num_subsets] - _clade_H[_num_subsets]);
    _total_I[_num_subsets] += _I[_num_subsets];

    _Ipct[_num_subsets] = 100.0*_I[_num_subsets]/_total_entropy;

    // Compute clade-specific D (note: not percent of maximum D)
    double D = (sum_I/_num_subsets) - _I[_num_subsets];

    saveDetailedInfoForClade(detailedinfostr, D);

    // store numbers for this clade in clade_info vector
    std::vector<double> tmp;
    tmp.push_back(_Ipct[_num_subsets]);
    tmp.push_back(D);
    tmp.push_back(_clade_prob[_num_subsets]);
    tmp.push_back(_I[_num_subsets]);  // added for information profiling
    _total_D += D;

    return tmp;
    }

double GalaxData::calcNaiveEntropyForMerged(SubsetTreeMapType & treeMap, TreeMapType & m)
    {
    // construct merged tree map m and return entropy for merged tree sample

    typedef std::pair<TreeIDType, unsigned> TreeMapPair;

    // first, create merged tree map and count total number of trees sampled
    unsigned total_count = 0;
    unsigned num_subsets = (unsigned)treeMap.size();
    for (unsigned subset_index = 0; subset_index < num_subsets; ++subset_index)
        {
        BOOST_FOREACH(TreeMapPair p, treeMap[subset_index])
            {
            efficientAddTo(m, p.first, p.second);
            total_count += p.second;
            }
        }
    double log_total_count = log(total_count);

    // now calculate naive entropy
    double naive_entropy = 0.0;
    BOOST_FOREACH(TreeMapPair p, m)
        {
        double prob = double(p.second)/total_count;
        double log_prob = log(p.second) - log_total_count;
        naive_entropy -= prob*log_prob;
        }
    return naive_entropy;
    }

double GalaxData::calcNaiveEntropyForSubset(SubsetTreeMapType & treeMap, unsigned subset_index)
    {
    // calculate for specific subset indicated

    // first calculate total number of trees sampled
    unsigned total_count = 0;
    typedef std::pair<TreeIDType, unsigned> TreeMapPair;
    BOOST_FOREACH(TreeMapPair p, treeMap[subset_index])
        {
        total_count += p.second;
        }
    double log_total_count = log(total_count);

    // now calculate naive entropy
    double naive_entropy = 0.0;
    BOOST_FOREACH(TreeMapPair p, treeMap[subset_index])
        {
        double prob = double(p.second)/total_count;
        double log_prob = log(p.second) - log_total_count;
        naive_entropy -= prob*log_prob;
        }
    return naive_entropy;
    }

std::pair<unsigned,double> GalaxData::estimateCoverageForSubset(CCDMapType & ccdmap, SubsetTreeSetType & treeCCD, SubsetTreeMapType & treeMap, std::string & details, unsigned subset_index)
    {
    TreeManip<Node>::TreeManipShPtr tm(new TreeManip<Node>());
    SplitSet nontrivial_splits;

    double max_log_product = 0.0;
    std::vector< double > log_products;

    const TreeIDSetType & tree_set = treeCCD[subset_index];
    unsigned num_unique_trees = (unsigned)tree_set.size();

    if (_show_details)
        {
        details += boost::str(boost::format("Distinct tree topologies from file %s:\n\n") % _tree_file_names[subset_index]);
        details += "Key to columns:\n";
        details += "  tree      = (arbitrary) index of a specific tree topology\n";
        details += "  count     = number of trees sampled corresponding to topology\n";
        details += "  log(prob) = log of the probability of topology computed using conditional clade distribution\n";
        details += "  newick    = Newick-formatted tree description of topology\n\n";
        details += boost::str(boost::format("%20s %12s %12s  %s\n") % "tree" % "count" % "log(prob)" % "newick");
        }

    unsigned tree_number = 1;
    for (TreeIDSetType::iterator tree_iter = tree_set.begin(); tree_iter != tree_set.end(); ++tree_iter, ++tree_number)
        {
        const TreeIDType & tree_id = *tree_iter;

        nontrivial_splits.clear();
        double log_prob = 0.0;

        double log_product = 0.0;
        for (TreeIDType::const_iterator ccd_iter = tree_id.begin(); ccd_iter != tree_id.end(); ++ccd_iter)
            {
            const SplitVector & v = *ccd_iter;
            const CountVector & cv = ccdmap[v];
            double numer_count = cv[subset_index];
            SplitVector mother;
            mother.push_back(v[0]);
            double denom_count = ccdmap[mother][subset_index];

            log_product += log(numer_count);
            log_product -= log(denom_count);

            if (_show_details)
                {
                if (v[0].countOnBits() > 1 && v[0].countOffBits() > 1)
                    nontrivial_splits.insert(v[0]);
                if (v[1].countOnBits() > 1 && v[1].countOffBits() > 1)
                    nontrivial_splits.insert(v[1]);
                if (v[2].countOnBits() > 1 && v[2].countOffBits() > 1)
                    nontrivial_splits.insert(v[2]);
                log_prob += log(numer_count);
                log_prob -= log(denom_count);
                }
            }
        if (log_product > max_log_product)
            max_log_product = log_product;
        log_products.push_back(log_product);

        if (_show_details)
            {
            //std::cerr << "***** Calling TreeManip::buildFromSplitVector with _outgroup = "<< _outgroup << std::endl;
            tm->buildFromSplitVector(SplitVector(nontrivial_splits.begin(), nontrivial_splits.end()), _outgroup);
            unsigned count = treeMap[subset_index][tree_id];
            std::string newick = tm->makeNewick(0, false);
            details += boost::str(boost::format("%20d %12d %12.5f  %s\n") % tree_number % count % log_prob % newick);
            }
        }

    // compute coverage for this subset
    double sum_exp_diffs = 0.0;
    BOOST_FOREACH(double logprod, log_products)
        {
        sum_exp_diffs += exp(logprod - max_log_product);
        }
    double log_coverage = max_log_product + log(sum_exp_diffs);
    double exp_log_coverage = exp(log_coverage);
    if (exp_log_coverage - 1.0 > 1.e-8)
        {
        std::cerr << "*** this can't be right! exp_log_coverage = " << exp_log_coverage << std::endl;
        }

    if (_show_details)
        {
        double naive_entropy = calcNaiveEntropyForSubset(treeMap, subset_index);
        details += boost::str(boost::format("\nFrom tree file %s: unique = %d  naive entropy = %.5f coverage = %.5f\n") % _tree_file_names[subset_index] % num_unique_trees % naive_entropy % exp_log_coverage);
        } 

    return std::pair<unsigned,double>(num_unique_trees, exp_log_coverage);
    }

std::pair<unsigned,double> GalaxData::estimateMergedCoverage(CCDMapType & ccdmap, SubsetTreeSetType & treeCCD, SubsetTreeMapType & treeMap, std::string & details)
    {
    // This function is very similar to GalaxData::estimateCoverageForSubset, but decided to
    // separate it out because of the need to construct the merged tree set and the fact that
    // ccdmap stores only counts for individual subsets, not the merged case, so the sums must
    // be done ad hoc. All this would lead to a very difficult-to-understand combined function.

    TreeManip<Node>::TreeManipShPtr tm(new TreeManip<Node>());
    SplitSet nontrivial_splits;

    double max_log_product = 0.0;
    std::vector< double > log_products;

    // Compute merged tree set
    TreeIDSetType tree_set;
    for (unsigned i = 0; i < _num_subsets; ++i)
        {
        TreeIDSetType & treeset_subseti = treeCCD[i];
        for (TreeIDSetType::iterator treeid_iter = treeset_subseti.begin(); treeid_iter != treeset_subseti.end(); ++treeid_iter)
            {
            const TreeIDType & tree_id = (*treeid_iter);
            tree_set.insert(tree_id);
            }
        }

    unsigned num_unique_trees = (unsigned)tree_set.size();

    TreeMapType merged_tree_map;
    double naive_entropy = 0.0;
    if (_show_details)
        {
        details += "Distinct tree topologies in merged tree sample\n\n";
        details += "Key to columns:\n";
        details += "  tree      = (arbitrary) index of a specific tree topology\n";
        details += "  count     = number of trees sampled corresponding to topology\n";
        details += "  log(prob) = log of the probability of topology computed using conditional clade distribution\n";
        details += "  newick    = Newick-formatted tree description of topology\n\n";
        details += boost::str(boost::format("%20s %12s %12s  %s\n") % "tree" % "count" % "log(prob)" % "newick");
        naive_entropy = calcNaiveEntropyForMerged(treeMap, merged_tree_map);
        }

    unsigned tree_number = 1;
    for (TreeIDSetType::iterator tree_iter = tree_set.begin(); tree_iter != tree_set.end(); ++tree_iter, ++tree_number)
        {
        const TreeIDType & tree_id = *tree_iter;

        nontrivial_splits.clear();
        double log_prob = 0.0;

        double log_product = 0.0;
        for (TreeIDType::const_iterator ccd_iter = tree_id.begin(); ccd_iter != tree_id.end(); ++ccd_iter)
            {
            const SplitVector & v = *ccd_iter;
            const CountVector & cv = ccdmap[v];
            double numer_count = std::accumulate(cv.begin(), cv.end(), 0.0);

            SplitVector mother;
            mother.push_back(v[0]);
            const CountVector & cvmother = ccdmap[mother];
            double denom_count = std::accumulate(cvmother.begin(), cvmother.end(), 0.0);

            log_product += log(numer_count);
            log_product -= log(denom_count);

            if (_show_details)
                {
                nontrivial_splits.insert(v[0]);
                if (v[1].countOnBits() > 1)
                    nontrivial_splits.insert(v[1]);
                if (v[2].countOnBits() > 1)
                    nontrivial_splits.insert(v[2]);
                log_prob += log(numer_count);
                log_prob -= log(denom_count);
                }
            }

        if (log_product > max_log_product)
            max_log_product = log_product;
        log_products.push_back(log_product);

        if (_show_details)
            {
            tm->buildFromSplitVector(SplitVector(nontrivial_splits.begin(), nontrivial_splits.end()), _outgroup);
            unsigned count = merged_tree_map[tree_id];
            std::string newick = tm->makeNewick(0, false);
            details += boost::str(boost::format("%20d %12d %12.5f  %s\n") % tree_number % count % log_prob % newick);
            }
        }

    // compute coverage for this subset
    double sum_exp_diffs = 0.0;
    BOOST_FOREACH(double logprod, log_products)
        {
        sum_exp_diffs += exp(logprod - max_log_product);
        }
    double log_coverage = max_log_product + log(sum_exp_diffs);
    double exp_log_coverage = exp(log_coverage);
    if (exp_log_coverage - 1.0 > 1.e-8)
        {
        std::cerr << "*** this can't be right! exp_log_coverage = " << exp_log_coverage << std::endl;
        }

    if (_show_details)
        {
        details += boost::str(boost::format("\nFrom merged: unique = %d  naive entropy = %.5f coverage = %.5f\n") % num_unique_trees % naive_entropy % exp_log_coverage);
        }

    return std::pair<unsigned,double>(num_unique_trees, exp_log_coverage);
    }
    
void GalaxData::estimateCoverage(CCDMapType & ccdmap, SubsetTreeSetType & treeCCD, SubsetTreeMapType & treeMap, std::string & details)
    {
    std::pair<unsigned,double> unicov;

    // Individual subsets
    for (unsigned subset_index = 0; subset_index < _num_subsets; ++subset_index)
        {
        if (_show_details)
            details += boost::str(boost::format("\nTree topology posterior probabilities for tree file %d\n") % (subset_index + 1));
        unicov = estimateCoverageForSubset(ccdmap, treeCCD, treeMap, details, subset_index);
        _unique[subset_index] = unicov.first;
        _coverage[subset_index] = unicov.second;
        }

    // Merged case
    if (_show_details)
        details += "\nTree topology posterior probabilities for the merged tree file\n";
    unicov = estimateMergedCoverage(ccdmap, treeCCD, treeMap, details);
    _unique[_num_subsets] = unicov.first;
    _coverage[_num_subsets] = unicov.second;
    }

std::ostream & operator<<(std::ostream & os, const GalaxInfo & info)
  {
  os << "\n\n";
  BOOST_FOREACH(double d, info._value)
    {
    os << "  " << d;
    }
  os << std::endl;
  return os;
  }

void GalaxData::reportTotals(std::string & infostr, std::vector<GalaxInfo> & clade_info)
    {
    // Report totals for each subset
    std::vector<GalaxInfo> subset_info;
    for (unsigned subset_index = 0; subset_index < _num_subsets; ++subset_index)
        {
        std::vector<double> tmp;
        tmp.push_back(_total_I[subset_index]);
        tmp.push_back(_unique[subset_index]);
        tmp.push_back(_coverage[subset_index]);
        subset_info.push_back(GalaxInfo(_tree_file_names[subset_index], tmp));
        }
    //GalaxInfo::_sortby_index = 8;  // 0=Ipct, 1=unique, 2=coverage //wtf? 8 is the length of the column header "treefile"
    //std::sort(subset_info.begin(), subset_info.end(), std::greater<GalaxInfo>());

    unsigned max_length_treefile_name = 0;
    BOOST_FOREACH(GalaxInfo & c, subset_info)
        {
        if (c._name.size() > max_length_treefile_name)
            max_length_treefile_name = (unsigned)c._name.size();
        }

    // Create format strings to use with boost::format
    std::string column_headers = boost::str(boost::format("%%%ds %%12s %%12s %%12s %%12s %%12s %%12s %%12s %%12s\n") % max_length_treefile_name);
    std::string indiv_summary  = boost::str(boost::format("%%%ds %%12d %%12.5f %%12.5f %%12.5f %%12.5f %%12.5f %%12s %%12s\n") % max_length_treefile_name);
    std::string merged_summary = boost::str(boost::format("%%%ds %%12d %%12.5f %%12.5f %%12.5f %%12.5f %%12.5f %%12.5f %%12.5f\n") % max_length_treefile_name);

    // Create summary table showing coverage, entropy, and information for individual subsets
    // as well as averages across subsets and the corresponding results for the merged sample
    infostr += "\nTotals\n\n";
    infostr += "Key to columns:\n";
    infostr += "  unique   = number of distinct tree topologies\n";
    infostr += "  coverage = Larget (2013. Syst. Biol. 62:501–511) estimated posterior coverage\n";
    infostr += "  H        = entropy of marginal prior tree topology distribution\n";
    infostr += "  H*       = entropy of marginal posterior tree topology distribution\n";
    infostr += "  I        = Lindley information (H - H*)\n";
    infostr += "  Ipct     = I expressed as percent of maximum (100*I/H)\n";
    infostr += "  D        = dissonance (merged H* - average H*)\n";
    infostr += "  Dpct     = D expressed as percent of maximum (100*D/(merged H*))\n\n";
    infostr += boost::str(boost::format(column_headers)
        % "treefile"
        % "unique"
        % "coverage"
        % "H"
        % "H*"
        % "I"
        % "Ipct"
        % "D"
        % "Dpct");
    unsigned subset_index = 0;

    double average_unique = 0.0;
    double average_coverage = 0.0;
    double average_posterior_entropy = 0.0;
    double average_info = 0.0;
    double average_info_pct = 0.0;
    double prior_entropy = _total_entropy;

    BOOST_FOREACH(GalaxInfo & c, subset_info)
        {
        double info = c._value[0];
        double unique = c._value[1];
        double coverage = c._value[2];
        double posterior_entropy = prior_entropy - info;
        double pct = 100.0*info/_total_entropy;

        average_unique += unique;
        average_coverage += coverage;
        average_posterior_entropy += posterior_entropy;
        average_info += info;
        average_info_pct += pct;

        infostr += boost::str(boost::format(indiv_summary)
            % c._name
            % int(unique)
            % coverage
            % prior_entropy
            % posterior_entropy
            % info
            % pct
            % "---"
            % "---");
        ++subset_index;
        }

    if (subset_index > 0)
        {
        average_posterior_entropy /= subset_index;
        average_info /= subset_index;
        average_unique /= subset_index;
        average_coverage /= subset_index;
        average_info_pct /= subset_index;
        }

    // Line for averages
    infostr += boost::str(boost::format(indiv_summary)
        % "average"
        % average_unique
        % average_coverage
        % prior_entropy
        % average_posterior_entropy
        % average_info
        % average_info_pct
        % "---"
        % "---");

    // Sort clades for merged case by Ipct, D, and w and save majrule tree
    double totalIpct = 0.0;
    double cumpct = 0.0;
    double unique = _unique[_num_subsets];
    double coverage = _coverage[_num_subsets];
    if (_num_subsets > 1)
        {
        double Imerged = _total_I[_num_subsets];
        totalIpct = (100.0*Imerged/_total_entropy);
        double merged_entropy = prior_entropy - Imerged;
        double D = merged_entropy - average_posterior_entropy;
        double Dpct = 0.0;
        if (merged_entropy > 0.0)
            {
            Dpct = 100.0*D/merged_entropy;
            }
        infostr += boost::str(boost::format(merged_summary)
            % "merged"
            % int(unique)
            % coverage
            % prior_entropy
            % merged_entropy
            % _total_I[_num_subsets]
            % totalIpct
            % D
            % Dpct);

        // Uncomment next line to output sum of clade specific D estimates to compare with overall estimate output on merged line
        //infostr += boost::str(boost::format("\n_total_D = %.8f\n") % _total_D);

        // Report merged info for each clade sorted from highest to lowest
        GalaxInfo::_sortby_index = 0; // 0=Ipct, 1=D, 2=w, 3=I
        std::sort(clade_info.begin(), clade_info.end(), std::greater<GalaxInfo>());
        const char clade_summary[]  = "%12.5f %12.5f %12.5f %12.5f %12.5f %s\n";
        infostr += std::string("\nClades sorted by merged info (top 95% shown):\n");
        infostr += boost::str(boost::format("%12s %12s %12s %12s %12s %s\n")
            % "I"
            % "Ipct"
            % "cum. Ipct"
            % "D"
            % "w"
            % "clade");
        cumpct = 0.0;
        BOOST_FOREACH(GalaxInfo & c, clade_info)
            {
            double info = c._value[0];
            double diff = c._value[1];
            double w    = c._value[2];
            double pct  = 100.0*info/totalIpct;
            cumpct += pct;
            if (cumpct > 95)
                break;
            infostr += boost::str(boost::format(clade_summary)
                % info
                % pct
                % cumpct
                % diff
                % w
                % c._name);
            }

        // Report diff for each clade sorted from highest to lowest
        GalaxInfo::_sortby_index = 1; // 0=Ipct, 1=D, 2=w, 3=I
        std::sort(clade_info.begin(), clade_info.end(), std::greater<GalaxInfo>());
        infostr += std::string("\nClades sorted by D (top 95% shown):\n");
        infostr += boost::str(boost::format("%12s %12s %12s %12s %12s %s\n")
            % "D"
            % "Dpct"
            % "cum. Dpct"
            % "Ipct"
            % "w"
            % "clade");
        cumpct = 0.0;
        BOOST_FOREACH(GalaxInfo & c, clade_info)
            {
            double info = c._value[0];
            double diff = c._value[1];
            double w    = c._value[2];
            double pct  = 0.0;
            if (_total_D > 0.0)
                pct = 100.0*diff/_total_D;
            cumpct += pct;
            if (cumpct > 95)
                break;
            infostr += boost::str(boost::format(clade_summary)
                % diff
                % pct
                % cumpct
                % info
                % w
                % c._name);
            }
        }   //if (num_subsets > 1)

    // Report clade posterior for each clade sorted from highest to lowest
    GalaxInfo::_sortby_index = 2;   // 0=Ipct, 1=D, 2=w, 3=I
    std::sort(clade_info.begin(), clade_info.end(), std::greater<GalaxInfo>());
    infostr += std::string("\nClades sorted by merged clade posterior (w) (only those >= 50% shown):\n");
    infostr += boost::str(boost::format("%12s %12s %12s %s\n")
        % "w"
        % "Ipct"
        % "D"
        % "clade");
    cumpct = 0.0;
    //double log2 = log(2.0);
    BOOST_FOREACH(GalaxInfo & c, clade_info)
        {
        double info = c._value[0];
        double diff = c._value[1];
        double w    = c._value[2];

        if (w >= 0.5)
            {
            infostr += boost::str(boost::format("%12.5f %12.5f %12.5f %s\n")
                % w
                % info
                % diff
                % c._name);
            }
        }

    }
}

#endif

