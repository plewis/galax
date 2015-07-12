//
//  galaxdata.hpp
//  galax
//
//  Created by Paul O. Lewis on 7/11/15.
//  Copyright (c) 2015 Paul O. Lewis. All rights reserved.
//

#ifndef GALAX_GALAXDATA_HPP
#define GALAX_GALAXDATA_HPP

#include <vector>
#include "galaxinfo.hpp"
#include "galaxutil.hpp"
#include "split.hpp"

namespace galax
{

typedef std::vector< Split >                    SplitVector;
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
                                    GalaxData(const TreeCountsVector & tree_counts, const NameVector & tree_file_names, unsigned ntaxa);
                                    ~GalaxData() {}

        void                        saveDetailedInfoForClade(std::string & detailedinfostr, double D);
        void                        newClade(const SplitVector & v, CountVector & c);
        void                        conditionalClade(const SplitVector & v, CountVector & c);
        std::vector<double>         finalizeClade(std::string & detailedinfostr);
        std::pair<unsigned,double>  estimateCoverage(CCDMapType & ccdmap, SubsetTreeSetType & treeCCD, unsigned subset_index);
        void                        calcCoverage(CCDMapType & ccdmap, SubsetTreeSetType & treeCCD);
        void                        reportTotals(std::string & infostr, std::vector<GalaxInfo> & clade_info);

        double                      getTotalEntropy() const {return _total_entropy;}

        const NameVector &          _tree_file_names;
        const TreeCountsVector &    _tree_counts;
        unsigned                    _total_trees;
        double                      _total_entropy;
        double                      _total_D;

        // Reassigned or recalculated for each new clade
        double                      _sum_counts;
        unsigned                    _num_subsets;
        Split                       _clade;

        // These vectors store information for one particular clade for each subset plus the merged case
        // All of these vectors are reinitialized for every new clade: not currently storing all information
        // for every subset-clade combination
        std::vector<double>         _clade_H;
        std::vector<double>         _clade_Hp;
        std::vector<double>         _clade_denom;
        std::vector<double>         _clade_prob;
        std::vector<double>         _I;
        std::vector<double>         _Ipct;

        // These vectors hold cumulative I across clades, number of unique tree topologies, and coverage for each subset plus the merged case
        // These are not reset for each clade.
        std::vector<double>         _total_I;
        std::vector<double>         _unique;
        std::vector<double>         _coverage;
	};

GalaxData::GalaxData(const TreeCountsVector & tree_counts, const NameVector & tree_file_names, unsigned ntaxa)
  : _tree_file_names(tree_file_names),
    _tree_counts(tree_counts),
    _total_trees(0),
    _total_entropy(0.0),
    _total_D(0.0),
    _sum_counts(0.0),
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
    _total_entropy = lognrooted(ntaxa);
    }

void GalaxData::newClade(const SplitVector & v, std::vector<double> & c)
    {
    _sum_counts = std::accumulate(c.begin(), c.end(), 0.0);
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
    _clade_Hp[_num_subsets] = lognrooted(_clade.countOnBits());
    _clade_H[_num_subsets] = 0.0;
    _clade_denom[_num_subsets] = _sum_counts;
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
    double p = _sum_counts/_clade_denom[_num_subsets];
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
            % "D_KL"
            % "w"
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
            double D_KL = _clade_Hp[subset_index] - _clade_H[subset_index];
            detailedinfostr += boost::str(boost::format(indiv_template)
                % _tree_file_names[subset_index]
                % _tree_counts[subset_index]
                % D_KL
                % _clade_prob[subset_index]
                % _I[subset_index]
                % _Ipct[subset_index]
                % "---");
            }
        if (_num_subsets > 1)
            {
            double D_KL = _clade_Hp[_num_subsets] - _clade_H[_num_subsets];
            detailedinfostr += boost::str(boost::format(merged_template)
                % "merged"
                % _total_trees
                % D_KL
                % _clade_prob[_num_subsets]
                % _I[_num_subsets]
                % _Ipct[_num_subsets]
                % D);
            }
        }
    }

std::vector<double> GalaxData::finalizeClade(std::string & detailedinfostr)
    {
    double sum_Ipct = 0.0;
    for (unsigned subset_index = 0; subset_index < _num_subsets; ++subset_index)
        {
        // Compute I for previous clade
        _clade_prob[subset_index] = _clade_denom[subset_index]/_tree_counts[subset_index]; // marginal posterior clade probability

        _I[subset_index] = _clade_prob[subset_index]*(_clade_Hp[subset_index] - _clade_H[subset_index]);
        _total_I[subset_index] += _I[subset_index];

        _Ipct[subset_index] = 100.0*_I[subset_index]/_total_entropy;
        sum_Ipct += _Ipct[subset_index];
        }

    // handle merged case
    // Compute I for previous clade
    _clade_prob[_num_subsets] = _clade_denom[_num_subsets]/_total_trees; // marginal posterior clade probability

    _I[_num_subsets] = _clade_prob[_num_subsets]*(_clade_Hp[_num_subsets] - _clade_H[_num_subsets]);
    _total_I[_num_subsets] += _I[_num_subsets];

    _Ipct[_num_subsets] = 100.0*_I[_num_subsets]/_total_entropy;

    // Compute clade-specific D
    double D = (sum_Ipct/_num_subsets) - _Ipct[_num_subsets];

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

std::pair<unsigned,double> GalaxData::estimateCoverage(CCDMapType & ccdmap, SubsetTreeSetType & treeCCD, unsigned subset_index)
    {
    double max_log_product = 0.0;
    std::vector< double > log_products;
    const TreeIDSetType & tree_set = treeCCD[subset_index];
    unsigned num_unique_trees = (unsigned)tree_set.size();

    for (TreeIDSetType::iterator tree_iter = tree_set.begin(); tree_iter != tree_set.end(); ++tree_iter)
        {
        const TreeIDType & tree_id = *tree_iter;

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
            }
        if (log_product > max_log_product)
            max_log_product = log_product;
        log_products.push_back(log_product);
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

    return std::pair<unsigned,double>(num_unique_trees, exp_log_coverage);
    }
    
void GalaxData::calcCoverage(CCDMapType & ccdmap, SubsetTreeSetType & treeCCD)
    {
    // Estimate coverage and determine number of unique tree topologies for each subset
    for (unsigned subset_index = 0; subset_index < _num_subsets; ++subset_index)
        {
        std::pair<unsigned,double> unicov = estimateCoverage(ccdmap, treeCCD, subset_index);
        _unique[subset_index] = unicov.first;
        _coverage[subset_index] = unicov.second;
        }
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
    GalaxInfo::_sortby_index = 0;
    std::sort(subset_info.begin(), subset_info.end(), std::greater<GalaxInfo>());

    const char indiv_summary[]  = "%20s %12.5f %12.5f %12.5f %12.5f %12s\n";
    const char merged_summary[] = "%20s %12s %12s %12.5f %12.5f %12.5f\n";
    infostr += std::string("\nTotals\n");
    infostr += boost::str(boost::format("%20s %12s %12s %12s %12s %12s\n")
        % "treefile"
        % "unique"
        % "coverage"
        % "I"
        % "Ipct"
        % "D");
    unsigned subset_index = 0;
    BOOST_FOREACH(GalaxInfo & c, subset_info)
        {
        double info = c._value[0];
        double unique = c._value[1];
        double coverage = c._value[2];
        double pct = 100.0*info/_total_entropy;
        infostr += boost::str(boost::format(indiv_summary)
            % c._name
            % unique
            % coverage
            % info
            % pct
            % "---");
        ++subset_index;
        }

    // Sort clades for merged case by Ipct, D, and w and save majrule tree
    double totalIpct = 0.0;
    double cumpct = 0.0;
    if (_num_subsets > 1)
        {
        totalIpct = (100.0*_total_I[_num_subsets]/_total_entropy);
        infostr += boost::str(boost::format(merged_summary)
            % "merged"
            % "---"
            % "---"
            % _total_I[_num_subsets]
            % (100.0*_total_I[_num_subsets]/_total_entropy)
            % _total_D);

        // Report merged info for each clade sorted from highest to lowest
        GalaxInfo::_sortby_index = 0;
        std::sort(clade_info.begin(), clade_info.end(), std::greater<GalaxInfo>());
        const char clade_summary[]  = "%12.5f %12.5f %12.5f %12.5f %12.5f %s\n";
        infostr += std::string("\nClades sorted by merged info (top 95% shown):\n");
        infostr += boost::str(boost::format("%12s %12s %12s %12s %12s %s\n")
            % "I"
            % "%"
            % "cum. %"
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
        GalaxInfo::_sortby_index = 1;
        std::sort(clade_info.begin(), clade_info.end(), std::greater<GalaxInfo>());
        infostr += std::string("\nClades sorted by D (top 95% shown):\n");
        infostr += boost::str(boost::format("%12s %12s %12s %12s %12s %s\n")
            % "D"
            % "%"
            % "cum. %"
            % "Ipct"
            % "w"
            % "clade");
        cumpct = 0.0;
        BOOST_FOREACH(GalaxInfo & c, clade_info)
            {
            double info = c._value[0];
            double diff = c._value[1];
            double w    = c._value[2];
            double pct  = 100.0*diff/_total_D;
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
    GalaxInfo::_sortby_index = 2;
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

