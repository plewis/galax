//
//  galax.hpp
//  galax
//
//  Created by Paul O. Lewis on 10/25/14.
//  Copyright (c) 2014 Paul O. Lewis. All rights reserved.
//

#ifndef GALAX_GALAX_HPP
#define GALAX_GALAX_HPP

#include <string>
#include <vector>
#include <map>
#include <boost/date_time/posix_time/posix_time.hpp>
#include "treemanip.hpp"

namespace galax {

extern const double _smallest_edge_length;

class Galax
    {
    public:

                                            Galax(const std::string outfile_prefix, const std::string version);
                                            ~Galax();

        void                                run(std::string treefname, std::string listfname, unsigned skip, bool trees_rooted, std::string maptofname, bool mapto_trees_rooted, bool save_details, unsigned outgroup_taxon);

    private:

        void                                debugShowCCDMap(CCDMapType & ccdmap, unsigned subset_index);
        void                                saveDetailedInfoForClade(std::string & detailedinfostr, std::string pattern, unsigned num_subsets, unsigned total_trees, std::vector<double> & clade_Hp, std::vector<double> & clade_H, std::vector<double> & w, std::vector<double> & I, std::vector<double> & Ipct, double D);
        void                                estimateInfo(CCDMapType & ccdmap, std::string & summaryinfostr, std::string & detailedinfostr, std::vector<GalaxInfo> & clade_info);
        void                                buildMajorityRuleTree(std::vector<GalaxInfo> & majrule_info, std::vector<GalaxInfo> & annotate_info, std::string & majrule_newick);
        void                                writeMajruleTreefile(std::string fnprefix, std::string & majrule_newick);
        void                                mapToTree(std::string & maptofname, std::vector<GalaxInfo> & annotate_info, std::string & mapto_newick);
        void                                initTreeCCD(unsigned num_subsets);
        void                                processTrees(TreeManip<Node>::TreeManipShPtr tm, CCDMapType & ccdmap, unsigned subset_index, unsigned num_subsets);
        void                                storeTrees(std::string file_contents, unsigned skip, std::vector< std::string > & tree_descriptions);
        //void                              getTrees(std::string file_contents, unsigned skip);
        bool                                isNexusFile(const std::string & file_contents);
        unsigned                            taxonNumberFromName(const std::string taxon_name, bool add_if_missing);
        bool                                replaceTaxonNames(const std::string & newick_with_taxon_names, std::string & newick_with_taxon_numbers);
        void                                parseTranslate(const std::string & file_contents);
        void                                getPhyloBayesNewicks(std::vector< std::string > & tree_descriptions, const std::string & file_contents, unsigned skip);
        std::string                         standardizeNodeNumber(boost::smatch const & what);
        std::string                         standardizeTreeDescription(std::string & newick);
        void                                getNewicks(std::vector< std::string > & tree_descriptions, const std::string & file_contents, unsigned skip);
        std::vector<std::string>            getTreeFileList(std::string listfname);

        // profile
        void                                writeInfoProfile(TreeManip<Node>::TreeManipShPtr tm, std::vector<GalaxInfo> & clade_info);

    private:
        static const unsigned               _ALLSUBSETS;
        std::vector< std::string >          _treefile_names;
        std::vector< std::string >          _newicks;
        std::vector< std::string >          _merged_newicks;
        std::vector< unsigned >             _tree_counts;
        std::map< unsigned, std::string >   _translate;
        std::map< std::string, unsigned >   _taxon_map;
        boost::posix_time::ptime            _start_time;
        boost::posix_time::ptime            _end_time;
        double                              _total_seconds;
        bool                                _rooted;
        bool                                _input_rooted;
        bool                                _mapto_rooted;
        std::string                         _outfprefix;
        std::string                         _version;
        bool                                _show_details;

        std::ofstream                       _treef;     // file for storing majority rule consensus tree
        std::ofstream                       _outf;      // file used for output summary
        std::ofstream                       _detailsf;  // file used for details of clades

        unsigned                            _outgroup;  // 1-based index of outgroup taxon to use for rooting unrooted trees (1 assumed if not specified)

        CCDMapType                          _ccdlist;   // _ccdlist[i][j] stores conditional clade info for clade i, subset j, for trees specified by --listfile
        CCDMapType                          _ccdtree;   // _ccdtree[i][0] stores conditional clade info for clade i for trees specified by --treefile

        SubsetTreeSetType                   _treeCCD;  // _treeCCD[i] is a set of tree IDs (list of cond. clades) appearing in the sample for subset i
        SubsetTreeMapType                   _treeMap;  // _treeMap[i] is a map relating tree IDs (keys) to counts (values) of each tree topology sampled in subset i (only used if _show_details is true)
    };

}

#endif /* defined(GALAX_GALAX_HPP) */
