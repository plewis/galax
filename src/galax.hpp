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

                                            Galax(const std::string outfile_prefix);
                                            ~Galax();

        void                                run(std::string treefname, std::string listfname, unsigned skip, bool trees_rooted);

    private:

        void                                showCCDMap(unsigned subset_index);
        void                                estimateInfo(TreeManip<Node>::TreeManipShPtr tm, std::string & infostr, std::string & majrule_newick);
        void                                writeMajruleTreefile(std::string & majrule_newick);
        void                                processTrees(TreeManip<Node>::TreeManipShPtr tm, unsigned subset_index, unsigned num_subsets);
        void                                getTreesFromFile(std::string treefname, unsigned skip);
        std::vector<std::string>            getTreeFileList(std::string listfname);

    private:
        static const unsigned               ALLSUBSETS;
        std::vector< std::string >          _treefile_names;
        std::vector< std::string >          _newicks;
        std::vector< std::string >          _merged_newicks;
        std::vector< unsigned >             _tree_counts;
        std::map< unsigned, std::string >   _translate;
        boost::posix_time::ptime            _start_time;
        boost::posix_time::ptime            _end_time;
        double                              _total_seconds;
        bool                                _rooted;
        std::string                         _outfprefix;
        std::ofstream                       _treef;
        std::ofstream                       _outf;
        CCDMapType                          _ccdmap;
    };

}

#endif /* defined(GALAX_GALAX_HPP) */
