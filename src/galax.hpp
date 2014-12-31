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

        void run(std::string treefname, std::string listfname, unsigned skip, bool trees_rooted, bool save_details);

    private:
        double estimateInfo(TreeManip<Node>::TreeManipShPtr tm, unsigned subset_index, unsigned num_subsets, bool details);
        void processTrees(TreeManip<Node>::TreeManipShPtr tm, bool rooted, unsigned subset_index, unsigned num_subsets);
        void getTreesFromFile(std::string treefname, unsigned skip);
        std::vector<std::string> getTreeFileList(std::string listfname);

    private:
        std::vector< std::string > _newicks;
        std::vector< std::string > _merged_newicks;
        std::map< unsigned, std::string > _translate;
        boost::posix_time::ptime _start_time;
        boost::posix_time::ptime _end_time;
        double _total_seconds;
        std::string _outfprefix;
        std::ofstream _outf;
    };

}

#endif /* defined(GALAX_GALAX_HPP) */
