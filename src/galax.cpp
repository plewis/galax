//
//  galax.cpp
//  galax
//
//  Created by Paul O. Lewis on 10/25/14.
//  Copyright (c) 2014 Paul O. Lewis. All rights reserved.
//

#include <fstream>
#include <boost/format.hpp>
#include "galax.hpp"
#include "galaxutil.hpp"
#include "xgalax.hpp"

using namespace galax;

Galax::Galax(const std::string outfile_prefix)
    {
    _outfprefix = outfile_prefix;
    std::string outfname = outfile_prefix + ".txt";
    _outf.open(outfname.c_str());
    }

Galax::~Galax()
    {
    _outf.close();
    }

void Galax::getTreesFromFile(std::string treefname, unsigned skip)
    {
    _start_time = getCurrentTime();

    std::string file_contents;
    getFileContents(file_contents, treefname);
    parseTranslate(_translate, file_contents);

    std::vector< std::string > tree_descriptions;
    getNewicks(tree_descriptions, file_contents, skip);
    _tree_counts.push_back((unsigned)tree_descriptions.size());
    _newicks.clear();
    _newicks.insert(_newicks.end(), tree_descriptions.begin(), tree_descriptions.end());
    _merged_newicks.insert(_merged_newicks.end(), tree_descriptions.begin(), tree_descriptions.end());

    _end_time = getCurrentTime();
    _total_seconds += secondsElapsed(_start_time, _end_time);
    }

std::vector<std::string> Galax::getTreeFileList(std::string listfname)
    {
    _start_time = getCurrentTime();

    std::string file_contents;
    getFileContents(file_contents, listfname);

    std::vector<std::string> tree_file_names;
    extractAllWhitespaceDelimitedStrings(tree_file_names, file_contents);

    _end_time = getCurrentTime();
    _total_seconds += secondsElapsed(_start_time, _end_time);

    return tree_file_names;
    }


void Galax::processTrees(TreeManip<Node>::TreeManipShPtr tm, bool rooted, unsigned subset_index, unsigned num_subsets)
    {
    _start_time = getCurrentTime();

    for (std::vector< std::string >::const_iterator sit = _newicks.begin(); sit != _newicks.end(); ++sit)
        {
        try
            {
            tm->buildFromNewick(*sit, rooted);
            tm->addToCCDMap(subset_index, num_subsets);
            }
        catch(XGalax & x)
            {
            _outf << "ERROR: " << x.what() << std::endl;
            std::exit(1);
            }
        }

    _end_time = getCurrentTime();
    _total_seconds += secondsElapsed(_start_time, _end_time);
    }

double Galax::estimateInfo(TreeManip<Node>::TreeManipShPtr tm, unsigned subset_index, unsigned num_subsets, bool details)
    {
    std::string infostr;

    _start_time = getCurrentTime();

    double info = 0.0;
    unsigned which_subset = num_subsets + 1;
    if (subset_index == TreeManip<Node>::ALLSUBSETS)
        {
        info = tm->estimateMergedInfo(infostr, _tree_counts);
        }
    else
        {
        which_subset = subset_index + 1;
        info = tm->estimateSubsetInfo(infostr, (unsigned)_newicks.size(), subset_index);
        }

    if (details)
        {
        std::string outfname = boost::str(boost::format("%s-%d.txt") % _outfprefix % which_subset);
        std::ofstream tmpf(outfname.c_str());
        tmpf << infostr << std::endl;
        tmpf.close();
        }

    _end_time = getCurrentTime();
    _total_seconds += secondsElapsed(_start_time, _end_time);

    return info;
    }

void Galax::run(std::string treefname, std::string listfname, unsigned skip, bool rooted, bool details)
    {
    std::vector<std::string> treefiles;
    if (listfname.size() > 0)
        treefiles = getTreeFileList(listfname);
    else
        treefiles.push_back(treefname);

    unsigned ntreefiles = (unsigned)treefiles.size();
    if (ntreefiles == 1)
        _outf << "Read 1 tree file name from list file " << listfname << std::endl;
    else
        _outf << "Read " << ntreefiles << " tree file names from list file " << listfname << std::endl;

    unsigned n = (unsigned)treefiles.size();
    TreeManip<Node>::TreeManipShPtr tm(new TreeManip<Node>());
    std::string infostr;
    double info;
    double sum_info = 0.0;
    unsigned total_trees = 0;

    unsigned subset_index = 0;
    _outf << boost::str(boost::format("\n%12s %12s %12s       %s\n") % "Number" % "Info" % "Trees" % "Description");
    for (std::vector<std::string>::const_iterator it = treefiles.begin(); it != treefiles.end(); ++it)
        {
        std::string treefname = *it;

        getTreesFromFile(treefname, skip);
        processTrees(tm, rooted, subset_index, n);

        unsigned ntrees = (unsigned)_newicks.size();
        total_trees += ntrees;

        //tm->showCCDMap(subset_index);
        
        info = estimateInfo(tm, subset_index, ntreefiles, details);
        sum_info += info;

        _outf << boost::str(boost::format("%12d %12.5f %12d       %s\n") % (subset_index+1) % info % ntrees % treefname.c_str());
        ++subset_index;
        }

    double average_info = sum_info/ntreefiles;
    _outf << boost::str(boost::format("%12s %12.5f %12d       %s\n") % " " % average_info % total_trees % "average");

    info = estimateInfo(tm, TreeManip<Node>::ALLSUBSETS, ntreefiles, details);
    _outf << boost::str(boost::format("%12d %12.5f %12d       %s\n") % (ntreefiles + 1) % info % total_trees % "merged");

    double diff = average_info - info;
    _outf << boost::str(boost::format("%12s %12.5f %12d       %s\n") % " " % diff % total_trees % "difference");

    _outf << "\nRequired " << _total_seconds << " total seconds" << std::endl;
    }
