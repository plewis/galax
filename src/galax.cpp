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

Galax::Galax(const std::string outfname)
    {
    outf.open(outfname.c_str());
    }

Galax::~Galax()
    {
    outf.close();
    }

void Galax::getTreesFromFile(std::string treefname, unsigned skip)
    {
    _start_time = getCurrentTime();

    std::string file_contents;
    getFileContents(file_contents, treefname);
    parseTranslate(_translate, file_contents);

    std::vector< std::string > tree_descriptions;
    getNewicks(tree_descriptions, file_contents, skip);
    _newicks.insert(_newicks.end(), tree_descriptions.begin(), tree_descriptions.end());

    _end_time = getCurrentTime();
    double secs = secondsElapsed(_start_time, _end_time);
    outf << "Read " << tree_descriptions.size() << " trees from file " << treefname << " in " << secs << " seconds" << std::endl;
    _total_seconds = secs;
    }

std::vector<std::string> Galax::getTreeFileList(std::string listfname)
    {
    _start_time = getCurrentTime();

    std::string file_contents;
    getFileContents(file_contents, listfname);

    std::vector<std::string> tree_file_names;
    extractAllWhitespaceDelimitedStrings(tree_file_names, file_contents);

    _end_time = getCurrentTime();
    double secs = secondsElapsed(_start_time, _end_time);
    unsigned ntreefilenames = (unsigned)tree_file_names.size();
    if (ntreefilenames == 1)
        outf << "Read 1 tree file name from list file " << listfname << " in " << secs << " seconds" << std::endl;
    else
        outf << "Read " << ntreefilenames << " tree file names from list file " << listfname << " in " << secs << " seconds" << std::endl;
    _total_seconds = secs;

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
            outf << "ERROR: " << x.what() << std::endl;
            std::exit(1);
            }
        }

    _end_time = getCurrentTime();
    double secs = secondsElapsed(_start_time, _end_time);
    outf << "Processed " << _newicks.size() << " trees in " << secs << " seconds" << std::endl;
    _total_seconds += secs;

    //tm->showCCDMap(subset_index);
    
    _start_time = getCurrentTime();
    std::string infostr = tm->estimateInfo((unsigned)_newicks.size(), subset_index);
    _end_time = getCurrentTime();

    secs = secondsElapsed(_start_time, _end_time);
    _total_seconds += secs;

    outf << infostr << std::endl;
    outf << "Required " << secs << " seconds to compute information content" << std::endl;
    outf << "Required " << _total_seconds << " total seconds" << std::endl;
    }

void Galax::run(std::string treefname, std::string listfname, unsigned skip, bool rooted)
    {
    std::vector<std::string> treefiles;
    if (listfname.size() > 0)
        treefiles = getTreeFileList(listfname);
    else
        treefiles.push_back(treefname);

    unsigned i = 0;
    unsigned n = (unsigned)treefiles.size();
    TreeManip<Node>::TreeManipShPtr tm(new TreeManip<Node>());
    for (std::vector<std::string>::const_iterator it = treefiles.begin(); it != treefiles.end(); ++it)
        {
        _newicks.clear();
        getTreesFromFile(*it, skip);
        processTrees(tm, rooted, i++, n);
        }
    }
