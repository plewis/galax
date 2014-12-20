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
#include "treemanip.hpp"
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
    getNewicks(_newicks, file_contents, skip);

    _end_time = getCurrentTime();
    double secs = secondsElapsed(_start_time, _end_time);
    outf << "Read " << _newicks.size() << " trees from file " << treefname << " in " << secs << " seconds" << std::endl;
    _total_seconds = secs;
    }

void Galax::processTrees(bool rooted)
    {
    _start_time = getCurrentTime();

    TreeManip<Node>::TreeManipShPtr tm(new TreeManip<Node>());
    for (std::vector< std::string >::const_iterator sit = _newicks.begin(); sit != _newicks.end(); ++sit)
        {
        try
            {
            tm->buildFromNewick(*sit, rooted);
            tm->addToCCDMap();
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

    //tm->showCCDMap();
    
    _start_time = getCurrentTime();
    std::string infostr = tm->estimateInfo((unsigned)_newicks.size());
    _end_time = getCurrentTime();

    secs = secondsElapsed(_start_time, _end_time);
    _total_seconds += secs;

    outf << infostr << std::endl;
    outf << "Required " << secs << " seconds to compute information content" << std::endl;
    outf << "Required " << _total_seconds << " total seconds" << std::endl;
    }

void Galax::run(std::string treefname, unsigned skip, bool rooted)
    {
    getTreesFromFile(treefname, skip);
    processTrees(rooted);
    }
