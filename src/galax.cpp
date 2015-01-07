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

void Galax::estimateInfo(TreeManip<Node>::TreeManipShPtr tm, std::string & infostr, std::string & majrule_newick)
    {
    _start_time = getCurrentTime();

    tm->estimateMergedInfo(infostr, majrule_newick, _tree_counts, _treefile_names);

    _end_time = getCurrentTime();
    _total_seconds += secondsElapsed(_start_time, _end_time);
    }

void Galax::writeMajruleTreefile(std::string & majrule_newick)
    {
    _start_time = getCurrentTime();

    std::string treefname = _outfprefix + ".tre";
    _treef.open(treefname.c_str());

    _treef << "#nexus\n\nbegin trees;\n  translate\n";
    typedef std::map<unsigned, std::string> translate_map_type;
    translate_map_type::iterator it = _translate.begin();
    _treef << "  " << it->first << " '" << it->second << "'";
    for (; it != _translate.end(); ++it)
        {
        _treef << ",\n  " << it->first << " '" << it->second << "'";
        }
    _treef << ";\ntree majrule = " << majrule_newick << ";\nend;" << std::endl;
    _treef.close();

    _end_time = getCurrentTime();
    _total_seconds += secondsElapsed(_start_time, _end_time);
    }

void Galax::run(std::string treefname, std::string listfname, unsigned skip, bool rooted)
    {
	try
		{
        _treefile_names.clear();
        if (listfname.size() > 0)
            _treefile_names = getTreeFileList(listfname);
        else
            _treefile_names.push_back(treefname);

        unsigned ntreefiles = (unsigned)_treefile_names.size();
        if (ntreefiles == 1)
            _outf << "Read 1 tree file name from list file " << listfname << std::endl;
        else
            _outf << "Read " << ntreefiles << " tree file names from list file " << listfname << std::endl;

        unsigned n = (unsigned)_treefile_names.size();
        TreeManip<Node>::TreeManipShPtr tm(new TreeManip<Node>());

        unsigned subset_index = 0;
        for (std::vector<std::string>::const_iterator it = _treefile_names.begin(); it != _treefile_names.end(); ++it)
            {
            std::string treefname = *it;

            getTreesFromFile(treefname, skip);
            processTrees(tm, rooted, subset_index, n);

            ++subset_index;
            }

        std::string infostr;
        std::string majrule_newick;
        estimateInfo(tm, infostr, majrule_newick);
        writeMajruleTreefile(majrule_newick);
        _outf << "\n" << infostr;

        _outf << "\nRequired " << _total_seconds << " total seconds" << std::endl;
        }
	catch(XGalax x)
		{
		std::cerr << "Exception: " << x.what() << std::endl;
		}
    }
