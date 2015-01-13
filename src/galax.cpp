//
//  galax.cpp
//  galax
//
//  Created by Paul O. Lewis on 10/25/14.
//  Copyright (c) 2014 Paul O. Lewis. All rights reserved.
//

#include <fstream>
#include <limits>
#include <boost/format.hpp>
#include "galax.hpp"
#include "galaxutil.hpp"
#include "xgalax.hpp"

using namespace galax;

const unsigned Galax::_ALLSUBSETS = std::numeric_limits<unsigned>::max();

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
    bool complete = parseTranslate(file_contents);
    if (!complete)
        throw XGalax(boost::str(boost::format("File %s contains a different number of taxa than other files. Galax requires taxa to be the same for all tree files processed.") % treefname));

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

bool Galax::parseTranslate(const std::string & file_contents)
    {
    unsigned prev_ntaxa = (unsigned)_translate.size();

    // Will either add to an empty _translate or check values if _translate has been created previously

    // This set will contain the taxon index of all taxa encountered in translate statement. It should have the
    // same number of elements as _translate has keys when done, otherwise this file lacks one or more taxa
    std::vector<unsigned> taxa_seen;

    // First, separate out the contents of the translate statement from the rest of the tree file
    std::string translate_contents;
    boost::regex pattern("[Tt]ranslate(.+?);");
    boost::smatch what;
    bool regex_ok = boost::regex_search(file_contents, what, pattern);
    if (regex_ok)
        {
        // what[0] contains the whole string
        // what[1] contains the translate statement contents
        // Construct a string using characters in contents from what[1].first to what[1].second
        translate_contents.insert(translate_contents.begin(), what[1].first, what[1].second);
        }
    else
        {
        throw XGalax("regex failed to find translate statement in tree file");
        }

    // Now create the map by iterating through each element of the translate statement
    boost::regex re("(\\d+)\\s+'?(.+?)'?,?$");
    boost::sregex_iterator m1(translate_contents.begin(), translate_contents.end(), re);
    boost::sregex_iterator m2;
    for (; m1 != m2; ++m1)
        {
        const boost::match_results<std::string::const_iterator>& what = *m1;
        unsigned taxon_index = 0;
        try
            {
            taxon_index = boost::lexical_cast<unsigned>(what[1].str());
            }
        catch(const boost::bad_lexical_cast &)
            {
            throw XGalax("Could not interpret taxon index in translate statement as a number");
            }

        taxa_seen.push_back(taxon_index);

        // efficientAddOrCheck returns valid iterator on first insertion or if identical association already previously made
        if (efficientAddOrCheck(_translate, taxon_index, what[2].str()) == _translate.end())
            throw XGalax(boost::str(boost::format("Taxon name (%s) does not match name already associated (%s) with taxon %d") % what[2].str() % _translate[taxon_index] % taxon_index));
        }

    // Check to make sure this file does not have one or more additional taxa than files previously processed
    if (prev_ntaxa > 0 && prev_ntaxa < _translate.size())
        return false;

    // Check to make sure this file does not have one or more fewer taxa than files previously processed
    if (taxa_seen.size() < _translate.size())
       return false;

    return true;
    }

void Galax::getNewicks(std::vector< std::string > & tree_descriptions, const std::string & file_contents, unsigned skip)
    {
    // tree STATE_0 [&lnP=-4493.80476846934,posterior=-4493.80476846934] = [&R] (9:[&rate=0.49971158909783764]1851.4724462198697,((1:[&rate=0.5965730394621352]292.73199858783727,(10:[&rate=0.6588031360335018]30.21172743645451,5:[&rate=0.7098036299867017]30.21172743645451):[&rate=1.0146941544458208]262.52027115138276):[&rate=1.0649642758561977]510.452441519872,((8:[&rate=0.7554924641162211]145.1076605992074,(7:[&rate=0.7984750329147966]64.0435017480143,6:[&rate=0.8402528958963882]64.0435017480143):[&rate=1.1206854064651213]81.06415885119311):[&rate=1.1844450597679457]522.823827314411,((3:[&rate=0.8818808237868384]60.2962343089954,4:[&rate=0.9242396697890951]60.2962343089954):[&rate=1.260685743226102]12.793802911399744,2:[&rate=0.9681896872253556]73.09003722039515):[&rate=1.3582802932633053]594.8414506932232):[&rate=1.4999660689010508]135.25295219409088):[&rate=1.7907115550989796]1048.2880061121605);
    tree_descriptions.clear();
    boost::regex re("^\\s*[Tt]ree.+?(\\(.+?\\))\\s*;\\s*$");
    boost::sregex_iterator m1(file_contents.begin(), file_contents.end(), re);
    boost::sregex_iterator m2;  // empty iterator used only to detect when we are done
    unsigned n = 0;
    for (; m1 != m2; ++m1)
        {
        const boost::match_results<std::string::const_iterator>& what = *m1;
        if (n >= skip)
            {
            tree_descriptions.push_back( stripComments( what[1].str() ) );
            }
        n += 1;
        }
    }

void Galax::processTrees(TreeManip<Node>::TreeManipShPtr tm, CCDMapType & ccdmap, unsigned subset_index, unsigned num_subsets)
    {
    _start_time = getCurrentTime();

    for (std::vector< std::string >::const_iterator sit = _newicks.begin(); sit != _newicks.end(); ++sit)
        {
        try
            {
            tm->buildFromNewick(*sit, (_rooted ? 0 : _outgroup));
            tm->addToCCDMap(ccdmap, subset_index, num_subsets);
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

void Galax::showCCDMap(CCDMapType & ccdmap, unsigned subset_index)
	{
    std::cerr << "CCD map has " << ccdmap.size() << " elements." << std::endl;
    for (CCDMapType::iterator it = ccdmap.begin(); it != ccdmap.end(); ++it)
        {
        const SplitVector & v = it->first;
        double count = 0.0;
        if (subset_index < _ALLSUBSETS)
            count = it->second[subset_index];
        else
            count = std::accumulate(it->second.begin(), it->second.end(), 0.0);
        std::cerr << "\n+--------------------------" << std::endl;
        for (SplitVector::const_iterator vit = v.begin(); vit != v.end(); ++vit)
            {
            _outf << "| " << vit->createPatternRepresentation() << std::endl;
            }
        _outf << "| " << count << std::endl;
        }
    }

void Galax::estimateInfo(TreeManip<Node>::TreeManipShPtr tm, CCDMapType & ccdmap, std::string & infostr, std::vector<GalaxInfo> & clade_info)
    {
    _start_time = getCurrentTime();

    assert(ccdmap.size() > 0);

    unsigned num_subsets = (unsigned)_tree_counts.size();
    unsigned total_trees = std::accumulate(_tree_counts.begin(), _tree_counts.end(), 0);

    // Determine maximum possible entropy
    Split first_clade = (ccdmap.begin()->first)[0];
    unsigned ntaxa = first_clade.countOnBits() + first_clade.countOffBits();
    if (!_rooted)
        ntaxa -= 1;
    double total_entropy = lognrooted(ntaxa);

    infostr = "Order of taxa in split representations:\n";
    typedef std::pair<unsigned, std::string> translate_map_type;
    BOOST_FOREACH(const translate_map_type & translate_key_value, _translate)
        {
        infostr += boost::str(boost::format("%12d  %s\n") % translate_key_value.first % translate_key_value.second);
        }

    infostr += "\nLindley Information\n";
    infostr += boost::str(boost::format("  Number of taxa in ingroup: %d\n") % ntaxa);
    infostr += boost::str(boost::format("  Total prior entropy: %.5f\n") % total_entropy);

    Split clade;
    clade_info.clear();

    // These vectors store information for one particular clade for each subset plus the merged case
    // All of these vectors are reinitialized for every new clade: not currently storing all information
    // for every subset-clade combination
    std::vector<double> clade_H(num_subsets+1, 0.0);
    std::vector<double> clade_Hp(num_subsets+1, 0.0);
    std::vector<double> clade_denom(num_subsets+1, 0.0);
    std::vector<double> w(num_subsets+1, 0.0);
    std::vector<double> I(num_subsets+1, 0.0);
    std::vector<double> Ipct(num_subsets+1, 0.0);

    // This vector holds cumulative I across clades for each subset plus the merged case
    // It is not reset for each clade
    std::vector<double> total_I(num_subsets+1, 0.0);

    double total_D = 0.0;
    bool first = true;
    const char indiv_template[]  = "%20s %12d %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12s\n";
    const char merged_template[] = "%20s %12d %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f\n";
    for (CCDMapType::iterator it = ccdmap.begin(); it != ccdmap.end(); ++it)
        {
        // v contains either 1 split (unconditional) or 3 splits (conditional)
        const SplitVector & v = it->first;

        std::vector<double> & count = it->second;
        double sum_counts = std::accumulate(count.begin(), count.end(), 0.0);

        // Depends on ccdmap keys (split vectors) being sorted such that an entry for an unconditional clade (e.g. ABC)
        // precedes any conditional clade entries (e.g. A|BC) that have that clade as parent. Thus, when next unconditional
        // clade entry is encountered, we know that there are no more conditional clade entries for the previous parent
        // clade and we can at that point compute the information for the previous parent.
        if (first)
            {
            first = false;

            // Initialize data for first clade
            assert(v.size() == 1);  // first entry should be an unconditional clade entry
            clade = v[0];

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
                    w[subset_index] = clade_denom[subset_index]/_tree_counts[subset_index]; // marginal posterior clade probability

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

                if (I[num_subsets] > 0.0)   // don't bother outputting numbers for clades with no information
                    {
                    infostr += boost::str(boost::format("\nClade %s\n") % clade.createPatternRepresentation());
                    infostr += boost::str(boost::format("%20s %12s %12s %12s %12s %12s %12s %12s %12s\n")
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
                        assert(_treefile_names.size() > 0);
                        assert(_tree_counts.size() > 0);
                        assert(clade_Hp.size() > 0);
                        assert(clade_H.size() > 0);
                        assert(w.size() > 0);
                        assert(I.size() > 0);
                        assert(Ipct.size() > 0);
                        infostr += boost::str(boost::format(indiv_template)
                            % _treefile_names[subset_index]
                            % _tree_counts[subset_index]
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
                        infostr += boost::str(boost::format(merged_template)
                            % "merged"
                            % total_trees
                            % clade_Hp[num_subsets]
                            % clade_H[num_subsets]
                            % (clade_Hp[num_subsets] - clade_H[num_subsets])
                            % w[num_subsets]
                            % I[num_subsets]
                            % Ipct[num_subsets]
                            % D);
                        }
                    std::vector<double> tmp;
                    tmp.push_back(Ipct[num_subsets]);
                    tmp.push_back(D);
                    tmp.push_back(w[num_subsets]);
                    clade_info.push_back(GalaxInfo(clade.createPatternRepresentation(), tmp));
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
            w[subset_index] = clade_denom[subset_index]/_tree_counts[subset_index];

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
        infostr += boost::str(boost::format("\nClade %s\n") % clade.createPatternRepresentation());
        infostr += boost::str(boost::format("%20s %12s %12s %12s %12s %12s %12s %12s %12s\n")
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
            infostr += boost::str(boost::format(indiv_template)
                % _treefile_names[subset_index]
                % _tree_counts[subset_index]
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
            infostr += boost::str(boost::format(merged_template)
                % "merged"
                % total_trees
                % clade_Hp[num_subsets]
                % clade_H[num_subsets]
                % (clade_Hp[num_subsets] - clade_H[num_subsets])
                % w[num_subsets]
                % I[num_subsets]
                % Ipct[num_subsets]
                % D);
            }
            std::vector<double> tmp;
            tmp.push_back(Ipct[num_subsets]);
            tmp.push_back(D);
            tmp.push_back(w[num_subsets]);
            clade_info.push_back(GalaxInfo(clade.createPatternRepresentation(), tmp));
        }

    // Report totals for each subset
    std::vector<GalaxInfo> subset_info;
    for (unsigned subset_index = 0; subset_index < num_subsets; ++subset_index)
        {
        std::vector<double> tmp;
        tmp.push_back(total_I[subset_index]);
        subset_info.push_back(GalaxInfo(_treefile_names[subset_index], tmp));
        }
    GalaxInfo::_sortby_index = 0;
    std::sort(subset_info.begin(), subset_info.end(), std::greater<GalaxInfo>());

    const char indiv_summary[]  = "%20s %12.5f %12.5f %12s\n";
    const char merged_summary[] = "%20s %12.5f %12.5f %12.5f\n";
    infostr += std::string("\nTotals\n");
    infostr += boost::str(boost::format("%20s %12s %12s %12s\n")
        % "treefile"
        % "I"
        % "Ipct"
        % "D");
    BOOST_FOREACH(GalaxInfo & c, subset_info)
        {
        double info = c._value[0];
        double pct = 100.0*info/total_entropy;
        infostr += boost::str(boost::format(indiv_summary)
            % c._name
            % info
            % pct
            % "---");
        }

    // Sort clades for merged case by Ipct, D, and w and save majrule tree
    double totalIpct = 0.0;
    double cumpct = 0.0;
    if (num_subsets > 1)
        {
        totalIpct = (100.0*total_I[num_subsets]/total_entropy);
        infostr += boost::str(boost::format(merged_summary)
            % "merged"
            % total_I[num_subsets]
            % (100.0*total_I[num_subsets]/total_entropy)
            % total_D);

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
            double pct  = 100.0*diff/total_D;
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
    double log2 = log(2.0);
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

    infostr += "\nOrder of taxa in split representations:\n";
    typedef std::pair<unsigned, std::string> translate_map_type;
    BOOST_FOREACH(const translate_map_type & translate_key_value, _translate)
        {
        infostr += boost::str(boost::format("%12d  %s\n") % translate_key_value.first % translate_key_value.second);
        }

    _end_time = getCurrentTime();
    _total_seconds += secondsElapsed(_start_time, _end_time);
    }

void Galax::buildMajorityRuleTree(std::vector<GalaxInfo> & majrule_info, std::vector<GalaxInfo> & annotate_info, std::string & majrule_newick)
    {
    _start_time = getCurrentTime();

    majrule_newick.clear();

    // We can save a lot of computation if both info vectors are the same
    bool same_info = (&majrule_info == &annotate_info);

    // Ensure that both majrule_info and annotate_info are sorted so that the clades with highest posterior are first
    GalaxInfo::_sortby_index = 2;
    std::sort(majrule_info.begin(), majrule_info.end(), std::greater<GalaxInfo>());
    if (!same_info)
        std::sort(annotate_info.begin(), annotate_info.end(), std::greater<GalaxInfo>());

    // Go through all clades in majrule_info and store splits with w > 0.5
    double log2 = log(2.0);
    std::vector<Split> majrule_splits;
    BOOST_FOREACH(GalaxInfo & c, majrule_info)
        {
        double w    = c._value[2];

        Split split;
        split.createFromPattern(c._name);

        if (w < 0.5)
            {
            // majrule_info is sorted by w, so we can stop once we pass the 0.5 point
            break;
            }
        else
            {
            // Add split to vector used to construct majority-rule tree
            majrule_splits.push_back(split);
            }
        }

    // Have now found all the splits comprising the majority rule tree but still
    // need to annotate each of these splits.
    if (majrule_splits.size() == 0)
        return;

    // identify most probable conflicting splits (from annotate_info) for those in majority rule tree
    BOOST_FOREACH(Split & split, majrule_splits)
        {
        bool split_found = false;
        bool ic_computed = false;
        BOOST_FOREACH(GalaxInfo & c, annotate_info)
            {
            Split asplit;
            asplit.createFromPattern(c._name);
            asplit.setWeight(c._value[2]);
            if (split == asplit)
                {
                double info = c._value[0];
                double d    = c._value[1];
                double w    = c._value[2];

                split.setInfo(info);
                split.setWeight(w);
                split.setDisparity(d);

                split_found = true;
                }
            else if (split_found && !asplit.isCompatible(split))
                {
                // calculate internode certainty
                double p = split.getWeight();
                double q = asplit.getWeight();
                double ic = 1.0 + (p/(p+q))*log(p/(p+q))/log2 + (q/(p+q))*log(q/(p+q))/log2;
                split.setCertainty(ic);
                ic_computed = true;
                break;
                }
            }
        assert(split_found);
        if (!ic_computed)
            {
            split.setCertainty(1.0);    // no conflicting splits found
            }

        //temporary!
        //std::cerr << "  " << split.createPatternRepresentation() << " w = " << split.getWeight() << ", d = " << split.getDisparity() << ", ic = " << split.getCertainty() << std::endl;
        }

    //temporary!
    //std::cerr << std::endl;

    TreeManip<Node>::TreeManipShPtr tm(new TreeManip<Node>());
    tm->buildFromSplitVector(majrule_splits, (_rooted ? 0 : _outgroup));
    majrule_newick = tm->makeNewick(5, true);

    _end_time = getCurrentTime();
    _total_seconds += secondsElapsed(_start_time, _end_time);
    }

void Galax::writeMajruleTreefile(std::string fnprefix, std::string & majrule_newick)
    {
    _start_time = getCurrentTime();

    if (majrule_newick.length() == 0)
        {
        _outf << "Majority rule tree is the star tree: no tree file written." << std::endl;
        return;
        }

    std::string treefname = boost::str(boost::format("%s-%s.tre") % _outfprefix % fnprefix);
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

void Galax::run(std::string treefname, std::string listfname, unsigned skip, bool rooted, unsigned outgroup_taxon)
    {
    assert (_ALLSUBSETS == (unsigned)(-1));
	try
		{
        _outgroup = outgroup_taxon;
        _rooted = rooted;
        if (!_rooted && _outgroup == 0)
            throw XGalax("outgroup taxon specified must be a number greater than zero unless trees are rooted (in which case outgroup specification is ignored)");

        _treefile_names.clear();
        bool is_treefile = (treefname.length() > 0);
        bool is_listfile = (listfname.length() > 0);
        if (is_listfile)
            {
            _treefile_names = getTreeFileList(listfname);
            if (_treefile_names.size() == 0)
                throw XGalax("found no tree file names in specified listfile");
            else if (_treefile_names.size() == 1)
                throw XGalax("use --treefile (not --listfile) if there is only one tree file to process");
            }
        else
            _treefile_names.push_back(treefname);

        std::string infostr;
        TreeManip<Node>::TreeManipShPtr tm(new TreeManip<Node>());
        if (is_treefile && is_listfile)
            {
            //
            // Both --listfile and --treefile provided on command line
            //
            _outf << "Read " << _treefile_names.size() << " tree file names from list file " << listfname << "\n";
            _outf << "Read 1 tree file name from tree file " << treefname << "\n";

            // process --treefile and construct majrule consensus tree
            getTreesFromFile(treefname, skip);
            processTrees(tm, _ccdtree, 0, 1);
            std::vector<GalaxInfo> majrule_clade_info;
            std::vector<Split> majrule_splits;
            std::string majrule_tree;
            estimateInfo(tm, _ccdtree, infostr, majrule_clade_info);
            buildMajorityRuleTree(majrule_clade_info, majrule_clade_info, majrule_tree);
            writeMajruleTreefile("majrule", majrule_tree);

            // process --listfile and construct majrule consensus tree for merged trees
            // note: infostr is overwritten
            unsigned subset_index = 0;
            _tree_counts.clear();
            BOOST_FOREACH(std::string & tree_file_name, _treefile_names)
                {
                getTreesFromFile(tree_file_name, skip);
                processTrees(tm, _ccdlist, subset_index++, (unsigned)_treefile_names.size());
                }
            std::vector<GalaxInfo> merged_clade_info;
            std::vector<Split> merged_splits;
            std::string merged_tree;
            estimateInfo(tm, _ccdlist, infostr, merged_clade_info);
            buildMajorityRuleTree(merged_clade_info, merged_clade_info, merged_tree);
            writeMajruleTreefile("merged", merged_tree);

            // finally, add information content from --listfile tree files onto majority rule tree derived from --treefile
            std::string majrule_merged_tree;
            buildMajorityRuleTree(majrule_clade_info, merged_clade_info, majrule_merged_tree);
            writeMajruleTreefile("majrule-merged", majrule_merged_tree);
            }
        else if (is_treefile)
            {
            //
            // Only --treefile provided on command line
            //
            _outf << "Read 1 tree file name from tree file " << treefname << "\n";

            getTreesFromFile(treefname, skip);
            processTrees(tm, _ccdtree, 0, 1);

            std::vector<GalaxInfo> majrule_clade_info;
            std::vector<Split> majrule_splits;
            std::string majrule_tree;
            estimateInfo(tm, _ccdtree, infostr, majrule_clade_info);
            buildMajorityRuleTree(majrule_clade_info, majrule_clade_info, majrule_tree);
            writeMajruleTreefile("majrule", majrule_tree);
            }
        else
            {
            //
            // Only --listfile provided on command line
            //
            _outf << "Read " << _treefile_names.size() << " tree file names from list file " << listfname << "\n";

            unsigned subset_index = 0;
            BOOST_FOREACH(std::string & tree_file_name, _treefile_names)
                {
                getTreesFromFile(tree_file_name, skip);
                processTrees(tm, _ccdlist, subset_index++, (unsigned)_treefile_names.size());
                }
            std::vector<GalaxInfo> merged_clade_info;
            std::vector<Split> merged_splits;
            std::string merged_tree;
            estimateInfo(tm, _ccdlist, infostr, merged_clade_info);
            buildMajorityRuleTree(merged_clade_info, merged_clade_info, merged_tree);
            writeMajruleTreefile("merged", merged_tree);
            }

        if (_rooted)
            _outf << "Input trees assumed to be rooted\n";
        else
            {
            _outf << "Input trees assumed to be unrooted\n";
            _outf << boost::str(boost::format("Each input tree was rooted at outgroup taxon %d (\"%s\")\n") % _outgroup % _translate[_outgroup]);
            }

        _outf << "\n" << infostr;

        _outf << "\nRequired " << _total_seconds << " total seconds" << std::endl;
        }
	catch(XGalax x)
		{
		std::cerr << "ERROR: " << x.what() << std::endl;
		}
    }

