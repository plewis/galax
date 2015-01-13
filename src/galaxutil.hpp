//
//  galaxutil.hpp
//  galax
//
//  Created by Paul O. Lewis on 10/25/14.
//  Copyright (c) 2014 Paul O. Lewis. All rights reserved.
//

#ifndef GALAX_GALAXUTIL_HPP
#define GALAX_GALAXUTIL_HPP

#include <cmath>
#include <string>
#include <set>
#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include "xgalax.hpp"

namespace galax
{

// From http://stackoverflow.com/questions/6097927/is-there-a-way-to-implement-analog-of-pythons-separator-join-in-c
template <typename Iter>
std::string join(Iter curr, Iter stop, std::string const & separator)
    {
    std::ostringstream result;
    if (curr != stop)
        result << *curr++;
    while (curr != stop)
        result << separator << *curr++;
    return result.str();
    }

inline double logfactorial(unsigned n)
    {
    // Returns log(n!)
    if (n == 0 || n == 1)
        return 0.0;
    return lgamma(n+1.0);
    }

inline double choose(unsigned n, unsigned y)
    {
    // Returns n choose y
    double tmp = logfactorial(n) - logfactorial(y) - logfactorial(n-y);
    if (n-y == y)
        tmp -= log(2.0);
    return exp(tmp);
    }

inline double lognrooted(unsigned y)
    {
    // Returns number of rooted trees of y taxa
    double fy = (double)y;
    double log_num_rooted_trees = 0.0;
    if (y > 2)
        log_num_rooted_trees = logfactorial(2*y - 3) - (fy - 2.0)*log(2.0) - logfactorial(y - 2);
    return log_num_rooted_trees;
    }

// See item 24, p. 110, in Meyers, Scott. 2001. Effective STL: 50 specific ways to improve your
// use of the Standard Template Library. Addison-Wesley, Boston.
template< typename MapType, typename KeyArgType, typename ValueArgType >
typename MapType::iterator efficientAddOrCheck(MapType & m, const KeyArgType & k, const ValueArgType & v)
    {
    // Find where k is, or should be
    typename MapType::iterator lb = m.lower_bound(k);

    // If lb points to a pair whose key is equivalent to k,
    if (lb != m.end() && !(m.key_comp()(k, lb->first)))
        {
        // check the pair's value and return m.end() if v is not equal to currently-store value or lb if they are the same
        if (lb->second == v)
            return lb;
        else
            return m.end();
        }
    else
        {
        // else add pair(k,v) to m and return an iterator to the new map element.
        typedef typename MapType::value_type MVT;
        return m.insert(lb, MVT(k, v));
        }
    }

// Modified version of item 24, p. 110, in Meyers, Scott. 2001. Effective STL: 50 specific ways to improve
// your use of the Standard Template Library. Addison-Wesley, Boston.
template< typename MapType, typename KeyArgType >
typename MapType::iterator efficientIncrement(MapType & m, const KeyArgType & k, unsigned subset_index, unsigned num_subsets)
    {
    // Find where k is, or should be
    typename MapType::iterator lb = m.lower_bound(k);

    // If lb points to a pair whose key is equivalent to k,
    if (lb != m.end() && !(m.key_comp()(k, lb->first)))
        {
        // update the pair's value and return an iterator to that pair,
        lb->second[subset_index] += 1.0;
        return lb;
        }
    else
        {
        // else add pair(k,v) to m and return an iterator to the new map element.
        typedef typename MapType::value_type MVT;
        std::vector<double> v(num_subsets, 0.0);
        v[subset_index] = 1.0;
        return m.insert(lb, MVT(k, v));
        }
    }

inline boost::posix_time::ptime getCurrentTime()
    {
    return boost::posix_time::microsec_clock::local_time();
    }

inline double secondsElapsed(boost::posix_time::ptime a, boost::posix_time::ptime b)
    {
    boost::posix_time::time_period tp(a, b);
    boost::posix_time::time_duration td = tp.length();
    return 0.001*td.total_milliseconds();
    }

inline void getFileContents(std::string & file_contents, std::string filename)
    {
    // Read contents of treefname into string
    std::ifstream f(filename.c_str());

    if (f.good())
        {
        f.seekg(0, std::ios::end);
        //long unsigned nbytes = f.tellg();
		std::streamoff nbytes = f.tellg();
        if (nbytes == 0)
            throw XGalax(boost::str(boost::format("File specified (%s) is empty") % filename));
        f.seekg(0, std::ios::beg);

        file_contents.insert( file_contents.begin(), std::istreambuf_iterator<char>(f), std::istreambuf_iterator<char>() );
        }
    else
        {
        throw XGalax(boost::str(boost::format("Could not read file %s") % filename));
        }

    }

inline void extractAllWhitespaceDelimitedStrings(std::vector<std::string> & receiving_vector, const std::string & text_to_process)
    {
    receiving_vector.clear();

    // Use regular expression to divide text_to_process into whitespace-delimited strings, which are stored in receiving_vector
    boost::regex re("\\s+");
    boost::sregex_token_iterator i(text_to_process.begin(), text_to_process.end(), re, -1);
    boost::sregex_token_iterator j;

    while(i != j)
        {
        receiving_vector.push_back(*i++);
        }
    }

inline bool parseTranslate(std::map< unsigned, std::string > & translate_map, const std::string & file_contents)
    {
    unsigned prev_ntaxa = (unsigned)translate_map.size();

    // Will either add to an empty translate_map or check values if translate_map has been created previously

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
        if (efficientAddOrCheck(translate_map, taxon_index, what[2].str()) == translate_map.end())
            throw XGalax(boost::str(boost::format("Taxon name (%s) does not match name already associated (%s) with taxon %d") % what[2].str() % translate_map[taxon_index] % taxon_index));
        }

    // Check to make sure this file does not have one or more additional taxa than files previously processed
    if (prev_ntaxa > 0 && prev_ntaxa < translate_map.size())
        return false;

    // Check to make sure this file does not have one or more fewer taxa than files previously processed
    if (taxa_seen.size() < translate_map.size())
       return false;

    return true;
    }

inline std::string stripComments(const std::string & commented_nexus_text)
    {
    boost::regex re("\\[.+?\\]");
    return boost::regex_replace(commented_nexus_text, re, std::string());
    }

inline void getNewicks(std::vector< std::string > & tree_descriptions, const std::string & file_contents, unsigned skip)
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

inline unsigned countNewickLeaves(const std::string newick)
    {
    boost::regex taxonexpr("[(,]\\s*(\\d+|\\S+?|['].+?['])\\s*(?=[,):])");
    boost::sregex_iterator m1(newick.begin(), newick.end(), taxonexpr);
    boost::sregex_iterator m2;
    return (unsigned)std::distance(m1, m2);
    }

inline std::pair<bool, unsigned> convertStringToUnsignedInteger(const std::string & theString)
    {
    // attempt to convert theString to an unsigned integer
    std::string s = theString;
    bool worked = true;
    boost::trim(s);
    unsigned x = 0;
    try
        {
        x = boost::lexical_cast<unsigned>(s);
        }
    catch(boost::bad_lexical_cast &)
        {
        // node name could not be converted to an integer value
        worked = false;
        }
    return std::pair<bool,unsigned>(worked, x);
    }

inline bool numberAlreadyUsed(unsigned x, std::set<unsigned> & used)
    {
    std::pair<std::set<unsigned>::iterator, bool> insert_result = used.insert(x);
    if (insert_result.second)
        {
        // insertion was made, so x has NOT already been used
        return false;
        }
    else
        {
        // insertion was not made, so set already contained x
        return true;
        }
    }

}

#endif /* defined(GALAX_GALAXUTIL_HPP) */
