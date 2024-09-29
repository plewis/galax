//
//  galax.cpp
//  galax
//
//  Created by Paul O. Lewis on 10/25/14.
//  Copyright (c) 2014 Paul O. Lewis. All rights reserved.
//

#include <fstream>
#include <limits>
#include <algorithm>
#include <sstream>
#include <regex>
#include <boost/bind.hpp>
#include <boost/format.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/algorithm/string/join.hpp>
#include "conditionals.hpp"
#include "galax.hpp"
#include "galaxutil.hpp"
#include "galaxinfo.hpp"
#include "galaxdata.hpp"
#include "xgalax.hpp"

using namespace galax;

const unsigned Galax::_ALLSUBSETS = std::numeric_limits<unsigned>::max();

Galax::Galax(const std::string outfile_prefix, const std::string version)
    {
    _alt_nexus = false;
    _alt_next_taxon_index = 0;
    _outfprefix = outfile_prefix;
    _version = version;
    }

Galax::~Galax()
    {
    }

void Galax::storeTrees(std::string file_contents, unsigned skip, std::vector< std::string > & tree_descriptions)
    {
    _start_time = getCurrentTime();

    if (isNexusFile(file_contents))
        {
        parseTranslate(file_contents);
        getNewicks(tree_descriptions, file_contents, skip);
        }
    else if (isRevBayesFile(file_contents))
        {
        getNewicksRevBayes(tree_descriptions, file_contents, skip);
        }
    else
        {
        throw XGalax("Galax currently requires tree files to be in either NEXUS format or RevBayes format.");
        //getPhyloBayesNewicks(tree_descriptions, file_contents, skip);

        //TODO: this is ugly - try to get rid of _reverse_translate by converting _translate to map(string, int) rather than map(int, string)
        //_translate.clear();
        //for (std::map< std::string, unsigned>::iterator it = _reverse_translate.begin(); it != _reverse_translate.end(); ++it)
        //    {
        //    _translate[it->second] = it->first;
        //    }
        }

    _end_time = getCurrentTime();
    _total_seconds += secondsElapsed(_start_time, _end_time);
    }

//void Galax::getTrees(std::string file_contents, unsigned skip)
//    {
//    _start_time = getCurrentTime();
//
//    std::vector< std::string > tree_descriptions;
//    if (isNexusFile(file_contents))
//        {
//        parseTranslate(file_contents);
//        getNewicks(tree_descriptions, file_contents, skip);
//        }
//    else
//        {
//        throw XGalax("Galax currently requires tree files to be in NEXUS format.");
//        //getPhyloBayesNewicks(tree_descriptions, file_contents, skip);
//
//        //TODO: this is ugly - try to get rid of _reverse_translate by converting _translate to map(string, int) rather than map(int, string)
//        //_translate.clear();
//        //for (std::map< std::string, unsigned>::iterator it = _reverse_translate.begin(); it != _reverse_translate.end(); ++it)
//        //    {
//        //    _translate[it->second] = it->first;
//        //    }
//        }
//    _tree_counts.push_back((unsigned)tree_descriptions.size());
//    _newicks.clear();
//    _newicks.insert(_newicks.end(), tree_descriptions.begin(), tree_descriptions.end());
//    _merged_newicks.insert(_merged_newicks.end(), tree_descriptions.begin(), tree_descriptions.end());
//
//    _end_time = getCurrentTime();
//    _total_seconds += secondsElapsed(_start_time, _end_time);
//    }

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

bool Galax::isNexusFile(const std::string & file_contents)
    {
    std::regex pattern("#nexus|#Nexus|#NEXUS");
    std::smatch what;
    std::string first_six_letters = file_contents.substr(0,6);
    bool is_nexus_file = (bool)std::regex_match(first_six_letters, what, pattern);
    return is_nexus_file;
    }

bool Galax::isRevBayesFile(const std::string & file_contents)
    {
    //  Iteration	Posterior	Likelihood	Prior	psi
    //  0	-852.3526	-904.9879	52.6353	(...newick...)[&index=102];
    //  1	-844.5589	-897.5716	53.01273	(...newick...)[&index=102];
    //  ...
    //  10000	-839.9447	-893.1031	53.15838	(...newick...)[&index=102];
    bool is_revbayes_file = false;
    std::regex re("Iteration\\s+Posterior\\s+Likelihood\\s+Prior\\s+psi", std::regex_constants::ECMAScript);
    std::sregex_iterator m1(file_contents.begin(), file_contents.end(), re);
    std::sregex_iterator m2;  // empty iterator used only to detect when we are done
    if (m1 != m2)
        is_revbayes_file = true;
    return is_revbayes_file;
    }

unsigned Galax::taxonNumberFromName(const std::string taxon_name, bool add_if_missing)
    {
    unsigned taxon_number = 0;
    std::map< std::string, unsigned >::iterator it = _taxon_map.find(taxon_name);
    if (it != _taxon_map.end())
        {
        // taxon_name is an existing key, so return its value
        taxon_number = it->second;
        }
    else
        {
        if (add_if_missing)
            {
            // taxon_name is not an existing key, so create an entry for it
            taxon_number = (unsigned)(_taxon_map.size() + 1);
            _taxon_map[taxon_name] = taxon_number;
            }
        else
            throw XGalax(boost::str(boost::format("Could not find %s in _taxon_map keys") % taxon_name));
        }

    return taxon_number;
    }

bool Galax::replaceTaxonNames(const std::string & newick_with_taxon_names, std::string & newick_with_taxon_numbers)
    {
    newick_with_taxon_numbers.clear();
    std::regex taxonexpr("[(,]+\\s*([^,):]+)");
    std::sregex_iterator m1(newick_with_taxon_names.begin(), newick_with_taxon_names.end(), taxonexpr);
    std::sregex_iterator m2;
    std::sregex_iterator last;
    for (; m1 != m2; ++m1)
        {
        std::smatch what = *m1;
        newick_with_taxon_numbers.append(what.prefix());
        newick_with_taxon_numbers.append(what[0].first, what[1].first) ;
        unsigned n = taxonNumberFromName(what[1].str(), true);
        if (n > 0)
            {
            std::ostringstream oss;
            oss << n;
            newick_with_taxon_numbers.append(oss.str());
            }
        else
            {
            throw XGalax(boost::str(boost::format("Could not convert taxon name (%s) into a valid taxon number") % what[1].str()));
            }
        last = m1;
        }
    std::smatch what = *last;
    newick_with_taxon_numbers.append(what.suffix());
    return true;
    }

void Galax::parseTranslate(const std::string & file_contents)
    {
    _translate.clear();

    // The first time this function is called, the _translate command
    // will be used to establish a _taxon_map such that
    // _taxon_map[<taxon name>] = <taxon index>. Subsequent calls
    // will create _translate anew but will not modify _taxon_map;
    // however, an exception is thrown if a taxon is about to be
    // added to _translate that is not a key in _taxon_map.
    
    // If there is no translate statement in the tree file, assume
    // that the tree descriptions contain names and not numbers.
    // In this case, _translate will be empty and _taxon_map will be
    // populated using the taxon names from the first newick tree
    // description read.
    bool translate_statement_found = false;

    // First, separate out the contents of the translate
    // statement from the rest of the tree file
    std::string translate_contents;
    std::regex pattern("[Tt]ranslate([\\S\\s]+?;)");
    std::smatch what;
    bool regex_ok = std::regex_search(file_contents, what, pattern);
    if (regex_ok)
        {
        translate_statement_found = true;
        
        // what[0] contains the whole string
        // what[1] contains the translate statement contents
                
        // Construct a string using characters in contents from what[1].first to what[1].second
        translate_contents.insert(translate_contents.begin(), what[1].first, what[1].second);
        }
    //else
    //    {
    //    throw XGalax("regex failed to find translate statement in tree file");
    //    }

    if (translate_statement_found)
        {
        // Now create the map by iterating through each element of the translate statement
        //std::regex re("(\\d+)\\s+'?(.+?)'?,?$");
        std::regex re("(\\d+)\\s+'?([\\s\\S]+?)'?\\s*[,;]");
        std::sregex_iterator m1(translate_contents.begin(), translate_contents.end(), re);
        std::sregex_iterator m2;
        for (; m1 != m2; ++m1)
            {
            const std::match_results<std::string::const_iterator>& what = *m1;
            
            //std::ofstream doof("doof_translate_element.txt");
            //doof << "taxon index = " << what[1] << std::endl;
            //doof << "taxon name  = " << what[2] << std::endl;
            //doof.close();
            //std::cout << "taxon index = " << what[1] << std::endl;
            //std::cout << "taxon name  = " << what[2] << std::endl;
            
            unsigned taxon_index = 0;
            try
                {
                taxon_index = boost::lexical_cast<unsigned>(what[1].str());
                }
            catch(const boost::bad_lexical_cast &)
                {
                throw XGalax("Could not interpret taxon index in translate statement as a number");
                }

            // efficientAddOrCheck returns valid iterator on first insertion or if identical association already previously made
            if (efficientAddOrCheck(_translate, taxon_index, what[2].str()) == _translate.end())
                throw XGalax(boost::str(boost::format("Taxon name (%s) does not match name already associated (%s) with taxon %d") % what[2].str() % _translate[taxon_index] % taxon_index));
            }
        
        typedef std::map< unsigned, std::string > translate_t;
        typedef std::map< std::string, unsigned > taxonmap_t;
        if (_taxon_map.empty())
            {
            // Copy _translate to _taxon_map
            BOOST_FOREACH(translate_t::value_type & p, _translate)
                {
                _taxon_map[p.second] = p.first;
                }
            }
        else
            {
            // Check to make sure this file does not have more or fewer taxa than files previously processed
            if (_translate.size() != _taxon_map.size())
               throw XGalax("Galax requires taxa to be the same for all tree files processed.");

            // Check to make sure all taxa in _translate are keys in _taxon_map
            BOOST_FOREACH(translate_t::value_type & p, _translate)
                {
                taxonmap_t::iterator it = _taxon_map.lower_bound(p.second);
                if (it == _taxon_map.end())
                    throw XGalax(boost::str(boost::format("Taxon name (%s) does not match any name stored for tree files read previously") % p.second));
                }
            }
            _alt_nexus = false;
        }
        else {
            _alt_nexus = true;
        }
    }

void Galax::getPhyloBayesNewicks(std::vector< std::string > & tree_descriptions, const std::string & file_contents, unsigned skip)
    {
    tree_descriptions.clear();
    //std::regex re("^(.+?);");
    std::regex re("([\\s\\S]+?);", std::regex_constants::ECMAScript);
    std::sregex_iterator m1(file_contents.begin(), file_contents.end(), re);
    std::sregex_iterator m2;  // empty iterator used only to detect when we are done
    unsigned n = 0;
    for (; m1 != m2; ++m1)
        {
        const std::match_results<std::string::const_iterator>& what = *m1;
        if (n >= skip)
            {
            std::string newick_only;
            replaceTaxonNames(what[1].str(), newick_only);
            tree_descriptions.push_back(newick_only);
            }
        n += 1;
        std::cerr << n << "\n";
        }
    }

std::string Galax::standardizeNodeNumber(std::smatch const & what)
    {
    std::string s;
    if (_alt_nexus) {
        std::string x = what[2].str();
        if (_taxon_map.count(x) == 0) {
            _taxon_map[x] = ++_alt_next_taxon_index;
            s = boost::str(boost::format("%s%d:") % what[1] % _alt_next_taxon_index);
        }
        else {
            unsigned stored_taxon_index = _taxon_map[x];
            s = boost::str(boost::format("%s%d:") % what[1] % stored_taxon_index);
        }
    }
    else {
        int x = atoi(what[2].str().c_str());

        // ensure that node numbers correspond to taxon indices stored in taxon_map,
        // which are not necessarily the same as those used in the translate command that
        // was in the tree file
        std::string & stored_taxon_name = _translate[x];
        assert(stored_taxon_name != "");
        unsigned stored_taxon_index = _taxon_map[stored_taxon_name];
        s = boost::str(boost::format("%s%d:") % what[1] % stored_taxon_index);
    }
    
    return s;
    }

std::string Galax::standardizeTreeDescriptionRevBayes(std::string & newick_in)
    {
    std::string stdnewick;
    std::regex re("([(,])([^:(,]+)", std::regex_constants::ECMAScript);
    std::sregex_iterator m1(newick_in.begin(), newick_in.end(), re);
    std::sregex_iterator m2;
    std::sregex_iterator last;
    int last_node_number = 1;
    bool create_taxon_map = (_taxon_map.size() == 0);
    for (; m1 != m2; ++m1)
        {
        std::smatch what = *m1;
        
        // //temporary!
        // std::cerr << "-------------" << std::endl;
        // std::cerr << "prefix: " << what.prefix().str() << std::endl;
        // for (auto & w : what) {
        //     std::cerr << w.str() << std::endl;
        // }
        
        std::string prefix = what.prefix().str();
        std::string before = what[1].str();
        std::string node_name = what[2].str();
        int node_number = 0;
        if (create_taxon_map) {
            // Ensure that same node name does not appear twice
            assert(_taxon_map.count(node_name) == 0);
            _translate[last_node_number] = node_name;
            _taxon_map[node_name] = last_node_number;
            node_number = last_node_number++;
        }
        else {
            node_number = _taxon_map.at(node_name);
        }
        
        stdnewick.append(prefix);
        stdnewick.append(before);
        stdnewick.append(std::to_string(node_number));
        
        last = m1;
        }
    std::smatch what = *last;
    stdnewick.append(what.suffix());
    return stdnewick;
    }

std::string Galax::standardizeTreeDescription(std::string & newick_in)
    {
    std::regex re("([(,])(\\d+):");
    if (_alt_nexus) {
        // Taxa in newick_in are in the form of names rather than numbers
        re = std::regex("([(,])([^:]):");
    }
    
    std::string stdnewick;
    std::sregex_iterator m1(newick_in.begin(), newick_in.end(), re);
    std::sregex_iterator m2;
    std::sregex_iterator last;
    for (; m1 != m2; ++m1)
        {
        std::smatch what = *m1;
        
        //std::cerr << "what.prefix() = " << what.prefix() << std::endl;
        //std::cerr << "what[1]       = " << what[1] << std::endl;
        //std::cerr << "what[2]       = " << what[2] << std::endl;
        //std::cerr << "what.suffix() = " << what.suffix() << std::endl;
        
        stdnewick.append(what.prefix());
        stdnewick.append(standardizeNodeNumber(what));
        
        last = m1;
        }
    std::smatch what = *last;
    stdnewick.append(what.suffix());

    return stdnewick;
    }

//std::string Galax::standardizeTreeDescription(std::string & newick_in)
//    {
//    // Search for node numbers in newick_in and replace them with standard node numbers stored in _taxon_map
//    std::regex re("([(,])(\\d+)");
//    boost::function<std::string (std::smatch const &)> function = boost::bind(&Galax::standardizeNodeNumber, this, _1);
//    std::string newick_out = std::regex_replace(newick_in, re, function);
//    return newick_out;
//    }

void Galax::getNewicks(std::vector< std::string > & tree_descriptions, const std::string & file_contents, unsigned skip)
    {
    // tree STATE_0 [&lnP=-4493.80476846934,posterior=-4493.80476846934] = [&R] (9:[&rate=0.49971158909783764]1851.4724462198697,((1:[&rate=0.5965730394621352]292.73199858783727,(10:[&rate=0.6588031360335018]30.21172743645451,5:[&rate=0.7098036299867017]30.21172743645451):[&rate=1.0146941544458208]262.52027115138276):[&rate=1.0649642758561977]510.452441519872,((8:[&rate=0.7554924641162211]145.1076605992074,(7:[&rate=0.7984750329147966]64.0435017480143,6:[&rate=0.8402528958963882]64.0435017480143):[&rate=1.1206854064651213]81.06415885119311):[&rate=1.1844450597679457]522.823827314411,((3:[&rate=0.8818808237868384]60.2962343089954,4:[&rate=0.9242396697890951]60.2962343089954):[&rate=1.260685743226102]12.793802911399744,2:[&rate=0.9681896872253556]73.09003722039515):[&rate=1.3582802932633053]594.8414506932232):[&rate=1.4999660689010508]135.25295219409088):[&rate=1.7907115550989796]1048.2880061121605);
    // tree STATE_0 [&lnP=-222468.46405040708,posterior=-222468.46405040708] = [&R] (((((((13:[&rate=0.41948130885774293]6.2959114346539975,((7:[&rate=0.4835763823640333]0.13550093579331035,1:[&rate=0.5209684767630073]0.13550093579331035):[&rate=0.9635652545109973]0.5551042325454463,15:[&rate=0.5492827727849298]0.6906051683387566):[&rate=0.9755113666282312]5.605306266315241):[&rate=0.9876472559631831]10.322795910985544,(17:[&rate=0.5728193948810578]5.685111439869553,2:[&rate=0.5933667775296758]5.685111439869553):[&rate=0.9999937594617617]10.93359590576999):[&rate=1.0125732516497858]5.142489242495607,(5:[&rate=0.6118563212473699]0.15949221153379095,28:[&rate=0.6288406629410516]0.15949221153379095):[&rate=1.0254099216508532]21.601704376601358):[&rate=1.038530094409559]116.24314628707417,(((29:[&rate=0.6446773336086477]0.762443175497069,8:[&rate=0.6596124676296977]0.762443175497069):[&rate=1.051962607302118]22.93607986063065,(23:[&rate=0.6738236879118028]0.5047866135985948,21:[&rate=0.6874440397483901]0.5047866135985948):[&rate=1.0657392562843637]23.193736422529124):[&rate=1.0798953297143943]16.280424093309122,20:[&rate=0.7005762553101323]39.97894712943684):[&rate=1.0944702533868567]98.02539574577247):[&rate=1.1095083777014945]20.87666238489902,18:[&rate=0.713301711000359]158.88100526010834):[&rate=1.1250599481079797]25.369142391234874,(((((6:[&rate=0.7256862972167103]4.157046671443346,27:[&rate=0.7377844041897827]4.157046671443346):[&rate=1.1411823142956972]0.4216380394703174,31:[&rate=0.749641711811296]4.578684710913663):[&rate=1.157941453976706]2.7888871791074665,((14:[&rate=0.7612971942888511]0.4820008033421106,26:[&rate=0.7727845943689399]0.4820008033421106):[&rate=1.175413916572047]2.1006906770419045,33:[&rate=0.7841335302910268]2.582691480384015):[&rate=1.1936893354723344]4.784880409637115):[&rate=1.2128737226386077]28.320425658481167,(4:[&rate=0.7953703429930343]9.047384492551547,19:[&rate=0.8065187562334204]9.047384492551547):[&rate=1.2330938592149618]26.640613055950748):[&rate=1.254503252991522]144.44188595962413,(34:[&rate=0.8176003998695023]49.84530458500085,(((35:[&rate=0.8286352317599452]4.64783851774415,9:[&rate=0.8396418838264186]4.64783851774415):[&rate=1.2772903877797566]17.7189245645138,(12:[&rate=0.8506379510097704]9.214275524579357,(3:[&rate=0.8616402371306506]6.8925402326469944,32:[&rate=0.8726649683413712]6.8925402326469944):[&rate=1.301690414163583]2.321735291932362):[&rate=1.3280021656518173]13.152487557678594):[&rate=1.356613709910537]7.407288154720803,((30:[&rate=0.8837279825004343]3.4163490119965814,25:[&rate=0.8948449011281925]3.4163490119965814):[&rate=1.3880421571600987]13.604316647910641,24:[&rate=0.9060312894219462]17.020665659907223):[&rate=1.422998494768541]12.75338557707153):[&rate=1.4624990961768445]20.071253348022093):[&rate=1.5080711470564487]130.2845789231256):[&rate=1.5621665830644125]4.120264143216787):[&rate=1.6291050095698845]125.24494247082069,(16:[&rate=0.9173028089947355]35.446341908327426,((22:[&rate=0.9286753674700471]0.013463304338836088,11:[&rate=0.940165268760073]0.013463304338836088):[&rate=1.7176458014779916]7.67054782472139,10:[&rate=0.9517893784722757]7.6840111290602255):[&rate=1.8504611669408024]27.762330779267202):[&rate=2.133204264216236]274.04874821383646);
    // tree STATE_0 = ((((1:0.015159374158485823,6:0.015159374158485823):0.0064804882886747035,2:0.021639862447160527):0.104324463428508,5:0.12596432587566853):0.30645174794996116,(3:0.4084347373105321,4:0.4084347373105321):0.023981336515097595):0.0;
    
    tree_descriptions.clear();
    //std::regex re("^\\s*[Tt]ree([\\s\\S]+?);");
    std::regex re("[Uu]*[Tt]ree\\s([\\s\\S]+?);", std::regex_constants::ECMAScript);
    std::sregex_iterator m1(file_contents.begin(), file_contents.end(), re);
    std::sregex_iterator m2;  // empty iterator used only to detect when we are done
    unsigned n = 0;
    for (; m1 != m2; ++m1)
        {
        const std::match_results<std::string::const_iterator>& what = *m1;
        if (n >= skip)
            {
            //std::ofstream doof("doof_tree_match.txt");
            //doof << what[1].str() << std::endl;
            //doof.close();
            
            std::string stripped = stripComments( what[1].str() );

            //doof.open("doof_stripped.txt");
            //doof << stripped << std::endl;
            //doof.close();
            
            std::string newick_only = stripTreeName(stripped);

            //doof.open("doof_newick_only.txt");
            //doof << newick_only << std::endl;
            //doof.close();
            
            std::string std_newick = standardizeTreeDescription(newick_only);

            //doof.open("doof_std_newick.txt");
            //doof << std_newick << std::endl;
            //doof.close();
            
            tree_descriptions.push_back(std_newick);
            }
        n += 1;
        }
    }

void Galax::getNewicksRevBayes(std::vector< std::string > & tree_descriptions, const std::string & file_contents, unsigned skip)
    {
    //  Iteration	Posterior	Likelihood	Prior	psi
    //  0	-852.3526	-904.9879	52.6353	(...newick...)[&index=102];
    //  1	-844.5589	-897.5716	53.01273	(...newick...)[&index=102];
    //  ...
    //  10000	-839.9447	-893.1031	53.15838	(...newick...)[&index=102];
    tree_descriptions.clear();
    std::regex re("\\d+\\s+[-.0-9e]+\\s+[-.0-9e]+\\s+[-.0-9e]+\\s+(\\S+);", std::regex_constants::ECMAScript);
    std::sregex_iterator m1(file_contents.begin(), file_contents.end(), re);
    std::sregex_iterator m2;  // empty iterator used only to detect when we are done
    unsigned n = 0;
    for (; m1 != m2; ++m1)
        {
        const std::match_results<std::string::const_iterator>& what = *m1;
        if (n >= skip)
            {
            //std::cerr << what[1].str() << std::endl;
            
            std::string newick_only = stripComments( what[1].str() );
            //std::cerr << newick_only << std::endl;
            
            std::string std_newick = standardizeTreeDescriptionRevBayes(newick_only);
            //std::cerr << std_newick << std::endl;
            
            tree_descriptions.push_back(std_newick);
            }
        n += 1;
        }
    }

void Galax::processTrees(TreeManip<Node>::TreeManipShPtr tm, CCDMapType & ccdmap, unsigned subset_index, unsigned num_subsets)
    {
    _start_time = getCurrentTime();

#if defined(POLNEW)
    // Create tree set for use in computing raw posterior entropy
    std::set<Split> tree_id;
    std::map<std::set<Split>, unsigned> topo_map;
#endif

    for (std::vector< std::string >::const_iterator sit = _newicks.begin(); sit != _newicks.end(); ++sit)
        {
        try
            {
            //std::ofstream doof("doof_newick.txt");
            //doof << *sit << std::endl;
            //doof.close();
            
            tm->buildFromNewick(*sit, (_rooted ? 0 : _outgroup));
#if defined(POLNEW)
            tm->calcTreeID(tree_id);
            if (topo_map.count(tree_id) > 0) {
                topo_map[tree_id]++;
            }
            else {
                topo_map[tree_id] = 1;
            }
#endif
            tm->addToCCDMap(ccdmap, _treeCCD, _treeMap, _show_details, subset_index, num_subsets);
            }
        catch(XGalax & x)
            {
            std::cerr << "ERROR: " << x.what() << std::endl;
            std::exit(1);
            }
        }

#if defined(POLNEW)
    // Calculate raw entropy for this subset
    double xlogx = 0.0;
    double n = 0.0;
    for (auto & p : topo_map) {
        double x = (double)p.second;
        n += x;
        xlogx += x*log(x);
    }
    double max_entropy = log(n);
    double raw_entropy = log(n) - xlogx/n;
    double raw_Ipct = 100.0*(max_entropy - raw_entropy)/max_entropy;
    
    unsigned n2 = (unsigned)_translate.size();
    double max_entropy2 = lognrooted(_rooted ? n2 : n2 - 1);
    double raw_Ipct2 = 100.0*(max_entropy2 - raw_entropy)/max_entropy2;
    
    _outf << std::endl;
    _outf << "Subset " << subset_index << std::endl;
    _outf << boost::str(boost::format("  %12.5f raw posterior entropy\n") % raw_entropy);
    _outf << boost::str(boost::format("  %12d sample size\n") % (unsigned)n);
    _outf << boost::str(boost::format("  %12.5f maximum entropy given sample size\n") % max_entropy);
    _outf << boost::str(boost::format("  %12.5f maximum entropy given total num. topologies\n") % max_entropy2);
    _outf << boost::str(boost::format("  %12.5f percent raw information given sample size\n") % raw_Ipct);
    _outf << boost::str(boost::format("  %12.5f percent raw information given total num. topologies\n") % raw_Ipct2);
    _outf << std::endl;
#endif

    _end_time = getCurrentTime();
    _total_seconds += secondsElapsed(_start_time, _end_time);
    }

void Galax::writeInfoProfile(TreeManip<Node>::TreeManipShPtr tm, GalaxInfoVector & clade_info)
    {
    _start_time = getCurrentTime();

    double sample_size = (double)_newicks.size();

    std::vector< double > profile_times;
    std::vector< double > profile_weights;

    // Create a map from clade_info that provides the information content for any given clade
    std::map< std::string, double > clade_map;
    for (GalaxInfoVector::const_iterator gi = clade_info.begin(); gi != clade_info.end(); ++gi)
        {
        const GalaxInfo & ginfo = *gi;
        double Iprop =ginfo._value[0]/100.0;
        double cladeprob = ginfo._value[2];
        double wt = Iprop/(cladeprob*sample_size);
        clade_map[ginfo._name] = wt;
        }
    std::cerr << "clade_map.size() = " << clade_map.size() << std::endl;

    // Build all trees, storing the time and information content (proportion of maximum possible) for each clade in profile_times and profile_weights, respectively
    for (std::vector< std::string >::const_iterator sit = _newicks.begin(); sit != _newicks.end(); ++sit)
        {
        try
            {
            tm->buildFromNewick(*sit, (_rooted ? 0 : _outgroup));
            tm->addToProfile(profile_times, profile_weights, clade_map);
            }
        catch(XGalax & x)
            {
            std::cerr << "ERROR: " << x.what() << std::endl;
            std::exit(1);
            }
        //break;  //temporary!
        }

    // normalize profile weights
    double sum_weights = std::accumulate(profile_weights.begin(), profile_weights.end(), 0.0);
    std::cerr << "sum_weights = " << sum_weights << std::endl;
    std::transform(profile_weights.begin(), profile_weights.end(), profile_weights.begin(), boost::lambda::_1/sum_weights);

    /*
    This is an example in R of weighted kernel density estimation:

    x <- seq(-1,6,.01)  # this is a list of 7*100=700 points along the x axis to use for plotting lines
    y1 <- 1             # first data point
    y2 <- 4             # second data point
    y3 <- 5             # third data point
    w1 <- 2/3           # weight for first data point
    w2 <- 1/6           # weight for second data point
    w3 <- 1/6           # weight for third data point
    h <- 2.0            # the bandwidth (standard deviation of component kernel densities)

    # plot the kernel density estimate using the black-box function "density"
    plot(density(c(y1,y2,y3), bw=h, kernel="gaussian", weights=c(w1,w2,w3)), xlim=c(-1,6), ylim=c(0,0.5), lwd=10, col="gray", lty="solid")

    # plot the three component kernel densities using blue, dotted lines
    curve(dnorm((x-y1)/h), col="blue", lty="dotted", add=T)
    curve(dnorm((x-y2)/h), col="blue", lty="dotted", add=T)
    curve(dnorm((x-y3)/h), col="blue", lty="dotted", add=T)

    # plot the kernel density estimate again (but do not use the "density" function this time) using a thicker, red, dotted line
    #curve(w1*dnorm((x-y1)/h)/h + w2*dnorm((x-y2)/h)/h + w3*dnorm((x-y3)/h)/h, add=T, lty="dotted", col="red", lwd=2)
    curve(w1*dnorm(x,y1,h) + w2*dnorm(x,y2,h) + w3*dnorm(x,y3,h), add=T, lty="dotted", col="red", lwd=2)
    */

    std::ofstream profile("profile.R");
    unsigned n = (unsigned)profile_times.size();
    assert(n == (unsigned)profile_weights.size());
    profile << "t <- c(" << join(profile_times.begin(), profile_times.end(), ",") << ")" << std::endl;
    profile << "w <- c(" << join(profile_weights.begin(), profile_weights.end(), ",") << ")" << std::endl;
    profile << "density(t, weights=w, from=0, give.Rkern=T)" << std::endl;
    profile << "plot(density(t, weights=w, from=0, adjust=1), main=\"Information Profile\", xlab=\"Node Height\", ylab=\"\")" << std::endl;
    profile << "rug(t)" << std::endl;
    profile.close();

    std::ofstream pro_file("profile.txt");
    for (unsigned i = 0; i < n; ++i)
        {
        pro_file << profile_times[i] << '\t' << profile_weights[i] << std::endl;
        }
    profile.close();

    _end_time = getCurrentTime();
    _total_seconds += secondsElapsed(_start_time, _end_time);
    }

void Galax::debugShowCCDMap(CCDMapType & ccdmap, unsigned subset_index)
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

void Galax::estimateInfo(CCDMapType & ccdmap, std::string & summaryinfostr, std::string & detailedinfostr, GalaxInfoVector & clade_info)
    {
    // ccdmap: key = SplitVector defining (un)conditional clade, value = CountVector (counts for each subset)
    // summaryinfostr: holds summary output over all clades
    // detailedinfostr: holds detailed output for each clade
    // clade_info: vector of GalaxInfo objects that is filled here and used to build and label majority rule tree later
    
    _start_time = getCurrentTime();

    assert(ccdmap.size() > 0);
    Split clade;
    clade_info.clear();
    bool first = true;
    std::vector<double> tmp;

    //std::cerr << "***** In Galax::estimateInfo, _outgroup = " << (_rooted ? 0 : _outgroup) << std::endl;

    // Create a GalaxData object to do most of the work
    unsigned ntaxa = (unsigned)_translate.size();
    GalaxData gd(_tree_counts, _treefile_names, (_rooted ? ntaxa : ntaxa - 1), (_rooted ? 0 : _outgroup));
    if (_show_details)
        gd.setShowDetails(true);
    else
        gd.setShowDetails(false);

    //unsigned total_trees = std::accumulate(_tree_counts.begin(), _tree_counts.end(), 0);

    // Provide a key to the taxa in split representations in both output files
    std::string tmpstr = "Order of taxa in split representations:\n";
    typedef std::pair<unsigned, std::string> translate_map_type;
    BOOST_FOREACH(const translate_map_type & translate_key_value, _translate)
        {
        tmpstr += boost::str(boost::format("%12d  %s\n") % translate_key_value.first % translate_key_value.second);
        }
    summaryinfostr = tmpstr;
    detailedinfostr = tmpstr;
    detailedinfostr += "\nLindley Information\n";
    detailedinfostr += boost::str(boost::format("  Number of taxa in ingroup: %d\n") % ntaxa);
    detailedinfostr += boost::str(boost::format("  Total prior entropy: %.5f\n") % gd.getTotalEntropy());
    detailedinfostr += "\nKey to columns:\n\n";
    detailedinfostr += "  treefile = file containing sampled trees from one data subset\n";
    detailedinfostr += "  trees    = number of trees sampled\n";
    detailedinfostr += "  K        = Kullback-Leibler divergence for current clade\n";
    detailedinfostr += "  p*       = marginal posterior probability of current clade\n";
    detailedinfostr += "  I        = information for current clade (p_C*KL_c)\n";
    detailedinfostr += "  Ipct     = information for current clade as percent of maximum information\n";
    detailedinfostr += "  D        = component of dissonance specific to current clade\n\n";

    for (CCDMapType::iterator it = ccdmap.begin(); it != ccdmap.end(); ++it)
        {
        // *it comprises a key (SplitVector) and a value (CountVector)
        // it->first is the key: it is a vector containing either
        //       a) 1 split (represents an unconditional clade probability) or
        //       b) 3 splits (represents a conditional clade probability, with
        //          first split being parent clade and remaining two splits
        //          representing the bipartition of that parent split)
        // it->second is the value: it is a vector of counts, with one element for each subset
        //const SplitVector   & v     = it->first;
        //std::vector<double> & count = it->second;
        //double sum_counts = std::accumulate(count.begin(), count.end(), 0.0);

        // The following depends on ccdmap keys (split vectors) being sorted such that an
        // entry for an unconditional clade, e.g. (ABC), precedes any conditional clade entries,
        // e.g. (ABC,A,BC) that have that clade as parent. Thus, when next unconditional
        // clade entry is encountered, we know that there are no more conditional clade entries
        // for the previous parent clade and we can at that point compute the information for
        // the previous parent.
        if (first)
            {
            first = false;
            clade = (it->first)[0];
            gd.newClade(it->first, it->second);
            }
        else
            {
            if ((it->first)[0] == clade)
                {
                // Conditional clade entry whose parent is clade
                gd.conditionalClade(it->first, it->second);
                }
            else
                {
                // Just found next unconditional clade entry, so now is the time to finish computing
                // information content for the previous clade

                tmp = gd.finalizeClade(detailedinfostr);
                gd.newClade(it->first, it->second);

                // tmp[0] = Ipct
                // tmp[1] = D (not Dpct)
                // tmp[2] = P
                // tmp[3] = I
                clade_info.push_back(GalaxInfo(clade.createPatternRepresentation(), tmp));
                clade = (it->first)[0];
                }
            }
        }

    // Compute I for final clade
    tmp = gd.finalizeClade(detailedinfostr);
    clade_info.push_back(GalaxInfo(clade.createPatternRepresentation(), tmp));

    gd.estimateCoverage(ccdmap, _treeCCD, _treeMap, detailedinfostr);
    gd.reportTotals(summaryinfostr, clade_info);

    _end_time = getCurrentTime();
    _total_seconds += secondsElapsed(_start_time, _end_time);
    }

void Galax::mapToTree(std::string & maptofname, GalaxInfoVector & annotate_info, std::string & mapto_newick)
    {
    _start_time = getCurrentTime();

    _rooted = _mapto_rooted;

    // Read tree file maptofname and store newick tree descriptions therein
    mapto_newick.clear();
    std::string file_contents;
    getFileContents(file_contents, maptofname);
    std::vector< std::string > newicks;
    storeTrees(file_contents, 0, newicks);

    // Build the first newick tree description
    TreeManip<Node>::TreeManipShPtr tm(new TreeManip<Node>());
    tm->buildFromNewick(newicks[0], (_rooted ? 0 : _outgroup));

    // Ensure that annotate_info is sorted so that the clades with highest posterior are first
    // This is important only for calculation of IC
    GalaxInfo::_sortby_index = 2;   // 0=Ipct, 1=D, 2=w, 3=I
    std::sort(annotate_info.begin(), annotate_info.end(), std::greater<GalaxInfo>());

    // Walk through tree calling setInfo, setWeight, setDisparity, and setCertainty for each internal node
    tm->annotateTree(annotate_info);

    // Create annotated newick tree description
    mapto_newick = tm->makeNewick(5, true);

    _end_time = getCurrentTime();
    _total_seconds += secondsElapsed(_start_time, _end_time);
    }

void Galax::buildMajorityRuleTree(GalaxInfoVector & majrule_info, GalaxInfoVector & annotate_info, std::string & majrule_newick)
    {
    _start_time = getCurrentTime();

    majrule_newick.clear();

    // We can save a lot of computation if both info vectors are the same
    bool same_info = (&majrule_info == &annotate_info);

    // Ensure that both majrule_info and annotate_info are sorted so that the clades with highest posterior are first
    GalaxInfo::_sortby_index = 2;   // 0=Ipct, 1=D, 2=w, 3=I
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

        if (w <= 0.5)
            {
            // majrule_info is sorted by w, so we can stop once we pass the 0.5 point
            // Note the less-than-or-equals sign above: using just less-than will get you into trouble
            // if there are muliple conflicting clades that have exactly 50% support.
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
                //TODO: this can be done more efficiently using the conditional clade distribution
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
        }

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

    // Output the translate command
    _treef << "#nexus\n\nbegin trees;\n  translate\n";
    typedef std::map<std::string, unsigned > taxon_map_type;
#if 0
    taxon_map_type::iterator it = _taxon_map.begin();
    _treef << "  " << it->second << " '" << it->first << "'";
    for (; it != _taxon_map.end(); ++it)
        {
        _treef << ",\n  " << it->second << " '" << it->first << "'";
        }
#else
    std::vector<std::string> translate_elements(_taxon_map.size());
    BOOST_FOREACH(taxon_map_type::value_type & p, _taxon_map)
        {
        translate_elements[p.second - 1] = boost::str(boost::format("    %d '%s'") % p.second % p.first);
        }
    _treef << boost::algorithm::join(translate_elements, ",\n") << ";" << std::endl;
#endif

    // Output the tree command
    _treef << "tree majrule = " << majrule_newick << ";\nend;" << std::endl;
    _treef.close();

    _end_time = getCurrentTime();
    _total_seconds += secondsElapsed(_start_time, _end_time);
    }

void Galax::initTreeCCD(unsigned num_subsets)
    {
    // initialize _treeCCD, which will hold tree IDs for all unique tree topologies in each subset
    _treeCCD.clear();
    _treeMap.clear();
    for (unsigned i = 0; i <= num_subsets; ++i)
        {
        TreeIDSetType v;
        _treeCCD.push_back(v);
        TreeMapType m;
        _treeMap.push_back(m);
        }
    }

void Galax::run(std::string treefname, std::string listfname, unsigned skip, bool trees_rooted, std::string maptofname, bool mapto_trees_rooted, bool save_details, unsigned outgroup_taxon)
    {
    assert (_ALLSUBSETS == (unsigned)(-1));
	try
		{
        std::string outfname = std::string(_outfprefix + ".txt");
        _outf.open(outfname.c_str());

        // Start by reporting settings used
        _outf << "Galax " << _version << "\n\n";
        _outf << "Options specified:\n";
        _outf << "   --treefile: " << treefname << "\n";
        _outf << "   --listfile: " << listfname << "\n";
        _outf << "       --skip: " << skip << "\n";
        _outf << "     --rooted: " << (trees_rooted ? "true" : "false") << "\n";
        _outf << "    --details: " << (save_details ? "true" : "false") << "\n";
        _outf << "   --outgroup: " << outgroup_taxon << "\n";
        _outf << "    --outfile: " << _outfprefix << "\n";
        _outf << "      --mapto: " << maptofname << "\n";
        _outf << "--maptorooted: " << (mapto_trees_rooted ? "true" : "false") << "\n";
        _outf << std::endl;

        _show_details = save_details;
        _outgroup = outgroup_taxon;
        _input_rooted = trees_rooted;
        _mapto_rooted = mapto_trees_rooted;
        _rooted = _input_rooted;
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

        std::string file_contents;
        std::string summaryinfostr;
        std::string detailedinfostr;
        TreeManip<Node>::TreeManipShPtr tm(new TreeManip<Node>());

        if (is_treefile)
            {
            //
            // --treefile provided on command line
            //
            initTreeCCD(1);
            _ccdtree.clear();
            getFileContents(file_contents, treefname);
            std::vector< std::string > tree_descriptions;
            storeTrees(file_contents, skip, tree_descriptions);
            _tree_counts.push_back((unsigned)tree_descriptions.size());
            _newicks.clear();
            _newicks.insert(_newicks.end(), tree_descriptions.begin(), tree_descriptions.end());
            _merged_newicks.insert(_merged_newicks.end(), tree_descriptions.begin(), tree_descriptions.end());
            processTrees(tm, _ccdtree, 0, 1);

            std::string msg = boost::str(boost::format("Read %d trees from tree file %s\n") % _newicks.size() % treefname);
            std::cout << msg;
            _outf << msg;
            
            if (_newicks.size() == 0)
                throw XGalax("Quitting because no trees were found to process.");

            GalaxInfoVector majrule_clade_info;
            std::vector<Split> majrule_splits;
            std::string majrule_tree;
            estimateInfo(_ccdtree, summaryinfostr, detailedinfostr, majrule_clade_info);
            buildMajorityRuleTree(majrule_clade_info, majrule_clade_info, majrule_tree);
            writeMajruleTreefile("majrule", majrule_tree);

            //std::string mapto_tree;
            //mapToTree(majrule_clade_info, majrule_clade_info, mapto_tree);
            //writeMajruleTreefile("@@@", mapto_tree);

            //if (_rooted)
            //    {
            //    // profiling informativeness only makes sense for ultrametric trees
            //    writeInfoProfile(tm, majrule_clade_info);
            //    }
            }
        else
            {
            //
            // --listfile provided on command line
            //
            std::cout << boost::str(boost::format("Reading tree file names from list file %s\n") % listfname);

            unsigned subset_index = 0;
            unsigned num_subsets = (unsigned)_treefile_names.size();
            initTreeCCD(num_subsets);
            _ccdlist.clear();
            BOOST_FOREACH(std::string & tree_file_name, _treefile_names)
                {
                file_contents.clear();
                getFileContents(file_contents, tree_file_name);
                std::vector< std::string > tree_descriptions;
                storeTrees(file_contents, skip, tree_descriptions);
                _tree_counts.push_back((unsigned)tree_descriptions.size());
                _newicks.clear();
                _newicks.insert(_newicks.end(), tree_descriptions.begin(), tree_descriptions.end());
                _merged_newicks.insert(_merged_newicks.end(), tree_descriptions.begin(), tree_descriptions.end());
                processTrees(tm, _ccdlist, subset_index++, (unsigned)_treefile_names.size());
                std::string msg = boost::str(boost::format("Read %d trees from tree file %s\n") % _newicks.size() % tree_file_name);
                std::cout << msg;
                _outf << msg;
                }
            GalaxInfoVector merged_clade_info;
            std::vector<Split> merged_splits;
            std::string merged_tree;
            estimateInfo(_ccdlist, summaryinfostr, detailedinfostr, merged_clade_info);
            buildMajorityRuleTree(merged_clade_info, merged_clade_info, merged_tree);
            writeMajruleTreefile("merged", merged_tree);

            if (maptofname.size() > 0)
                {
                std::string mapto_tree;
                mapToTree(maptofname, merged_clade_info, mapto_tree);
                writeMajruleTreefile("mapped", mapto_tree);
                }

            }

        if (_rooted)
            std::cout << "Input trees assumed to be rooted\n";
        else
            {
            std::cout << "Input trees assumed to be unrooted\n";
            std::cout << boost::str(boost::format("Each input tree was rooted at outgroup taxon %d (\"%s\")\n") % _outgroup % _translate[_outgroup]);
            }

        _outf << "\n" << summaryinfostr;
        _outf.close();

        if (save_details)
            {
            std::string detailsfname = std::string(_outfprefix + "-details.txt");
            _detailsf.open(detailsfname.c_str());
            _detailsf << "\n" << detailedinfostr;
            _detailsf.close();
            }

        std::cout << "\nRequired " << _total_seconds << " total seconds" << std::endl;
        }
	catch(XGalax x)
		{
		std::cerr << x.what() << std::endl;
		}
    }

