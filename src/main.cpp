//
//  main.hpp
//  galax
//
//  Created by Paul O. Lewis on 10/25/14.
//  Copyright (c) 2014 Paul O. Lewis. All rights reserved.
//

#include <iostream>
#include <string>
#include <fstream>
#include <streambuf>
#include <boost/program_options.hpp>
#include <boost/format.hpp>
#include <boost/regex.hpp>
#include "galaxinfo.hpp"
#include "galax.hpp"

using namespace galax;

unsigned GalaxInfo::_sortby_index = 0;
std::string program_name = "galax";
unsigned major_version = 1;
unsigned minor_version = 2;
unsigned bugfix_version = 0;
std::string tree_file_name = "";
std::string list_file_name = "";
std::string output_file_name = "output-galax";
unsigned skipped_newicks = 0;
unsigned outgroup_taxon = 1;
bool trees_rooted = false;

void processCommandLineOptions(int argc, const char * argv[])
    {
    // Get user-supplied configuration settings from command line and/or config file
    boost::program_options::variables_map vm;
    boost::program_options::options_description desc("Allowed options");
    std::vector<std::string> to_pass_further;
    desc.add_options()
        ("help,h", "produce help message")
        ("version,v", "show program version")
        ("rooted,r", boost::program_options::bool_switch()->default_value(false), "expect trees to be rooted (leave out this option to assume unrooted)")
        ("details,d", boost::program_options::bool_switch()->default_value(false), "save information content details for each tree file")
        ("skip,s", boost::program_options::value<unsigned>(), "number of tree descriptions to skip at beginning of tree file")
        ("treefile,t", boost::program_options::value<std::string>(), "name of tree file in NEXUS format (used to generated a majority-rule consensus tree)")
        ("listfile,l", boost::program_options::value<std::string>(), "name of file listing whitespace-separated, NEXUS-formatted tree file names to be processed")
        ("outfile,o", boost::program_options::value<std::string>(), "file name prefix of output file to create (.txt extension will be added)")
        ("outgroup,g", boost::program_options::value<unsigned>(), "number of taxon to use as the outgroup (where first taxon listed in treefile translate statement is 1)")
        ;
    boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);
    try
        {
        const boost::program_options::parsed_options & parsed = boost::program_options::parse_config_file< char >("galax.conf", desc, false);  // 3rd. argument = allow unrecognized options
        to_pass_further = collect_unrecognized(parsed.options, boost::program_options::include_positional);
        boost::program_options::store(parsed, vm);
        }
    catch(boost::program_options::reading_file &)
        {
        std::cout << "Note: configuration file (galax.conf) not found" << std::endl;
        }
    boost::program_options::notify(vm);

    if (to_pass_further.size() > 1)
        {
        for (unsigned i = 0; i < to_pass_further.size() - 1; ++i)   // - 1 because this vector comprises alternating keys and values, so last element is a value
            {
            // If option name contains "unused.", then it is inside a group named "unused" in the config file, in which case the user wishes it to be ignored
            // Anything else that did not match, however, should generate a warning
            bool unused = !(to_pass_further[i].find("unused.") == std::string::npos);
            if (!unused)
                {
                std::cout << boost::str(boost::format("Warning: unrecognized option \"%s=%s\" (\"%s\" misspelled?)") % to_pass_further[i] % to_pass_further[i + 1] % to_pass_further[i]) << std::endl;
                }
            ++i;    // skip value field
            }
        }

    // If user used --help on command line, output usage summary and quit
    if (vm.count("help"))
        {
        std::cout << desc << "\n";
        std::exit(1);
        }

    // If user used --version on command line, output version and quit
    if (vm.count("version"))
        {
        std::cout << boost::str(boost::format("This is %s version %d.%d.%d") % program_name % major_version % minor_version % bugfix_version) << std::endl;
        std::exit(1);
        }

    // If user used --rooted on command line, expect trees in specified tree file to be rooted
    if (vm.count("rooted"))
        {
        trees_rooted = vm["rooted"].as<bool>();
        }

    // If user used --skip on command line, store number of tree descriptions to skip
    if (vm.count("skip"))
        {
        skipped_newicks = vm["skip"].as<unsigned>();
        }

    // If user used --outgroup on command line, store number of taxon to use as the outgroup
    if (vm.count("outgroup"))
        {
        outgroup_taxon = vm["outgroup"].as<unsigned>();
        if (outgroup_taxon < 1)
            {
            std::cout << "Outgroup taxon must be an integer greater than 0 and less than or equal to the number of taxa\n";
            std::exit(1);
            }
        }

    // Assuming user used --outfile on command line or outfile=<file name> in config file, store supplied file name prefix
    if (vm.count("outfile") > 0)
        {
        output_file_name = vm["outfile"].as<std::string>();
        }

    // Assuming user used --treefile on command line or treefile=<file name> in config file, store supplied file name
    if (vm.count("treefile") > 0)
        {
        tree_file_name = vm["treefile"].as<std::string>();
        }

    // Assuming user used --listfile on command line or listfile=<file name> in config file, store supplied file name
    if (vm.count("listfile") > 0)
        {
        list_file_name = vm["listfile"].as<std::string>();
        }

    // Check to make sure user specified either --treefile or --listfile
    if (list_file_name.size() == 0 && tree_file_name.size() == 0)
        {
        std::cout << "Must specify either --treefile or --listfile (or both) on command line or in config file\n";
        std::cout << desc << std::endl;
        std::exit(1);
        }

    // Output feedback to reassure user

    if (list_file_name.length() > 0 && tree_file_name.length() > 0)
        {
        std::cout << "Trees will be read from tree files specified in file " << list_file_name << std::endl;
        std::cout << "Majority-rule consensus will be constructed from trees in the file " << tree_file_name << std::endl;
        }
    else if (list_file_name.length() > 0)
        {
        std::cout << "Trees will be read from tree files specified in file " << list_file_name << std::endl;
        std::cout << "Majority-rule consensus will be constructed from the merged tree set" << std::endl;
        }
    else
        {
        std::cout << "Trees will be read from file " << tree_file_name << std::endl;
        std::cout << "Majority-rule consensus will be constructed from these trees" << std::endl;
        }

    std::cout << "Output will be stored in file " << output_file_name << std::endl;

    if (trees_rooted)
        std::cout << "Trees assumed to be rooted" << std::endl;
    else
        {
        std::cout << "Trees assumed to be unrooted" << std::endl;
        std::cout << "Outgroup for polarizing splits will be taxon number " << outgroup_taxon << std::endl;
        }

    if (skipped_newicks == 0)
        std::cout << "No trees will be skipped" << std::endl;
    else if (skipped_newicks == 1)
        std::cout << "The first tree will be skipped" << std::endl;
    else
        std::cout << "The first " << skipped_newicks << " trees will be skipped" << std::endl;
    }

int main(int argc, const char * argv[])
    {
    processCommandLineOptions(argc, argv);

    Galax(output_file_name).run(tree_file_name, list_file_name, skipped_newicks, trees_rooted, outgroup_taxon);

    return 0;
    }
