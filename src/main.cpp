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
#include "galax.hpp"

using namespace galax;

std::string program_name = "galax";
unsigned major_version = 1;
unsigned minor_version = 0;
std::string tree_file_name = "";
std::string output_file_name = "output-galax.txt";
unsigned skipped_newicks = 0;
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
        ("skip,s", boost::program_options::value<unsigned>(), "number of tree descriptions to skip at beginning of tree file")
        ("treefile,t", boost::program_options::value<std::string>(), "name of tree file in NEXUS format")
        ("outfile,t", boost::program_options::value<std::string>(), "name of output file to create")
        ;
    boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);
    try
        {
        const boost::program_options::parsed_options & parsed = boost::program_options::parse_config_file< char >("galax.conf", desc, false);  // 3rd. argument = allow unrecognized options
        to_pass_further = collect_unrecognized(parsed.options, boost::program_options::include_positional);
        boost::program_options::store(parsed, vm);
        }
    catch(boost::program_options::reading_file & x)
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
        std::cout << boost::str(boost::format("This is %s version %d.%d") % program_name % major_version % minor_version) << std::endl;
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

    // Assuming user used --outfile on command line or outfile=<file name> in config file, store supplied file name
    if (vm.count("outfile") > 0)
        {
        output_file_name = vm["outfile"].as<std::string>();
        }

    // Assuming user used --treefile on command line or treefile=<file name> in config file, store supplied file name
    if (vm.count("treefile") > 0)
        {
        tree_file_name = vm["treefile"].as<std::string>();
        }
    else
        {
        std::cout << "Must specify --treefile on command line or in config file\n";
        std::cout << desc << std::endl;
        std::exit(1);
        }

    std::cout << "Output will be stored in file " << output_file_name << std::endl;
    std::cout << "Trees will be read from file " << tree_file_name << std::endl;

    if (trees_rooted)
        std::cout << "Trees assumed to be rooted" << std::endl;
    else
        std::cout << "Trees assumed to be unrooted" << std::endl;

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

    Galax(output_file_name).run(tree_file_name, skipped_newicks, trees_rooted);

    return 0;
    }