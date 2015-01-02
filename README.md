#  Galax

Galax is a program that reads NEXUS-formatted tree files containing tree samples from a
Bayesian phylogenetic posterior distribution and estimates the total Lindley Information
assuming that a non-informative prior was used. It also partitions the total information
among clades.

Galax was written to accompany the following paper:

> Lewis, PO, MH Chen, L Kuo, LA Lewis, K Fucikova, S Neupane, YB Wang, D Shi. Estimating
Bayesian phylogenetic information content. Systematic Biology manuscript USYB-2014-293
currently under review.

Usage summary:
```
Must specify either --treefile or --listfile on command line or in config file
Allowed options:
  -h [ --help ]         produce help message
  -v [ --version ]      show program version
  -r [ --rooted ]       expect trees to be rooted (leave out this option to
                        assume unrooted)
  -d [ --details ]      save information content details for each tree file
  -s [ --skip ] arg     number of tree descriptions to skip at beginning of
                        tree file
  -t [ --treefile ] arg name of tree file in NEXUS format
  -l [ --listfile ] arg name of file listing whitespace-separated,
                        NEXUS-formatted tree file names to be processed
  -o [ --outfile ] arg  file name prefix of output file to create (.txt
                        extension will be added)
```
Galax was written by Paul O. Lewis.

 Version  |    Released   | Description
:-------: | :-----------: | :----------------------------------------------------
1.0       |  22-Dec-2014  | Estimates and partitions I across clades
1.1       |   2-Jan-2015  | Added --listfile and partitioning of D across clades
