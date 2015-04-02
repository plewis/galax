#  Galax

Galax is a program that reads NEXUS-formatted tree files containing tree samples from a
Bayesian phylogenetic posterior distribution and estimates the total Lindley Information
assuming that a non-informative prior was used. It also partitions the total information
among clades.

Galax was written by Paul O. Lewis to accompany the following paper:

> Lewis, PO, MH Chen, L Kuo, LA Lewis, K Fucikova, S Neupane, YB Wang, D Shi. Estimating
Bayesian phylogenetic information content. Systematic Biology manuscript USYB-2014-293
currently under review.

## Usage summary
```
Must specify either --treefile or --listfile (or both) on command line or in config file
Allowed options:
  -h [ --help ]         produce help message
  -v [ --version ]      show program version
  -r [ --rooted ]       expect trees to be rooted (leave out this option to
                        assume unrooted)
  -d [ --details ]      save information content details for each tree file
  -s [ --skip ] arg     number of tree descriptions to skip at beginning of
                        tree file
  -t [ --treefile ] arg name of tree file in NEXUS format (used to generate a
                        majority-rule consensus tree)
  -l [ --listfile ] arg name of file listing whitespace-separated,
                        NEXUS-formatted tree file names to be processed
  -o [ --outfile ] arg  file name prefix of output file to create (.txt
                        extension will be added)
  -g [ --outgroup ] arg number of taxon to use as the outgroup (where first
                        taxon listed in treefile translate statement is 1)
```

## Usage notes

### Analyzing a single tree file

Running Galax as follows `galax --treefile trees.t` results in a majority rule tree being saved in the file *output-galax-majrule.tre* and a summary of phylogenetic information in the file *output-galax.txt*. The majority rule tree file contains information about each node that can be viewed in FigTree (http://tree.bio.ed.ac.uk/software/figtree/) by checking "Node Labels". For each clade seen in any tree in the file *trees.t*, Galax outputs a clade descriptor (sequence of hyphens and asterisks indicating which taxa are included in the clade) and values for each of the following quantities:

| Quantity |  Description                                                 |
|:---------|:-------------------------------------------------------------|
|   Hp     | Clade specific prior entropy                                 |
|   H      | Clade specific posterior entropy                             |
|   Hp - H | Clade specific unweighted information                        |
|   w      | Marginal clade posterior probability (clade weight)          |
|   Ipct   | I expressed as a percentage of the total                     |
|   D      | Disparity (not applicable for analyses of single tree files) |


### Analyzing a collection of tree files

Running Galax as follows `galax --listfile treefilelist.txt` results in a majority rule tree being saved in the file *output-galax-majrule.tre* and a summary of phylogenetic information in the file *output-galax.txt*. The majority rule tree file contains information about each node that can be viewed in FigTree (http://tree.bio.ed.ac.uk/software/figtree/) by checking "Node Labels". For each clade seen in any tree in the file *trees.t*, Galax outputs a clade descriptor (sequence of hyphens and asterisks indicating which taxa are included in the clade) and values for each of the following quantities:


## Version history

| Version  |    Released   | Description |
|:-------: | :-----------: | :---------------------------------------------------- |
|1.0       |  22-Dec-2014  | Estimates and partitions I across clades |
|1.1       |   2-Jan-2015  | Added --listfile and partitioning of D across clades |
