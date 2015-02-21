# hcc_clinde

A tool to infer causal regulatory network with delays from continuous
time series data, where the expression of a small but unknown number
of nodes is hidden.

====================================================================
Introduction to hcc_clinde
====================================================================

hcc_clinde is an algorithm to recover the causal GRN with delays in
the regulatory links from time series expression data, where a small
but unknown number of nodes are hidden, i.e. without expression data.

We assume that there is a sparse causal graph of of interest, where
the variables are continuous and each causal link has a delay
(possibly more than one time step). A small but unknown number of
variables are not observed. Each hidden variable has only observed
variables as children and parents, where the number of children is at
least two, and the number of parents is non-negative. Our objective is
to infer the causal graph with the delays, given the time series of
the observed variables. Since it is difficult to obtain very long time
series, hcc_clinde is also capable of utilizing multiple short time
series, which are not necessarily of the same length (e.g. obtained
from replicate experiments). The programs have been tested in Debian
and Ubuntu Linux distributions, and should build and work with minimal
tweaking in other Linux distributions.

==== Building gsl

Before buliding hcc_clinde and CLINDE/clinde, you need to build gsl,
which is the GNU Scientific Library, which contains routines used by
hcc_clinde and CLINDE.

gsl-1.16.tar.gz is included in CLINDE/ for convenience, and we have tested using
gsl-1.16. First extract it by:

    tar xzvf gsl-1.16.tar.gz

Then to build it, typically you can use:

    cd gsl-1.16/
    ./configure
    make
    make install

But refer to gsl-1.16/README for trouble shooting.

==== Building hcc_clinde

Having built gsl, you could build hcc_clinde by using the provided makefile:

    make

You may also build it directly using:

    gcc -Wall -O2 -static hcc_clinde.c parse_option.c tsv.c -o hcc_clinde -lgsl -lgslcblas -lm

==== Usage:

The usage of hcc_clinde is:

    Usage: ./hcc_clinde [-?] -data item1 [item2 ...] [-st1 real] [-st2 real]
    [-max.delay int] [-max.n int] [-method str] [-pruning str] [-keep.zero.delay]
    [-one.delay] [-no.dup] [-normalize] [-ctL real] [-ctU real] [-eV real]
    [-eV.tolerance real]
    
    Description of the options:
      -?:  Showing the usage.
    
      -data [REQUIRED]:  File name(s) of the input expression data (tab/space
        separated), each row is a time point, each column is a gene. Each file
        should have the same number of columns, but may have different number of
        rows.
    
      -st1:  Score threshold for stage 1. For method=pcor, score is -log10(p),
        where p is the intended p-value. For method=mi, score is the mutual
        information. Default 2.
    
      -st2:  Score threshold for stage 2. Similar to st1. Default 2.
    
      -max.delay:  Maximum delay (lags) to use in the inference, default 4.
    
      -max.n:  Maximum number of parents to condition on in stage 2, default 4.
    
      -method:  Method of testing links in stage 1 and 2, can be pcor (partial
        correlation) or mi (mutual information). Default "pcor"
    
      -pruning:  Pruning strategy in stage 2, can be all (consider all neighbors
        of the two vertices when doing conditional test of the link between the
        two) or common (consider only common neighbors of the two vertices).
        Default "all"
    
      -keep.zero.delay:  To keep links with zero delay after stage 2 and before
        remove dups. Default false.
    
      -one.delay:  To keep only one delay (the one with best score, smallest
        delay) for each link after stage 1. Default false.
    
      -no.dup:  Remove duplicate links. To keep only one delay (the one with best
        score, smallest delay) for each link after stage 2. Default false.
    
      -normalize:  Center and normalize the expressions such that for each gene,
        the mean is 0, and s.d. is 1. Default false.
    
      -ctL:  Lower Correlation threshold for initial clustering. If ctL <= abs
        correlation < ctU, consider the two nodes for having hidden common node.
        Default 0.75.
    
      -ctU:  Upper Correlation threshold for initial clustering. If ctU <= abs
        correlation, assume a direct link but NOT hidden common node. Default
        0.86.
    
      -eV:  Expected Error Variance if > 0. If <= 0, the expected error variance
        will be estimated to be the median of the error variances of the initial
        GRN. If error variance of a gene is > eV*(1+eV.tolerance), that gene will
        be considered for having hidden common node. Default 1.
    
      -eV.tolerance:  If error variance of a gene is > eV*(1+eV.tolerance), that
        gene will be considered for having hidden common node. Default 0.1.
    

==== Example usage:

An example use of hcc_clinde is:

    ./hcc_clinde -eV 0 -data sample_data/r9_s*.txt -no.dup > output.txt

Now output.txt contains some messages and the GRN after various stages
of hcc_clinde. To extract the final GRN (e.g. for comparison), you may
use sed and grep, which are available in Linux:

    sed -n '/^==== Final GRN/,/^==== End Final GRN/p' output.txt | grep "^To:" > outgrn.txt

Now outgrn.txt contains the GRN in a format acceptable for the
comparison tool grn_cmp_hcc.

====================================================================
Introduction to CLINDE (for comparison purpose)
====================================================================

CLINDE/clinde is a tool to infer causal GRN with delays from time
series expression data, but does not handle hidden nodes.

==== Building

Before buliding CLINDE/clinde, you need to build gsl, which is the GNU
Scientific Library, which contains routines used by CLINDE.

gsl-1.16.tar.gz is included for convenience, and we have tested using
gsl-1.16. First extract it by:

    tar xzvf gsl-1.16.tar.gz

Then to build it, typically you can use:

    cd gsl-1.16/
    ./configure
    make
    make install

But refer to gsl-1.16/README for trouble shooting.

After that, you could build clinde by using the provided makefile in CLINDE:

    cd CLINDE/
    make

You may also build it directly using:

    cd CLINDE/
    gcc -Wall -O2 -static clinde.c parse_option.c tsv.c -o clinde -lgsl -lgslcblas -lm

==== Usage:

The usage of CLINDE/clinde is:

    Usage: ./clinde [-?] -data item1 [item2 ...] [-st1 real] [-st2 real]
    [-max.delay int] [-max.n int] [-method str] [-pruning str] [-one.delay]
    [-no.dup]
    
    Description of the options:
      -?:  Showing the usage.
    
      -data [REQUIRED]:  File name(s) of the input expression data (tab/space
        separated), each row is a time point, each column is a gene. Each file
        should have the same number of columns, but may have different number of
        rows.
    
      -st1:  Score threshold for stage 1. For method=pcor, score is -log10(p),
        where p is the intended p-value. For method=mi, score is the mutual
        information. Default 2.
    
      -st2:  Score threshold for stage 2. Similar to st1. Default 2.
    
      -max.delay:  Maximum delay (lags) to use in the inference, default 4.
    
      -max.n:  Maximum number of parents to condition on in stage 2, default 4.
    
      -method:  Method of testing links in stage 1 and 2, can be pcor (partial
        correlation) or mi (mutual information). Default "pcor"
    
      -pruning:  Pruning strategy in stage 2, can be all (consider all neighbors
        of the two vertices when doing conditional test of the link between the
        two) or common (consider only common neighbors of the two vertices).
        Default "all"
    
      -one.delay:  To keep only one delay (the one with best score, smallest
        delay) for each link after stage 1. Default false.
    
      -no.dup:  Remove duplicate links. To keep only one delay (the one with best
        score, smallest delay) for each link after stage 2. Default false.
    

==== Example usage:

An example use for inferring initial GRN is:

    CLINDE/clinde -data sample_data/r9_s*.txt -no.dup > output_CL.txt

Now output_CL.txt contains some messages and the GRN after stage 1 and
stage 2 of CLINDE. To extract the GRN after stage 2 for use in the
next step, you may use sed and grep, which are available in Linux:

    sed -n '/^==== Stage 2/,/^==== End Stage 2/p' output_CL.txt | grep "^To:" > outgrn_CL.txt

Now outgrn_CL.txt contains the GRN in a format acceptable for the
comparison tool grn_cmp_hcc.

====================================================================
Comparison of two GRNs (if applicable)
====================================================================

If the true GRN is known, you may assess how close the predicted GRN
is to the true GRN by using grn_cmp_hcc.

==== Building grn_cmp_hcc:

You may build it by using the provided makefile:

    make grn_cmp_hcc

Or you may build it directly by:

    gcc -Wall -O3 grn_cmp_hcc.c parse_option.c -o grn_cmp_hcc -lm

==== Usage of grn_cmp_hcc:

    Usage: ./grn_cmp_hcc [-?] -p str -t str [-n int] [-v]
    
    Description of the options:
      -?:  Showing the usage.
    
      -p [REQUIRED]:  File name of the predicted GRN. Each line in the file
        represents an edge, which consists of a 0-based 'from' index, a 0-based
        'to' index, a delay, and the effect, separated by space..
    
      -t [REQUIRED]:  File name of the true GRN. Each line in the file represents
        an edge, which consists of a 0-based 'from' index, a 0-based 'to' index, a
        delay, and the effect, separated by space..
    
      -n:  Number of non-hidden genes, unspecified if <= 0. If an index is >=
        this, it is regarded as a hidden node. Default 0.
    
      -v:  Verbose mode.


The parameter -n specifies the number of observed genes, which is the
number of columns in the sample data.


==== Example usage:

An example use of grn_cmp_hcc is:

    ./grn_cmp_hcc -p outgrn.txt -t sample_data/r9_grn.txt -n 50

====================================================================
R Scripts for Generation of Synthetic Data (Optional)
====================================================================

The following functions in synthetic_small_cases_e.R and
synthetic_misc_cases.R can help generate a lot of different types of
synthetic data for testing purposes. Of course R need to be installed
first.

Useful functions:

synthetic_small_cases_e.R:

    "gen.hidden.cases()": small case with one hidden node, using one long time series.

    "gen.non.hidden.cases()": small case with no hidden node, using one long time series.

    "gen.seg.hidden.cases()": small case with one hidden node, using multiple short time series.

    "gen.seg.non.hidden.cases()": small case with no hidden node, using multiple short time series.

synthetic_misc_cases.R:

    "gen.misc.cases()": large case with 5 to 10 hidden nodes, using one long time series.

    "gen.seg.misc.cases()": large case with 5 to 10 hidden nodes, using multiple short time series.

Be warned that a LOT of synthetic data will be generated, and it may
take some time (up to a few hours on a slow computer) because R is not
very fast.

You may tweak the functions to generate cases with other parameters.

====================================================================
