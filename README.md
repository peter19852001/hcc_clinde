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
Sample Running Time
====================================================================

To give an idea of the running time required by hcc_clinde and clinde,
we provide a sample running time of running hcc_clinde (on incomplete
data), clinde on both complete and incomplete data, on a synthetic GRN
with n=100 genes, and 800 time points, for 40 random replicates. The
sample data for timing is also contained in
"sample_data_for_timing.zip". The sample is run on a Linux machine
with an Intel(R) Core(TM) i7-3770 CPU @ 3.40GHz CPU with 8 cores (but
only one core is used) with 32GB ram (only a small fraction of which
is used).

The overall time is:
    real    0m35.561s
    user    0m33.650s
    sys     0m0.268s

The timings of running each replicate are as follows, where
"data=hidden" is hcc_clinde with incomplete data, "data=hiddenCL" is
clinde with incomplete data, and "data=complete" is clinde with
complete data.

n=100, e=2, alpha=1.0, nps=800, st=2, data=hidden, r=1
0.28user 0.00system 0:00.29elapsed 98%CPU (0avgtext+0avgdata 12928maxresident)k
0inputs+584outputs (0major+1410minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hidden, r=2
0.35user 0.00system 0:00.37elapsed 94%CPU (0avgtext+0avgdata 12944maxresident)k
1552inputs+632outputs (0major+1328minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hidden, r=3
0.30user 0.00system 0:00.31elapsed 96%CPU (0avgtext+0avgdata 12912maxresident)k
1552inputs+688outputs (0major+1464minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hidden, r=4
0.32user 0.00system 0:00.33elapsed 95%CPU (0avgtext+0avgdata 12912maxresident)k
1552inputs+568outputs (0major+1334minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hidden, r=5
0.32user 0.00system 0:00.32elapsed 97%CPU (0avgtext+0avgdata 12928maxresident)k
1552inputs+584outputs (0major+1341minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hidden, r=6
0.32user 0.00system 0:00.33elapsed 97%CPU (0avgtext+0avgdata 12912maxresident)k
1552inputs+624outputs (0major+1298minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hidden, r=7
0.27user 0.00system 0:00.28elapsed 97%CPU (0avgtext+0avgdata 12912maxresident)k
1552inputs+608outputs (0major+1320minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hidden, r=8
0.37user 0.00system 0:00.37elapsed 98%CPU (0avgtext+0avgdata 12944maxresident)k
1552inputs+608outputs (0major+1438minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hidden, r=9
0.29user 0.00system 0:00.30elapsed 96%CPU (0avgtext+0avgdata 12928maxresident)k
1560inputs+480outputs (0major+1306minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hidden, r=10
0.29user 0.00system 0:00.30elapsed 96%CPU (0avgtext+0avgdata 12944maxresident)k
1552inputs+584outputs (0major+1407minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hidden, r=11
0.29user 0.00system 0:00.29elapsed 97%CPU (0avgtext+0avgdata 12928maxresident)k
1552inputs+520outputs (0major+1256minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hidden, r=12
0.26user 0.00system 0:00.27elapsed 98%CPU (0avgtext+0avgdata 12944maxresident)k
1552inputs+560outputs (0major+1334minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hidden, r=13
0.35user 0.00system 0:00.35elapsed 98%CPU (0avgtext+0avgdata 12928maxresident)k
1560inputs+520outputs (0major+1392minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hidden, r=14
0.32user 0.00system 0:00.33elapsed 96%CPU (0avgtext+0avgdata 12928maxresident)k
1552inputs+624outputs (0major+1418minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hidden, r=15
0.33user 0.00system 0:00.33elapsed 98%CPU (0avgtext+0avgdata 13072maxresident)k
1552inputs+568outputs (0major+1444minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hidden, r=16
0.31user 0.00system 0:00.38elapsed 82%CPU (0avgtext+0avgdata 12944maxresident)k
1552inputs+648outputs (0major+1453minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hidden, r=17
0.28user 0.00system 0:00.29elapsed 97%CPU (0avgtext+0avgdata 12928maxresident)k
1552inputs+632outputs (0major+1423minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hidden, r=18
0.28user 0.00system 0:00.28elapsed 98%CPU (0avgtext+0avgdata 12912maxresident)k
1552inputs+528outputs (0major+1326minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hidden, r=19
0.27user 0.00system 0:00.27elapsed 98%CPU (0avgtext+0avgdata 12912maxresident)k
1552inputs+520outputs (0major+1252minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hidden, r=20
0.28user 0.00system 0:00.29elapsed 98%CPU (0avgtext+0avgdata 13056maxresident)k
1552inputs+624outputs (0major+1345minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hidden, r=21
0.29user 0.00system 0:00.30elapsed 97%CPU (0avgtext+0avgdata 12944maxresident)k
1552inputs+472outputs (0major+1248minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hidden, r=22
0.30user 0.00system 0:00.30elapsed 97%CPU (0avgtext+0avgdata 12944maxresident)k
1552inputs+520outputs (0major+1328minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hidden, r=23
0.35user 0.00system 0:00.36elapsed 98%CPU (0avgtext+0avgdata 13072maxresident)k
1552inputs+624outputs (0major+1288minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hidden, r=24
0.32user 0.00system 0:00.37elapsed 88%CPU (0avgtext+0avgdata 12944maxresident)k
1552inputs+608outputs (0major+1420minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hidden, r=25
0.31user 0.00system 0:00.33elapsed 94%CPU (0avgtext+0avgdata 12944maxresident)k
1552inputs+680outputs (0major+1405minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hidden, r=26
0.32user 0.00system 0:00.33elapsed 96%CPU (0avgtext+0avgdata 13056maxresident)k
1552inputs+536outputs (0major+1322minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hidden, r=27
0.31user 0.00system 0:00.32elapsed 98%CPU (0avgtext+0avgdata 12912maxresident)k
1552inputs+616outputs (0major+1309minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hidden, r=28
0.32user 0.00system 0:00.32elapsed 98%CPU (0avgtext+0avgdata 12912maxresident)k
1552inputs+616outputs (0major+1283minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hidden, r=29
0.34user 0.00system 0:00.35elapsed 97%CPU (0avgtext+0avgdata 12944maxresident)k
1552inputs+464outputs (0major+1249minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hidden, r=30
0.28user 0.00system 0:00.29elapsed 96%CPU (0avgtext+0avgdata 12944maxresident)k
1560inputs+480outputs (0major+1280minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hidden, r=31
0.32user 0.00system 0:00.32elapsed 98%CPU (0avgtext+0avgdata 12928maxresident)k
1552inputs+640outputs (0major+1418minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hidden, r=32
0.27user 0.00system 0:00.30elapsed 93%CPU (0avgtext+0avgdata 13072maxresident)k
1552inputs+608outputs (0major+1323minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hidden, r=33
0.29user 0.00system 0:00.30elapsed 97%CPU (0avgtext+0avgdata 12928maxresident)k
1552inputs+560outputs (0major+1311minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hidden, r=34
0.30user 0.00system 0:00.31elapsed 97%CPU (0avgtext+0avgdata 12912maxresident)k
1552inputs+672outputs (0major+1404minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hidden, r=35
0.30user 0.00system 0:00.31elapsed 97%CPU (0avgtext+0avgdata 13072maxresident)k
1552inputs+496outputs (0major+1323minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hidden, r=36
0.30user 0.00system 0:00.30elapsed 98%CPU (0avgtext+0avgdata 13072maxresident)k
1552inputs+616outputs (0major+1329minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hidden, r=37
0.38user 0.00system 0:00.40elapsed 97%CPU (0avgtext+0avgdata 13056maxresident)k
1552inputs+536outputs (0major+1357minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hidden, r=38
0.25user 0.00system 0:00.32elapsed 78%CPU (0avgtext+0avgdata 12944maxresident)k
1552inputs+568outputs (0major+1300minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hidden, r=39
0.31user 0.00system 0:00.32elapsed 96%CPU (0avgtext+0avgdata 12944maxresident)k
1552inputs+648outputs (0major+1431minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hidden, r=40
0.30user 0.00system 0:00.31elapsed 98%CPU (0avgtext+0avgdata 12912maxresident)k
1552inputs+592outputs (0major+1298minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=complete, r=1
0.25user 0.00system 0:00.25elapsed 99%CPU (0avgtext+0avgdata 14896maxresident)k
0inputs+168outputs (0major+1302minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=complete, r=2
0.28user 0.00system 0:00.32elapsed 88%CPU (0avgtext+0avgdata 14896maxresident)k
1712inputs+176outputs (0major+1302minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=complete, r=3
0.26user 0.00system 0:00.27elapsed 93%CPU (0avgtext+0avgdata 14896maxresident)k
1712inputs+184outputs (0major+1302minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=complete, r=4
0.24user 0.00system 0:00.26elapsed 96%CPU (0avgtext+0avgdata 14896maxresident)k
1712inputs+168outputs (0major+1302minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=complete, r=5
0.28user 0.00system 0:00.29elapsed 95%CPU (0avgtext+0avgdata 15024maxresident)k
1712inputs+176outputs (0major+1310minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=complete, r=6
0.26user 0.00system 0:00.27elapsed 96%CPU (0avgtext+0avgdata 14896maxresident)k
1712inputs+176outputs (0major+1302minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=complete, r=7
0.25user 0.00system 0:00.26elapsed 97%CPU (0avgtext+0avgdata 15024maxresident)k
1712inputs+168outputs (0major+1311minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=complete, r=8
0.28user 0.00system 0:00.38elapsed 73%CPU (0avgtext+0avgdata 14880maxresident)k
1712inputs+184outputs (0major+1301minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=complete, r=9
0.28user 0.00system 0:00.29elapsed 96%CPU (0avgtext+0avgdata 15008maxresident)k
1712inputs+176outputs (0major+1309minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=complete, r=10
0.25user 0.00system 0:00.26elapsed 95%CPU (0avgtext+0avgdata 14864maxresident)k
1712inputs+176outputs (0major+1300minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=complete, r=11
0.26user 0.00system 0:00.27elapsed 96%CPU (0avgtext+0avgdata 14864maxresident)k
1712inputs+176outputs (0major+1300minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=complete, r=12
0.22user 0.00system 0:00.23elapsed 96%CPU (0avgtext+0avgdata 14880maxresident)k
1712inputs+152outputs (0major+1301minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=complete, r=13
0.29user 0.00system 0:00.30elapsed 97%CPU (0avgtext+0avgdata 15024maxresident)k
1712inputs+168outputs (0major+1310minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=complete, r=14
0.28user 0.00system 0:00.28elapsed 99%CPU (0avgtext+0avgdata 14896maxresident)k
1712inputs+168outputs (0major+1302minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=complete, r=15
0.24user 0.00system 0:00.25elapsed 96%CPU (0avgtext+0avgdata 15024maxresident)k
1712inputs+176outputs (0major+1310minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=complete, r=16
0.24user 0.00system 0:00.24elapsed 97%CPU (0avgtext+0avgdata 14880maxresident)k
1712inputs+176outputs (0major+1301minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=complete, r=17
0.24user 0.00system 0:00.25elapsed 97%CPU (0avgtext+0avgdata 14896maxresident)k
1712inputs+176outputs (0major+1311minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=complete, r=18
0.26user 0.00system 0:00.26elapsed 96%CPU (0avgtext+0avgdata 14896maxresident)k
1712inputs+168outputs (0major+1303minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=complete, r=19
0.24user 0.00system 0:00.25elapsed 98%CPU (0avgtext+0avgdata 14896maxresident)k
1712inputs+168outputs (0major+1302minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=complete, r=20
0.22user 0.00system 0:00.23elapsed 97%CPU (0avgtext+0avgdata 15024maxresident)k
1712inputs+168outputs (0major+1310minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=complete, r=21
0.28user 0.00system 0:00.29elapsed 96%CPU (0avgtext+0avgdata 14880maxresident)k
1712inputs+168outputs (0major+1301minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=complete, r=22
0.24user 0.00system 0:00.24elapsed 97%CPU (0avgtext+0avgdata 14896maxresident)k
1712inputs+160outputs (0major+1303minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=complete, r=23
0.28user 0.00system 0:00.30elapsed 95%CPU (0avgtext+0avgdata 15024maxresident)k
1712inputs+192outputs (0major+1311minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=complete, r=24
0.29user 0.00system 0:00.30elapsed 97%CPU (0avgtext+0avgdata 14896maxresident)k
1712inputs+184outputs (0major+1303minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=complete, r=25
0.26user 0.00system 0:00.27elapsed 96%CPU (0avgtext+0avgdata 14864maxresident)k
1712inputs+184outputs (0major+1310minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=complete, r=26
0.25user 0.00system 0:00.34elapsed 74%CPU (0avgtext+0avgdata 15024maxresident)k
1712inputs+184outputs (0major+1310minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=complete, r=27
0.25user 0.00system 0:00.25elapsed 98%CPU (0avgtext+0avgdata 14880maxresident)k
1712inputs+168outputs (0major+1301minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=complete, r=28
0.26user 0.00system 0:00.27elapsed 97%CPU (0avgtext+0avgdata 14864maxresident)k
1712inputs+168outputs (0major+1309minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=complete, r=29
0.24user 0.00system 0:00.26elapsed 95%CPU (0avgtext+0avgdata 14896maxresident)k
1712inputs+168outputs (0major+1302minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=complete, r=30
0.25user 0.00system 0:00.26elapsed 96%CPU (0avgtext+0avgdata 14864maxresident)k
1712inputs+176outputs (0major+1300minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=complete, r=31
0.26user 0.00system 0:00.28elapsed 98%CPU (0avgtext+0avgdata 14848maxresident)k
1704inputs+184outputs (0major+1299minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=complete, r=32
0.24user 0.00system 0:00.25elapsed 97%CPU (0avgtext+0avgdata 15008maxresident)k
1712inputs+168outputs (0major+1309minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=complete, r=33
0.23user 0.00system 0:00.25elapsed 94%CPU (0avgtext+0avgdata 14880maxresident)k
1712inputs+160outputs (0major+1301minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=complete, r=34
0.24user 0.00system 0:00.25elapsed 95%CPU (0avgtext+0avgdata 14864maxresident)k
1712inputs+176outputs (0major+1301minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=complete, r=35
0.26user 0.00system 0:00.27elapsed 96%CPU (0avgtext+0avgdata 15024maxresident)k
1712inputs+184outputs (0major+1311minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=complete, r=36
0.23user 0.00system 0:00.24elapsed 97%CPU (0avgtext+0avgdata 15008maxresident)k
1704inputs+168outputs (0major+1309minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=complete, r=37
0.34user 0.00system 0:00.35elapsed 98%CPU (0avgtext+0avgdata 15024maxresident)k
1712inputs+184outputs (0major+1310minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=complete, r=38
0.22user 0.00system 0:00.24elapsed 95%CPU (0avgtext+0avgdata 14896maxresident)k
1712inputs+168outputs (0major+1303minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=complete, r=39
0.24user 0.00system 0:00.25elapsed 96%CPU (0avgtext+0avgdata 14864maxresident)k
1712inputs+184outputs (0major+1301minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=complete, r=40
0.26user 0.00system 0:00.28elapsed 96%CPU (0avgtext+0avgdata 14896maxresident)k
1712inputs+192outputs (0major+1302minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hiddenCL, r=1
0.24user 0.00system 0:00.24elapsed 98%CPU (0avgtext+0avgdata 12880maxresident)k
0inputs+144outputs (0major+1145minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hiddenCL, r=2
0.30user 0.00system 0:00.31elapsed 99%CPU (0avgtext+0avgdata 12896maxresident)k
0inputs+160outputs (0major+1154minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hiddenCL, r=3
0.26user 0.00system 0:00.26elapsed 98%CPU (0avgtext+0avgdata 12880maxresident)k
0inputs+152outputs (0major+1145minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hiddenCL, r=4
0.28user 0.00system 0:00.28elapsed 99%CPU (0avgtext+0avgdata 12912maxresident)k
0inputs+144outputs (0major+1147minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hiddenCL, r=5
0.27user 0.00system 0:00.27elapsed 99%CPU (0avgtext+0avgdata 12896maxresident)k
0inputs+144outputs (0major+1146minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hiddenCL, r=6
0.28user 0.00system 0:00.28elapsed 99%CPU (0avgtext+0avgdata 12912maxresident)k
0inputs+152outputs (0major+1147minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hiddenCL, r=7
0.23user 0.00system 0:00.24elapsed 99%CPU (0avgtext+0avgdata 12896maxresident)k
0inputs+144outputs (0major+1146minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hiddenCL, r=8
0.30user 0.00system 0:00.30elapsed 98%CPU (0avgtext+0avgdata 12896maxresident)k
0inputs+160outputs (0major+1146minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hiddenCL, r=9
0.26user 0.00system 0:00.26elapsed 99%CPU (0avgtext+0avgdata 12864maxresident)k
0inputs+152outputs (0major+1144minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hiddenCL, r=10
0.24user 0.00system 0:00.24elapsed 99%CPU (0avgtext+0avgdata 12912maxresident)k
0inputs+152outputs (0major+1148minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hiddenCL, r=11
0.26user 0.00system 0:00.26elapsed 98%CPU (0avgtext+0avgdata 12880maxresident)k
0inputs+152outputs (0major+1146minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hiddenCL, r=12
0.21user 0.00system 0:00.22elapsed 98%CPU (0avgtext+0avgdata 12912maxresident)k
0inputs+128outputs (0major+1148minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hiddenCL, r=13
0.30user 0.00system 0:00.30elapsed 99%CPU (0avgtext+0avgdata 12880maxresident)k
0inputs+144outputs (0major+1145minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hiddenCL, r=14
0.28user 0.00system 0:00.28elapsed 99%CPU (0avgtext+0avgdata 12896maxresident)k
0inputs+144outputs (0major+1146minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hiddenCL, r=15
0.27user 0.00system 0:00.27elapsed 98%CPU (0avgtext+0avgdata 13040maxresident)k
0inputs+160outputs (0major+1156minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hiddenCL, r=16
0.26user 0.00system 0:00.27elapsed 98%CPU (0avgtext+0avgdata 12880maxresident)k
0inputs+152outputs (0major+1145minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hiddenCL, r=17
0.23user 0.00system 0:00.23elapsed 99%CPU (0avgtext+0avgdata 12880maxresident)k
0inputs+152outputs (0major+1145minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hiddenCL, r=18
0.24user 0.00system 0:00.24elapsed 99%CPU (0avgtext+0avgdata 12880maxresident)k
0inputs+136outputs (0major+1145minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hiddenCL, r=19
0.24user 0.00system 0:00.24elapsed 98%CPU (0avgtext+0avgdata 12880maxresident)k
0inputs+144outputs (0major+1145minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hiddenCL, r=20
0.24user 0.00system 0:00.25elapsed 98%CPU (0avgtext+0avgdata 13024maxresident)k
0inputs+144outputs (0major+1154minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hiddenCL, r=21
0.26user 0.00system 0:00.27elapsed 98%CPU (0avgtext+0avgdata 12896maxresident)k
0inputs+144outputs (0major+1146minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hiddenCL, r=22
0.26user 0.00system 0:00.26elapsed 98%CPU (0avgtext+0avgdata 12912maxresident)k
0inputs+136outputs (0major+1148minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hiddenCL, r=23
0.31user 0.00system 0:00.31elapsed 99%CPU (0avgtext+0avgdata 13024maxresident)k
0inputs+160outputs (0major+1154minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hiddenCL, r=24
0.28user 0.00system 0:00.28elapsed 99%CPU (0avgtext+0avgdata 12896maxresident)k
0inputs+160outputs (0major+1154minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hiddenCL, r=25
0.27user 0.00system 0:00.27elapsed 99%CPU (0avgtext+0avgdata 12912maxresident)k
0inputs+152outputs (0major+1148minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hiddenCL, r=26
0.28user 0.00system 0:00.28elapsed 98%CPU (0avgtext+0avgdata 13024maxresident)k
0inputs+152outputs (0major+1154minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hiddenCL, r=27
0.27user 0.00system 0:00.27elapsed 99%CPU (0avgtext+0avgdata 12896maxresident)k
0inputs+144outputs (0major+1156minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hiddenCL, r=28
0.28user 0.00system 0:00.28elapsed 98%CPU (0avgtext+0avgdata 12912maxresident)k
0inputs+144outputs (0major+1155minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hiddenCL, r=29
0.31user 0.00system 0:00.31elapsed 99%CPU (0avgtext+0avgdata 12880maxresident)k
0inputs+144outputs (0major+1154minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hiddenCL, r=30
0.24user 0.00system 0:00.25elapsed 98%CPU (0avgtext+0avgdata 12880maxresident)k
0inputs+144outputs (0major+1146minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hiddenCL, r=31
0.25user 0.00system 0:00.25elapsed 99%CPU (0avgtext+0avgdata 12896maxresident)k
0inputs+160outputs (0major+1148minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hiddenCL, r=32
0.24user 0.00system 0:00.24elapsed 99%CPU (0avgtext+0avgdata 13040maxresident)k
0inputs+144outputs (0major+1155minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hiddenCL, r=33
0.24user 0.00system 0:00.25elapsed 99%CPU (0avgtext+0avgdata 12896maxresident)k
0inputs+136outputs (0major+1146minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hiddenCL, r=34
0.25user 0.00system 0:00.25elapsed 99%CPU (0avgtext+0avgdata 12912maxresident)k
0inputs+152outputs (0major+1147minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hiddenCL, r=35
0.27user 0.00system 0:00.27elapsed 99%CPU (0avgtext+0avgdata 13024maxresident)k
0inputs+152outputs (0major+1154minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hiddenCL, r=36
0.25user 0.00system 0:00.26elapsed 99%CPU (0avgtext+0avgdata 13040maxresident)k
0inputs+144outputs (0major+1155minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hiddenCL, r=37
0.34user 0.00system 0:00.35elapsed 99%CPU (0avgtext+0avgdata 13024maxresident)k
0inputs+152outputs (0major+1154minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hiddenCL, r=38
0.22user 0.00system 0:00.22elapsed 99%CPU (0avgtext+0avgdata 12912maxresident)k
0inputs+136outputs (0major+1147minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hiddenCL, r=39
0.26user 0.00system 0:00.26elapsed 99%CPU (0avgtext+0avgdata 12864maxresident)k
0inputs+152outputs (0major+1144minor)pagefaults 0swaps
n=100, e=2, alpha=1.0, nps=800, st=2, data=hiddenCL, r=40
0.26user 0.00system 0:00.26elapsed 99%CPU (0avgtext+0avgdata 12896maxresident)k
0inputs+160outputs (0major+1146minor)pagefaults 0swaps

====================================================================
