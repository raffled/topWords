## topWords
Most of the files to be run by topWords are going to be TSV files
from divitext, so it's optimized for that.

For a very basic example, assume that the input file is in your
working directory (the directory you're running R from) and you are
leaving all the parameters as their defaults.  The only three
parameters you NEED to specify are the input file, number of clades to
cut the dendrogram into, and number of top distance values you want
(this will default to a Kruskal-Wallis test), e.g.:

> topWords(input.file = "DAN450_AZ.tsv", n.clades = 2, n.dist.values  = 1)

If you have the parameters in that order, it can be even shorter:

> topWords("DAN450_AZ.tsv", 2, 1)

Both of these will create an output directory with the default name
"dist-output" in your working directory with the summary text files,
color coded dendrogram, and graphs of the relative frequencies of the
most distinct words across the top 2 clades.

Obviously you may want to an input file from another directory, in
which case you just append the pathname to the file name, e.g.:

> topWords("/home/documents/DAN450_AZ.tsv", 2, 1)

You may specify the output folder's name, but it will always be in the
working directory [which you can change in R via the command:
setwd("directory")].  The data set name is also used in the titles of
the graphs and in the text files.

>topWords("DAN450_AZ.tsv", 2, 1, dataset.name="DAZ_KW_2C_1L")

The above line will cause all the output to be in a folder named
"DAZ_KW_2C_1L" with that string as part of the graph titles & summary
files.  This is obviously useful for keeping track of what files and
folders are part of what analysis.

If you want to create the dendro with different metrics and linkages,
the parameters are "metric", "p", and "linkage", which take the same
values as dist() and hclust().

The program defaults to treating input as TSV files, but you may wish
to use any other whitespace delimited files.  For a text file (from
the Lexomics command line Perl scripts, for example) you could call
the function as:

> topWords("DAZ.txt", 2, 1, is.tab.sep=FALSE, dataset.name="DAZ")

topWords only works when the data is in a format where the chunks are
rows and the words are the columns.  Because of this the function
automatically detects the dimensions of the data being read in and if
there are more rows than columns, tranposes it.  This only works if
there are more words than chunks, but except in a few boundary cases
this will almost always be the case.  If you need to run a file with
more chunks than words, you will need to edit the code.

The above examples have all assumed the default distance metric,
Kruskal-Wallis.  As of now, there are two other options, Beta-Binomial
Anova-F and Absolute Distance.  Distance metrics are specified using a
numeric code 1-3.  1 is Kruskal-Wallis, 2 is Beta-Binomial, and 3 is
the Absoulte Distance.  The metric used is also included in the graphs
and summary files to make referencing them easier.  Because the
Beta-Binomial will require some explaining, we'll start with the
Absolute Distance.

To run the function using the Absolute Distance metric:

> topWords("DAN450_AZ.tsv", 2, 1, dataset.name="DAZ_AD_2C_1L", d.metric=3).

The output structure will still be the same as if it had been
Kruskal-Wallis, but with the values and words calculated with the
different metric.

One caveaut of the Absolte-Distance metric is words that appear in one
clade but none of the others cause a divide-by-zero error.  To account
for this, the code will give you a warning when this happens and
assign another score to these words.  The score used is the next highest+1.

The Beta-Binomial Anova-F statistic is computationally intensive (much
moreso than the other metrics), so topWords has functionality to make
use of parallel processing to come up with faster results (if you
have the "multicore" package installed).  The parameters are defaulted
for sequential analysis, so to run it sequentially:

> topWords("DAN450_AZ.tsv", 2, 1, dataset.name="DAZ_BB_2C_1L", d.metric=2)

To run the Beta-Binomial Anova-F statistic in parallel, the
run.betabin.par flag needs to be set to TRUE.  The function will run
on two cores by default when running parallel, but you can run it
across as many cores as your system has.  For example:

> topWords("DAN450_AZ.tsv", 2, 1, dataset.name="DAZ_BB_2C_1L", d.metric=2, num.cores=12)



