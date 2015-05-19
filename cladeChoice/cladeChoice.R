## R Code to allow you to choose sets of chunks for distance scores
## so you can look at the inner clade structure.

## to avoid repitition, the parts of topWords that don't need change have
## been converted into two functions, one to get the input data so you
## can cut it, and the other to take the custom cut object to get the
## distance scores for it.  dist.Source contains these functions, and
## should be in your working dir.


##source file w/ function definitions
source("distSource.R")


## set the variables that would be input parameters for topWords.
## modify as needed.
input.file <- "/home/raffle/Documents/Lexomics/TSV/merge_transpose_Bedeat1000.tsv"
##number of clades to split data into
n.clades <- 2 ##irrelevent to this, but topWords needs in
##number of top level of scores to view
n.dist.values <- 5
##options for clustering
metric <- "euclidean"
p <- 2
linkage <- "ave"
##is input file a tsv?
is.tab.sep <- TRUE
##options for output directory, chart headings, etc
dataset.name <- "Bede-AD"
##which distance metric (1 = KW, 2 = BetaBin Anova-F, 3 = Absolute)
d.metric <- 3
##Options for BetaBin.
run.batabin.par <- TRUE
num.cores <- 2
#store current working directory
start.dir <- getwd()
#create output directory name
out.dir <- paste("./", dataset.name, sep="")
print(out.dir)
#if out.dir doesn't exist, create it (it shouldn't exist, but just to be safe)
if(!file.exists(out.dir)){
    dir.create(out.dir)
}
#set working directory to output directory
setwd(out.dir)


##process the input & get gen.data
gen.data <- processInput(is.tab.sep, input.file, metric, p, linkage)

#### cut gen.data into custom groups
##get number of rows
n.rows <- dim(gen.data)[1]
##print rownames w/ index to choose from
for(i in 1:n.rows)
{
    print(paste(i, rownames(gen.data)[i], sep=": "))
}
##create an empty vector the size of the number of chunks
cut <- numeric(n.rows)
##fill w/ 1's so we only need to change the inner clade
cut[1:n.rows] <- 1
##give rownames to cut
names(cut) <- rownames(gen.data)
##store a vector with the indices of the rownames to be used in
##the second clade.
c2.inds <- c(40, 72, 73)
c2.names <- rownames(gen.data)[c2.inds]
cut[c2.names] <- 2

####If you want more than two clades, just copy the three preceeding
####lines and change values from 2 to 3, 4, or whatever.

runDist(cut, gen.data, n.clades, n.dist.values, dataset.name, run.batabin.par, num.cores)

#set working directory back to starting directory
setwd(start.dir)






