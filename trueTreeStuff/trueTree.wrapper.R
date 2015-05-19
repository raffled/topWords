## trueTree.wrapper.R

## This program provides an easy-to-use wrapper for running trueTree
## locally.  Variables that are listed here are the most common ones
## you may need to run trueTree for a given dataset on one machine.
## For more advanced options, see the documentation for trueTree itself.

## To run trueTree on a file, set the parameters and run the code
## below in R.  Basic descriptions of the parameters are provided.
## For more details, see the original documentation in trueTree.R


##########################################################################
############## - UPDATE AND RUN ALL OF THE FOLLOWING CODE - ##############
##########################################################################

#################### - Get & Source TrueTree Files - #####################

##this is the directory the trueTree files are in.
trueTree.dir <- "/home/raffle/Documents/School_Archive/Lexomics/Wheaton/trueTreeStuff/Lexomics-Cluster-Validation"
setwd(trueTree.dir)
source("trueTree.r")

######################## - Input Parameters - ############################

##If the input file is in the trueTree directory, a path is not
##needed.  If they are in separate directories, provide the full path.
input.file <- "/directoryPath/inputFilename"

##### Output (omit extension from output file, it's added automatically)####
#use this if you just want the output in the working (trueTree) directory
outputFilename <- "outputFile"

#or use this if you want to specify the output directory
output.dir <- "/outputPath/"
if(!file.exists(output.dir)){
    dir.create(output.dir)
}
output.file.name <- "name of the output file"
outputFileName <- paste(output.dir, output.file.name, sep="")
#### end output section ####

##If you don't want to title the dendro, set main to NULL
main <- "Graph Title"

##Options for creating the dendro
distMetric <- "euclidean"
clustMethod <- "average"

##If your input file has chunks as rows, set this to false.  Otherwise
##keep true.  TSV files from divitext should leave this as TRUE
input.transposed <- TRUE

##number of bootstrap aggregations to perform.  10,000 has the best
##precision/performance ratio, but may take a while to compute.
nboot <- 10000

##Set this to TRUE if you want to run the program over multiple cores,
##or leave as FALSE to run sequentially
runParallel <- FALSE

##This is the number of cores to use to compute the bootstrap.  If
##you're running sequentially you can just ignore this.
numCPUs <- 2

##height and width of dendro (the default is too narrow to print with
##the au values readable.
height <- 900
width <- 2500

##################### - call the trueTree function - #####################
##set the result to a variable in case you need to re-plot it
result <- trueTree(input.file=input.file,
                   outputFileName=outputFileName, main=main,
                   distMetric=distMetric, clustMethod=clustMethod,
                   input.transposed=input.transposed, nboot=nboot,
                   runParallel=runParallel, numCPUs=numCPUSs,
                   height=height, width=width)
##code for re-plot if you need it.
plot.trueTree(result, outputFilename=outputFilename, main=main, height=height, width=width)

##########################################################################
############################ - Example  - ################################
trueTree.dir <- "/home/raffle/work/raffle_work/trueTreeStuff/Lexomics-Cluster-Validation"
setwd(trueTree.dir)
source("trueTree.r")

input.file <- "/home/raffle/Documents/School_Archive/Lexomics/TSV/merge_transpose_Bedeat1000.tsv"
#input.file <- "/home/raffle/work/TSVs/DAZ/DAN450_AZ.tsv"
output.dir <- "/home/raffle/Documents/School_Archive/Lexomics/TSV/"
if(!file.exists(output.dir)){
    dir.create(output.dir)
}
output.file.name <- "bede"
outputFilename <- paste(output.dir, output.file.name, sep="")


main <- "Bede"
distMetric <- "euclidean"
clustMethod <- "average"
input.transposed <- TRUE
nboot <- 100
runParallel <- TRUE
numCPUs <- 4
height <- 800
width <- 800

result <- trueTree(input.file=input.file,
                   outputFilename=outputFilename, main=main,
                   distMetric=distMetric, clustMethod=clustMethod,
                   input.transposed=input.transposed, nboot=nboot,
                   runParallel=runParallel, numCPUs=numCPUs,
                   height=height, width=width)

##Reset window size & re-plot
height <- 1000
width <- 1200
plot.trueTree(result, outputFilename=outputFilename, main=main, height=height, width=width)
