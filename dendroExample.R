## These are a few lines of code to get Pat started making dendros

library(stats)     ##load the stats library

##set the input file
input.file <- "the path/file.tsv"

##read it in

test.data <- read.table(as.character(input.file), header=T, comment.char="", row.names=1, sep="\t", quote="", fill = TRUE)


#move test.data to gen.data. Check to see based on dimensions if
#transposition is necessary
if(dim(test.data)[1] < dim(test.data)[2]){
    gen.data <- test.data
}
else{
    gen.data <- t(test.data)
}

# Get the row sums for the gen.data
rsums <- apply(gen.data,1,sum)

######test to see if scrubber & divitext options are included in
##first column of bottom two rows of .tsv (added at divitext
##through php) to facilitate subtitle on dendro listing all
##options used when generating data
##find number of rows
n.rows <- dim(gen.data)[1]
##if the last one or two rows have NA values for word counts, recreate
##gen.data without those rows
if(sum(is.na(rsums)[(n.rows-1):n.rows]) == 2){
    gen.data <- gen.data[-c(n.rows-1,n.rows),]
    rsums <- rsums[-c(n.rows-1,n.rows)]
}else if(sum(is.na(rsums)[(n.rows-1):n.rows]) == 1){
    gen.data <- gen.data[-c(n.rows),]
    rsums <- rsums[-c(n.rows)]
}

##store counts for betabin processing
gen.counts <- gen.data

####Normalize gen.data by converting word counts to relative frequencies
# create a matrix that gives the row sums as denominators at each
# element of the data matrix to allow conversion to relative frequencies
denoms <- matrix(rep(rsums,dim(gen.data)[2]),
                 byrow=F,ncol=dim(gen.data)[2])
# divide each member of gen.data by it's computed denominator to get
# relative frequency
gen.data <- gen.data/denoms


# compute distances between all vectors, using input parameters
dist.gen <- dist(gen.data, method="euclidean", p=2)
#generate cluster object of dist.gen using input parameter
hc <- hclust(dist.gen, method="ave")

##plot w/ just the default graphics output
plot(hc)

##or do it w/ a pdf (or any image type)
pdf(file="output.pdf")
plot(hc)
