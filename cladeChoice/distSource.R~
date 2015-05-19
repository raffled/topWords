#load libraries needed
library(stats)
currentNode = 1 #global variable to keep track of position while
                #recursively travelling dendrogram during the plotting
                #code


processInput <- function(is.tab.sep, input.file, metric, p, linkage)
{
    # read in the data, changing comment character to nothing so # isn't confused
    # check if read.table should use default separators or tabs based on parameter
    if(is.tab.sep == TRUE){
      test.data <- read.table(input.file,
                              header=T,comment.char="",sep="\t",row.names=1, check.names=FALSE, fill = TRUE)
    }
    else{
      test.data <- read.table(input.file, header=T,comment.char="",row.names=1)
    }

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
    ##if the last one or two rows have NA values for word counts, recreate gen.data without
    ##those rows
    if(sum(is.na(rsums)[(n.rows-1):n.rows]) == 2){
        gen.data <- gen.data[-c(n.rows-1,n.rows),]
        rsums <- rsums[-c(n.rows-1,n.rows)]
    }else if(sum(is.na(rsums)[(n.rows-1):n.rows]) == 1){
        gen.data <- gen.data[-c(n.rows),]
        rsums <- rsums[-c(n.rows)]
    }

    return(gen.data)
}

    ###################################################################################################
    #Distance Test
    ###################################################################################################



runDist <- function(cut, gen.data, n.clades, n.dist.values, dataset.name, run.betabin.par, num.cores)
{
    #get order of indexes for reordering gen.data based on which chunks
    #are in which clade of the cut
    reorderrows.index <- order(cut)

    #reorder gen.data based on indexes from above.  After executing this command,
    #chunks are regrouped so that the first N1 rows are in clade1, the next N2
    #rows are in clade2, etc.
    reorder.gen.data <- gen.data[reorderrows.index,]

    #get the table of how many chunks are in each cut
    n.clade <- table(cut)

    ##store counts for betabin processing
    gen.counts <- gen.data

    rsums <- apply(gen.data,1,sum)

    ####Normalize gen.data by converting word counts to relative frequencies
    # create a matrix that gives the row sums as denominators at each element of the data matrix
    # to allow conversion to relative frequencies
    denoms <- matrix(rep(rsums,dim(gen.data)[2]),
                     byrow=F,ncol=dim(gen.data)[2])
    # divide each member of gen.data by it's computed denominator to get relative frequency
    gen.data <- gen.data/denoms


    # compute distances between all vectors, using input parameters
    dist.gen <- dist(gen.data, method=metric, p=p)

    #generate cluster object of dist.gen using input parameter
    hc <- hclust(dist.gen, method=linkage)
    #switch on d.metric to get d.metric.name
    switch(d.metric,
           d.metric.name <- "KRUSKAL-WALLIS",
           d.metric.name <- "BETABIN-ANOVA",
           d.metric.name <- "ABSOLUTE-DISTANCE"
      )


    ##print list of elements to elements.txt
    #generate file heading and write to file
    heading <- paste("Chunks in each clade after cutting into ", n.clades, " clades:",sep="")
    write(heading, "./clade-listing.txt")
    #for each of the unique clades from cut, write which chunks are in that clade
    for(i in unique(cut))
      {
        write(paste("\nElements of Clade ", i, ":", sep=""), "./clade-listing.txt", append=TRUE)
        write(names(cut[cut==i]), "./clade-listing.txt",ncolumns=3,append=TRUE, sep = ", ")
      }

    #creates strings of the clade names based on the the clade number for distance test
    fdoc.collection <- factor(cut, levels=1:n.clades)
    levels(fdoc.collection) <- paste("clade",1:n.clades,sep="")

    ##Plot the dendro with each clade color-coded
    plot.trueTree(hc, outputFilename = paste(dataset.name, "_dendro", sep=""),
                  dataset.name=dataset.name, fdoc.collection=fdoc.collection)

    #create an array to store distance stat for all of the words
    compare.dist <- vector(length = number.words <- dim(gen.data)[2])

    #switch on d.metric to decide which distance metric to use to fill
    #compare.dist
    switch(d.metric,
           compare.dist <- getKW(number.words, fdoc.collection,
                                 gen.data, compare.dist),
           compare.dist <- getBetaBinAnova(number.words, fdoc.collection,
                                    gen.counts, compare.dist, rsums, run.betabin.par, num.cores),
           compare.dist <- getSydneyScore(gen.data,fdoc.collection, n.clades)
          )
    #let's order the results (distance values) from the biggest to the smallest
    ordered.diff <- order(compare.dist, decreasing = TRUE)

    #sort & store the values of compare.dist from largest to smallest
    sorted.uniquevals.dist <- rev(sort(unique(compare.dist)))

    #stores each of the top n.dist.values stats from distance test
    topscores.dist <- unique(sorted.uniquevals.dist[1:n.dist.values])

    #create list to store words at each level of KW score
    topwords.dist <- NULL

    #loop through each level of distance scores.  Words at each level are stored as list within the
    # list of lists topwords.dist, where each of the sublists are the list of all the words
    # at each of the distance values by index, where
    # index == 1 corresponds to the words with largest distance stat, index == 2
    # corresponds to the words with next largest distance stat,  and so
    # on... Can retrieve these using
    # topwords.dist[[1]], topwords.dist[[2]], etc.
    for(i in 1:n.dist.values)
      {
        topwords.dist <- c(topwords.dist,list(colnames(gen.data)[compare.dist == sorted.uniquevals.dist[i]]))
      }

    ##write words at each level to a file
    #create and write heading
    heading <- paste("LIST OF WORDS FOR EACH ", d.metric.name, "  DISTANCE SCORE.\n")
    write(heading, "./words-list.txt")
    #loop for however many top distance scores the user specifies
    #printing all the words that occur at that level
    for(i in 1:n.dist.values)
      {
        section.heading <- paste("\n - Words at ", getOrder(i), " ", d.metric.name, " score (",
                                 round(sorted.uniquevals.dist[i],2), ") - \n", sep="")
        write(section.heading, "./words-list.txt", append=TRUE)
        cat(topwords.dist[[i]], file="./words-list.txt",append=TRUE, fill=TRUE, sep=", ")
      }

    #####create plots for each word at each level of distance score
    #define x as the list of which chunks are in which clade
    x <- rep(1:n.clades, n.clade)
    #define n as the number of chunks being plotted
    n <- length(x)

    #create graphs directory at tempDirectory/graphs/
    if(!file.exists("./graphs")){
      dir.create("./graphs")
    }
    #set working directory to graphs directory
    setwd("./graphs")

    #### Print graphs
    #outer loop loops for each level i where i corresponds to the top levels
    #of distance scores for however many levels user wants to look at

#### commented this at 2:37, 6/20
    for(i in 1:n.dist.values)
      {
        #inner loop loops through every word [j] at each level of scores [i]
        for(j in 1:length(topwords.dist[[i]]))
          {
            #open pdf device and set file name for file
            pdf(file=paste(dataset.name, "_", d.metric.name,"_", round(topscores.dist[i],2), "_",
                  topwords.dist[[i]][j], ".pdf", sep=""))
            #create plot for word
            plot(x + (.1)*runif(n, -1, 1), y <-
                 reorder.gen.data[,topwords.dist[[i]][j]],
                 xlim=c(min(x)*.8, max(x)*1.2), ylim=c(0, 1.1*max(y)),
                 #set plot character and x,y labels
                 pch="*", xlab="Cluster", ylab="Relative Frequency",
                 #give the graph title some details about the word so the graph is useful
                 main=paste(dataset.name, " | ",  getOrder(i), " Distance \n", d.metric.name, "  Score: ",
                   round(topscores.dist[i], 2), " | Word[", j, "/",
                   length(topwords.dist[[i]]), "]: ", topwords.dist[[i]][j], sep = ""),
                 axes=FALSE)
            #specify parameters for the axes so they're not full of useless info
            axis(side=1, at=(ux <- unique(x)), labels=as.character(ux))
            axis(side=2, at=(uy <- seq(min(y),max(y), length=5)), labels=signif(uy, digits=2))
            box(which = "plot")

            #close pdf device(s) so R doesn't shit itself from too many open devices
            graphics.off()
          }

      }

}

##this function takes an integer and returns a string w/ the ordinal
##for more useful printing.  Only works on 1-20, if they want more than
##that they can deal with the 21TH, etc for the time being.
getOrder <- function(n){
  order <- NULL
  #switch on n for 1-3, which are unique
  switch(n,
         order <- "LARGEST",
         order <- "2ND LARGEST",
         order <- "3RD LARGEST"
         )
  #if n is > 3, order is set to null. From these just appending TH to n
  #takes care of most cases.
  if(is.null(order)){
    order <- paste(n, "TH LARGEST", sep="")
  }
  return(order)
}

## This function takes gen.data and fdoc.collection to find the
## Kruskal-Wallis scores for gen.data.  The function returns compare.dist
getKW <- function(number.words, fdoc.collection, gen.data, compare.dist){
    #fill compare.dist with the results of KW analysis
    for(j in 1:number.words)
      {
        compare.dist[j] <-
          as.numeric(kruskal.test(gen.data[,j]~fdoc.collection)$statistic)
      }
    return(compare.dist)
}

## ## This function takes gen.data and fdoc.collection to find the
## ## Anova-F scores for gen.data.  The function returns compare.dist
## getAnova <- function(number.words, fdoc.collection, gen.data, compare.dist){
##   ##tell the user if anova is a stupid choice.
##   for(j in 1:number.words)
##     {
##       compare.dist[j] <- as.numeric(anova(lm(gen.data[,j]~as.factor(fdoc.collection)))$F[1])
##     }
##   return(compare.dist)
## }

##New anova
getBetaBinAnova <- function(number.words, fdoc.collection, gen.counts,
                            compare.dist, rsums, run.betabin.par, num.cores){

    library(aod)

    #start system time, as this can take a fuckin' long time.
    start.t <- Sys.time()
    n.clades <- length(unique(fdoc.collection))
    #check if they're doing this parallel or sequential
    if(run.betabin.par == TRUE){
        library(multicore)

        #coerce gen.counts into a list for multicore processing
        gen.list <- as.list(gen.counts, all.names=TRUE)

        #names(gen.list) <- names(compare.dist) <- colnames(gen.counts)
        #print(gen.list[1:5])
        #print(names(gen.list["the"]))
        #function for parallel processing.  Takes the counts for each
        #chunk one at a time and returns the betaBin score one at a time.
        betaBinParallel <- function(gen.word, rsums, fdoc.collection, gen.counts){

            gen.word <- data.frame(gen.word)
            L1 <- attributes(m1 <- betabin(cbind(gen.word, rsums - gen.word)
                                           ~fdoc.collection - 1, ~1, data=gen.word,
                                           control = list(maxit = 10, reltol = .001),
                                           fixpar=list(n.clades+1, 0)))$logL
            L0 <- attributes(m0 <- betabin(cbind(gen.word, rsums - gen.word)
                                           ~1, ~1, data=gen.word ,
                                           control = list(maxit = 10, reltol = .001),
                                           fixpar=list(2, 0)))$logL
            compare.word <- 2 * (L1-L0)
            gc()
            return(compare.word)
        }


        #mclapply() returns a list of all the values from
        #betaBinParallel.  So a vector is created from the unlisted values
        #returned and stored as compare.dist
        compare.dist <- as.vector(unlist(mclapply(gen.list,
                                                  betaBinParallel,
                                                  rsums,
                                                  fdoc.collection, gen.counts,
                                                  mc.cores=num.cores)))

    }
    #process sequentially
    else{
        for(j in 1:number.words)
        {
            L1 <- attributes(m1 <- betabin(cbind(gen.counts[,j], rsums - gen.counts[,j])
                                           ~fdoc.collection - 1, ~1, data=gen.counts,
                                           control = list(maxit = 10, reltol = .001),
                                           fixpar=list(n.clades+1, 0)))$logL
            L0 <- attributes(m0 <- betabin(cbind(gen.counts[,j], rsums - gen.counts[,j])
                                           ~1, ~1, data=gen.counts,
                                           control = list(maxit = 10, reltol = .001),
                                           fixpar=list(2, 0)))$logL

            compare.dist[j] <- 2 * (L1 - L0)
            #let them know how it's doing
            percent.complete <- j / number.words * 100
            #speed up printing, etc.
            if(j %% 2 == 0){
                print(paste(round(percent.complete, digits=2), "% complete | ",
                            (format(Sys.time() - start.t, usetz = TRUE)), sep=""))
                                        #gc()
            }
        }
    }

    #tell the user how long this shit took.
    print(format(Sys.time() - start.t, usetz = TRUE))
    return(compare.dist)
}





##Revamped version of Sydney's function to analyze nclade > 2
getSydneyScore <- function(gen.data, fdoc.collection, n.clades){

  #find number of words based on number of cols in gen.data
  n.words <- dim(gen.data)[2]
  #create vector for means of each clade
  med.clades <- NULL
  #create matrix w/ relative frequencies & a col showing which
  #clade which chunk is in
  syd.data <- data.frame(gen.data, fdoc.collection)
  #set the names of syd.data
  names(syd.data)[length(colnames(syd.data))] <- "clade"

  #create an array of length number.words to store distances for each word
  syd.dist <- array(dim=n.words)



  #loop for each word in gen.data
  for(i in 1:n.words){
    med.clades <- NULL
    #loop for each clade in each word
    for(j in levels(fdoc.collection)){
      #get the median when the i = clade
      med.clades <- c(med.clades, median(syd.data[syd.data$clade == j, i]))
    }

    #get the medians of all of the clade medians
    grand.median <- median(med.clades)
    #calculate clade effects
    clade.effects <- med.clades - grand.median
    gen.data.resids <- gen.data[,i] - med.clades[syd.data$clade]
    #find sum of treatment groups
    sum.trt <-sum(abs(clade.effects[syd.data$clade]))
    #get sum of error between treatment groups
    sum.error <- sum(abs(gen.data.resids))
    #check to see if error is close to 0 (since it's a denom)
    if(sum.error < .0000000001){
      syd.dist[i] <- -10.1
      warning(paste("Not enough variation in data set within clades for word", names(gen.data)[i],
                 "to calculate ABSOLUTE DISTANCE"))
    }
    else{
      #store the distance score for the word
      syd.dist[i] <- (sum.trt/sum.error)*((length(fdoc.collection)-n.clades)/(n.clades-1))
    }
  }

  syd.dist[syd.dist==-10.1] <- max(syd.dist)+1

  #return matrix of all distance scores to calling function
  return(syd.dist)
}


lineColor <- function(x, colorOfNodes)
{
	attr(x, "nodePar") <- list("pch"  = NA, "lab.col"  = colorOfNodes[[currentNode]])
	attr(x, "edgePar") <- list("col"  = colorOfNodes[[currentNode]])
	assign("currentNode",  currentNode + 1, envir = .GlobalEnv)

	x #this line is necessary for some reason for the dendrapply function that calls this to work
}

#get the color for a given leaf based on it's label
getColor <- function(label, specialLabels, metaTable = NULL, colors)
{

    return(colors[label])
}

#generates an list containing the color for every node in the tree
#rules are a node is red if it only contains Archeia, green if it only contains Bacteria, Gold if it was specially selected for highlighting, and blue if it contains mulitple of the previous categories'
#The list is ordered in the order that nodes are visited by dendrapply.
generateLineColorList <- function(x, mergeTableRow, specialLabels,
                                  metaTable = NULL, colors)
{
	colorlist <- list()

	#color the left half of the clade
	if(x$merge[mergeTableRow,1] < 0) #if the left node is a chunk determine the chunk's color
	{
		leftColor <-
	getColor(x$labels[-x$merge[mergeTableRow,1]],
	specialLabels=specialLabels, metaTable = metaTable, colors) #the color of the chunk
		leftList <- list(leftColor) #list of the colors of all the nodes to the left
	}

	else #if the left node is a clade recursively run the function on that clade
	{
		result <- generateLineColorList(x,
	x$merge[mergeTableRow,1], specialLabels=specialLabels,
	metaTable = metaTable, colors)
		leftColor <- result$color #the overall color of the subclade
		leftList <- result$colorList #list of the colors of all the nodes to the left
	}

	#color the right half of the clade
	if(x$merge[mergeTableRow,2] < 0) #if the right node is a chunk determine the chunk's color
	{
		rightColor <-
	getColor(x$labels[-x$merge[mergeTableRow,2]],
	specialLabels=specialLabels, metaTable = metaTable, colors) #the color of the chunk
		rightList <- list(rightColor) #list of the colors of all the nodes to the right
	}

	else #if the right node is a clade recursively run the function on that clade
	{
		result <- generateLineColorList(x,
	x$merge[mergeTableRow,2], specialLabels=specialLabels,
	metaTable = metaTable, colors)
		rightColor <- result$color #the overall color of the subclade
		rightList <- result$colorList #list of the colors of all the nodes to the right
	}

	if(leftColor == rightColor) #check if the colors of the two subclades of the current clade are the same
	{
		color <- leftColor #if so use the color they share
	}

	else #if the colors are different the subclades have different contents
	{
		color <- "black" #set the clade to blue to mark it's mixed contents
	}

	#the colors found need to be put together in the proper order. The current clade has one node for each of it's childern which contains a clade instead of just a chunk.
	#Those nodes need to be given the color of the current clade, but only if they exist. These nodes will appear in the list of colors before all the colors for the nodes in the respective
	#subclades

	if(x$merge[mergeTableRow,1] > 0  && x$merge[mergeTableRow,2] > 0) #if both childern are subclades
	{
		colorList <- c(color, leftList, color, rightList) #both nodes in the current clade exist so add them into the color list
	}

	else if(x$merge[mergeTableRow,1] > 0) #if the right child is a chunk
	{
		colorList <- c(color, leftList, rightList) #there is only a node for the left clade so add that to the color list
	}

	else if(x$merge[mergeTableRow,2] > 0) #if the left child is a chunk
	{
		colorList <- c(leftList, color, rightList) #there is only a node for the right clade so add that to the color list
	}

	else #both children are individual chunks
	{
	colorList <- append(leftColor, rightColor)
	}

	result <- list(colorList=colorList, color=color)
	return(result)
}

#plots a pvclust object
plot.trueTree <- function(x, outputFilename = NULL, print.pv=TRUE, print.num=TRUE, float=0.01,
                         col.pv=c(2,3,8), cex.pv=0.8, font.pv=NULL,
                         col=NULL, cex=NULL, font=NULL, lty=NULL, lwd=NULL,
                         main=NULL, sub=NULL, xlab=NULL, height=800,
                          width=800, specialLabels=NULL, showBP=FALSE,
                          dataset.name, fdoc.collection, ...)
{
    gets.color <- length(levels(fdoc.collection)) <= 6
    if(gets.color){
        colors.vector <- c("red4", "dark green", "navy blue", "mediumpurple", "orange red", "darkgoldenrod4")
        colors <- ifelse(fdoc.collection=="clade1", colors.vector[1],fdoc.collection)
        colors <- ifelse(fdoc.collection=="clade2",colors.vector[2], colors)
        colors <- ifelse(fdoc.collection=="clade3",colors.vector[3],  colors)
        colors <- ifelse(fdoc.collection=="clade4",colors.vector[4],colors)
        colors <- ifelse(fdoc.collection=="clade5",colors.vector[5],  colors)
        colors <- ifelse(fdoc.collection=="clade6",colors.vector[6],  colors)
    }
    else{
        colors <- rep("black", length=length(fdoc.collection))
    }
    names(colors) <- names(fdoc.collection)


    pdf(file=paste(outputFilename, ".pdf", sep=""), paper="a4r", width=11, height=8.5)

    metaTable <- x$metaTable[[1]] #get metadata out of pvclust object

    main <- paste(dataset.name,  paste("Cluster method: ", x$method, sep=""), paste("Distance: ", x$dist.method), sep = "\n")


  if(is.null(sub))
    #sub=paste("Cluster method: ", x$hclust$method, sep="")

  if(is.null(xlab))
    #xlab=paste("Distance: ", x$hclust$dist.method)

  dend <- as.dendrogram(x) #convert the hclust object into a dendrogram object
  colorList <- generateLineColorList(x, dim(x$merge)[1],
  specialLabels=specialLabels, metaTable = metaTable, colors) #figure out what color each node should be
  colorList <- c(0, colorList$colorList) #the first node checked be dendrapply doesn't seem to be part of the dendrogram so add a dummy value at the start of the list

  assign("currentNode",  1, envir = .GlobalEnv) #currentNode is a global variable to keep track of where in the tree we are
  dend <- dendrapply(dend, lineColor, colorList) #add color to all the nodes in the tree

  #find length of longest chunk name
  maxL <- max( nchar( x$labels ))

  # set margins so there is just enough room for the labels
  # The numbers measure margin size in line units
  # The paramets are the size of the bottom,left,top,right margins
  # On average a margin one line wide seems to have room for about 2.5 characters)
  # so the margin on the bottom is set to the number of lines necessary to display
  # the longest label if there was only 2 characters per line which leave's a decent buffer
  par( mar=c((maxL / 2.0), 2.1, 4.1, 2.1), xpd=TRUE)



  plot(dend, main=main, sub=sub, xlab="", col=col, cex=cex,
       font=font, lty=lty, lwd=lwd, ...)

    if(gets.color){
        legend("topright", inset=c(0, -.1),
               legend = unique(levels(fdoc.collection)),
               fill = colors.vector
               )
    }
  if(!is.null(outputFilename)) #if writing to a file close the connection
  {
	dev.off()
  }
}
