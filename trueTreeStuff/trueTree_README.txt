----------------------------------------------------------------
--------------------------- trueTree ---------------------------

trueTree is computationally intensive.  It may be best to have a
machine dedicated to analysis, or that isn't used by anything else
while it is being run.  To help cut down on runtime, trueTree makes
use of the "snowfall" and "snow" packages in R for parallel
processing. In applications where this is impractical, it can be run
sequentially.

Because of the complicated nature of the analysis, a direct call to
trueTree can have up to 26 input parameters.  Despite this, however,
most Lexomics work only needs 11 of those.  To make it easier to use,
a simplified wrapper file was written.  The wrapper includes the code
for sourcing and calling the trueTree function, all you need to do is
set the basic parameters for it to use.  To make use of the more
complex features of trueTree, you'll have to call it directly.

Each step of running the wrapper is specified by comments in the
file.  You can either copy the commands to another .R file to modify
them or modify them directly in the wrapper function and run them from there.

Most of the commands in the wrapper are essentially strings to change
for where your files are locally, and where you want the output to
go.

Things you NEED to set to run the function (need to be set as strings):
       -The path to a folder containing all of the trueTree source files
       -The name of the input file (including the path to it if it's in
         another directory
       -The name you want to use for the output dendrogram and histogram OR
	 both the path to an output directory you want to use and the name for
	 the output files.
       -The title to use for the dendrogram (can be set to NULL)
       -The metrics and linkage method to use to create the distance
         matrix
Things you may need to change:
       -Distance metric and cluster method (if you know you need to
         change them, you know what these are.
       -Whether or not the input file is transposed (if chunks are
         rows and words are columns, leave this as TRUE, if chunks are
         columns and words are rows, set this to FALSE)
       -Number of times to run the nboot.  For quick & dirty analyses,
         1000 may give workable results.  10000 tends to have the best
         precision/performance ratio, though.
       -Whether or not to run the function in parallel (for large
         dataset, this is recommended).
       -number of CPUs to run on if running in parallel
       -width & height of the output dendrogram (you may need stretch
         the width so the AU and BP values are readable).

Whether or not you're changing the optional values, you need to
evaluate the code in R to set them to their defaults.

After you set all of the parameters, evaluate the code calling the
trueTree function.  If it's being run sequentially, there will be
output showing which r value is being executed.  If it's being run in
parallel, there will be no output until the code is finished being
evaluated (if you're really concerned, you can watch the CPUs in your
system's performance monitor).  The output files will automatically be
created.

The final dendrogram is stored to a trueTree object called "result,"
so if you find the dendrogram is hard to read, you can increase the height
and width values and call plot.trueTree directly to re-create it
larger or smaller (or with a different filename, title, etc).
