#!/usr/bin/env Rscript

# filter ONT fatsq reads based on elapsed time since beginning (since the min timestamp)
# timestamps in nanopore fastq are:
# start_time=2020-09-23T14:20:24Z

# input: pass 1 fastq file (path) and minutes (integer) as args
# output: write reads that are between min (start time) and minutes to a new file

# algorithm:
# - use sed to get min and max from a fastq, same as in get-times.sh
# - get a vector with timestamps from all reads
# filter vector based on minutes
# use seqkit grep to filter reads
 suppressMessages(
	require(lubridate)
 )
 
 arg <- commandArgs(trailingOnly = TRUE)
 if( length(arg) != 2 ) { stop(" Incorrect number of arguments") }

 fqfile_in <- file.path(arg[1])
 mins <- as.numeric(arg[2])
 fqfile_out <- paste(
 	paste(
 		basename(tools::file_path_sans_ext(fqfile_in)), 
 		paste(mins, "min", sep = ""), 
 		sep = "-"), 
 	"fastq",
 	sep = ".")
 
 # get all times strings from fastq file
 alltimes_cmd <- paste("sed -n 's/.*start_time=//p'", fqfile_in, "| cut -f 1 -d ' ' | sort")
 # as string
 alltimes_str <- system(alltimes_cmd, intern = TRUE)
 # as datetime
 alltimes_dt <- lubridate::as_datetime(alltimes_str)
 mintime <- min(alltimes_dt, na.rm = TRUE)
 cat("min time stamp is", as.character(mintime), "\n")
 
 maxtime <- max(alltimes_dt, na.rm = TRUE)
 cat("max time stamp is", as.character(maxtime), "\n")
 
cat(capture.output( round(difftime(maxtime, mintime, units = "auto"),digits = 2)), "\n")

 filtertime <- mintime + (mins * 60)

 # get strings to use in search pattern, write to temp file
 tempfile <- tempfile()
 writeLines( alltimes_str[alltimes_dt <= filtertime], con = tempfile )

 # and use in seqkit
 system2("seqkit", 
 				 args = c("grep", "-n", "-r", "-f", tempfile, fqfile_in) , 
 				 stdout = fqfile_out
 				 )

cat("Reads filtered to", mins, "minutes written to", fqfile_out, "\n")


