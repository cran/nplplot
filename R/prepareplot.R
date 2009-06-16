#   Mega2: Manipulation Environment for Genetic Analysis
#   Copyright (C) 1999-2009 Nandita Mukhopadhyay, Lee Almasy,
#            Mark Schroeder, William P. Mulvihill, Daniel E. Weeks
#  
#   This file is part of the Mega2 program, which is free software; you
#   can redistribute it and/or modify it under the terms of the GNU
#   General Public License as published by the Free Software Foundation;
#   either version 2 of the License, or (at your option) any later
#   version.
#  
#   Mega2 is distributed in the hope that it will be useful, but WITHOUT
#   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#   FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#   for more details.
#  
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#  
#   For further information contact:
#       Nandita Mukhopadhyay
#       e-mail: nandita@pitt.edu, Tel:412-624-7351
#       or
#       Daniel E. Weeks
#       e-mail: weeks@pitt.edu
# 
# ===========================================================================

# CREATE FILES (bed.data.# and/or GG.data.all) CONTAINING MARKER, POSITION AND SCORES
# READ IN MEGA2 OUTPUT FILE (rtable) AND MEGA2 ANNOTATED INPUT MAP FILE (mapfile)


prepareplot <- function(prefix, chrlist=1:24, mapfile, output="both")

{

# Delete file with the same name as genome data file "GG.data.all" and "bed.data.23"

unlink("GG.data.all")
unlink("bed.data.23")


# Read in map file

if (file.exists(mapfile) == TRUE) { map <- read.table(mapfile, header=T) } else {

	warning(paste("File", mapfile, "does not", "exist!", sep=" "))
	return(FALSE)

}

# Find the column of physical position in map file

if(length(grep(".p", names(map), fixed=T)) == 1) { 

	p <- names(map)[grep(".p", names(map), fixed=T)]

} else { 
	warning("Please make a map file with exactly one physical position!")
	return(FALSE)

}


# Read in R table

# Number of R tables
chrlist <- as.character(chrlist)
chr_reg <- chrlist[! chrlist %in% c("23", "24")]

# Process with regular chromosomes

for (i in chr_reg) {

	# Read in R table file

	# What is the R table name
	if (length(strsplit(i, split=NULL)[[1]]) == 1) r_n <- paste(prefix, ".0", i, sep="")
	if (length(strsplit(i, split=NULL)[[1]]) == 2) r_n <- paste(prefix, i, sep=".")

	# Check if R table file exists
	if (file.exists(r_n) == TRUE) { r <- read.table(r_n, header=T) } else {

		warning(paste("R table file", r_n, "does not", "exist!", sep=" "))
		return(FALSE)

	}

	r_bg <- r[r$Marker!="ltype" & r$Marker!="ptype", ]
	r_b <- r_bg[r_bg$Marker!="-", ]

	# Merge two tables

	w_bg <- merge(r_bg, map, by.x="Marker", by.y="Name", all.x=T)
	w_b <- merge(r_b, map, by.x="Marker", by.y="Name", all.x=T)

	# Linear interpolation for the bedgragh table

	f <- approxfun(w_b$Position, w_b[,p])
	w_bg[is.na(w_bg[,p]),p] <- round(f(w_bg$Position[is.na(w_bg[,p])]), 0)
	w_bg_ord <- w_bg[order(w_bg[,p]),]

	# Create the data for plotting bed and bedgraph

	bed.dat <- data.frame(as.character(w_bg_ord[,1]), w_bg_ord[,p], w_bg_ord[,3:length(r)])
	names(bed.dat) <- c("Marker", "Position", names(r)[3:length(r)])

	# Create the data for plotting genomegraph

	gg.dat <- data.frame(rep(i, length(w_bg_ord[,1])), as.character(w_bg_ord[,1]), w_bg_ord[,p], w_bg_ord[,3:length(r)])
	names(gg.dat) <- c("Chromosome", "Marker", "Position", names(r)[3:length(r)])

	if (output == "both" | output == "bed") {
	
		# Define the title of bed.data file

		if (length(strsplit(i, split=NULL)[[1]]) == 1) bed_t <- paste("bed.data", ".0", i, sep="")
		if (length(strsplit(i, split=NULL)[[1]]) == 2) bed_t <- paste("bed.data", i, sep=".")

		write.table(bed.dat, file=bed_t, append=F, quote=F, row.names=F, col.names=T, sep="\t")

	}

	if (output == "both" | output == "GG") {

		if (i == chr_reg[1]) {

			write.table(gg.dat, file="GG.data.all", append=F, quote=F, row.names=F, col.names=T, sep="\t")

		} else {

			write.table(gg.dat, file="GG.data.all", append=T, quote=F, row.names=F, col.names=F, sep="\t")

		}
		
	}

}

# Process with chromosomes 23 and 24

w_bg <- NULL

if ("23" %in% chrlist) {

	# Read in R table file

	r_n_23 <- paste(prefix, "23", sep=".")
	r_n_X <- paste(prefix, "X", sep=".")

	if (file.exists(r_n_23) == FALSE & file.exists(r_n_X) == FALSE) {

		warning(paste("Neither", r_n_23, "or", r_n_X, "exists!", sep=" "))
		return(FALSE)

	}

	if (file.exists(r_n_X) == TRUE & file.exists(r_n_23) == TRUE) {

		warning(paste("Both", r_n_23, "and", r_n_X, "exist!", sep=" "))
		return(FALSE)

	}

	if (file.exists(r_n_23) == TRUE & file.exists(r_n_X) == FALSE) r <- read.table(r_n_23, header=T)

	if (file.exists(r_n_X) == TRUE & file.exists(r_n_23) == FALSE) r <- read.table(r_n_X, header=T)

	r_bg <- r[r$Marker!="ltype" & r$Marker!="ptype", ]
	r_b <- r_bg[r_bg$Marker!="-", ]

	# Merge two tables

	w_bg <- merge(r_bg, map, by.x="Marker", by.y="Name", all.x=T)
	w_b <- merge(r_b, map, by.x="Marker", by.y="Name", all.x=T)

	# Linear interpolation for the bedgragh table

	f <- approxfun(w_b$Position, w_b[,p])
	w_bg[is.na(w_bg[,p]),p] <- round(f(w_bg$Position[is.na(w_bg[,p])]), 0)

}

	
if ("24" %in% chrlist) {
		
	# Read in R table file

	r_n_24 <- paste(prefix, "24", sep=".")
	r_n_XY <- paste(prefix, "XY", sep=".")

	if (file.exists(r_n_24) == FALSE & file.exists(r_n_XY) == FALSE) {

		warning(paste("Neither", r_n_24, "nor", r_n_XY, "exists!", sep=" "))
		return(FALSE)

	}

	if (file.exists(r_n_XY) == TRUE & file.exists(r_n_24) == TRUE) {

		warning(paste("Both", r_n_24, "and", r_n_XY, "exist!", sep=" "))
		return(FALSE)

	}

	if (file.exists(r_n_24) == TRUE & file.exists(r_n_XY) == FALSE) r24 <- read.table(r_n_24, header=T)

	if (file.exists(r_n_XY) == TRUE & file.exists(r_n_24) == FALSE) r24 <- read.table(r_n_XY, header=T)

	r24_bg <- r24[r24$Marker!="ltype" & r24$Marker!="ptype", ]
	r24_b <- r24_bg[r24_bg$Marker!="-", ]

	# Merge two tables

	w24_bg <- merge(r24_bg, map, by.x="Marker", by.y="Name", all.x=T)
	w24_b <- merge(r24_b, map, by.x="Marker", by.y="Name", all.x=T)

	# Linear interpolation for the bedgragh table

	f24 <- approxfun(w24_b$Position, w24_b[,p])
	w24_bg[is.na(w24_bg[,p]),p] <- round(f24(w24_bg$Position[is.na(w24_bg[,p])]), 0)

	# Combine records of chromosome 23 and 24
	w_bg <- rbind(w_bg, w24_bg)

}

if (length(w_bg[,1]) != 0) {

	w_bg_ord <- w_bg[order(w_bg[,p]),]

	# Create the data for plotting bed and bedgraph

	bed.dat <- data.frame(as.character(w_bg_ord[,1]), w_bg_ord[,p], w_bg_ord[,3:length(r)])
	names(bed.dat) <- c("Marker", "Position", names(r)[3:length(r)])

	# Create the data for plotting genomegraph

	gg.dat <- data.frame(rep("X", length(w_bg_ord[,1])), as.character(w_bg_ord[,1]), w_bg_ord[,p], w_bg_ord[,3:length(r)])
	names(gg.dat) <- c("Chromosome", "Marker", "Position", names(r)[3:length(r)])

	if (output == "both" | output == "bed") {
	
		write.table(bed.dat, file="bed.data.23", append=F, quote=F, row.names=F, col.names=T, sep="\t")

	}

	if (output == "both" | output == "GG") {

		if (length(chr_reg) == 0) {

			write.table(gg.dat, file="GG.data.all", append=F, quote=F, row.names=F, col.names=T, sep="\t")

		} else {
			
			write.table(gg.dat, file="GG.data.all", append=T, quote=F, row.names=F, col.names=F, sep="\t")

		}
		
	}

}


return(TRUE)

}
