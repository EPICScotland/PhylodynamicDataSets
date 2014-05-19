# example script for initial processing of avian influenza sequences
# S J Lycett
# 3 July 2013

###########################################################
# LOADING LIBRARIES AND CUSTOM FUNCTIONS

# load ape
library(ape)

# use helper functions
# assumes that working directory is C://Users//Samantha Lycett//Documents//
# and that there is a subdirectory "Rcode"
# alternatively just load and compile these files into R first

source("Rcode//getEl.R")
source("Rcode//birdSpecies_and_continents.R")
source("Rcode//calcDecimalDate.R")

###########################################################
# LOAD SEQUENCES

# load the fasta format sequences
# note the sequences names must be defined on the NCBI fasta def line like this:
# >{serotype}_{host}_{accession}_{strain}_{country}_{year}/{month}/{day}_{segname}

rootname 	<- "C://Users//Samantha Lycett//Documents//Sema//h7n9_avian_3_july_2013"
fext		<- ".fa"
fname		<- paste(rootname,fext,sep="")
seqs 		<- read.dna( fname, format="fasta", as.matrix=FALSE)

###########################################################
# EXTRACT AND PROCESS SEQUENCE NAMES

# get the sequence names
taxa	<- as.matrix( attributes(seqs)$names )

# break the sequence names to component parts
subtype <- apply(taxa, 1, getEl, ind=1, sep="_")
host    <- apply(taxa, 1, getEl, ind=2, sep="_")
accn	  <- apply(taxa, 1, getEl, ind=3, sep="_")

# isolate name - note this can go wrong if there are underscores in the isolate name its-self
# note MEGA will replace white spaces with underscores - i.e. it will screw up the names if you are not careful
# in this example this has not happened, but if it does you will need to reconstruct the isolate names
isName  <- apply(taxa, 1, getEl, ind=4, sep="_")

# reconstruct something like this if necessary
isName2 <- apply(taxa, 1, getEl, ind=4, fromEnd=TRUE, sep="_")
inds	  <- which(isName != isName2 )
if (length(inds) > 0) {
	isName[inds] <- paste(isName[inds], isName2[inds], sep="_")
	print("attempted to repair isolate names")
} else {
	print("isolate names ok")
}

country <- apply(taxa, 1, getEl, ind=3, sep="_", fromEnd=TRUE)
dateTxt <- apply(taxa, 1, getEl, ind=2, sep="_", fromEnd=TRUE)

seg	  	<- apply(taxa, 1, getEl, ind=1, sep="_", fromEnd=TRUE)
# seg is of form "4 (HA)"
segNumber 	<- as.integer( apply(as.matrix(seg), 1, getEl, ind=1, sep=" ") )


# calculate decimal date
decDate <- apply(as.matrix(dateTxt), 1, calcDecimalDate_from_yymmdd)

# extract details from isolate names
inds	  	<- which(host != "Human")
bird	 	<- array("human", length(isName))		# human sequences are just A/New York/1/2010 etc

# avian and other animal sequences are A/chicken/Hong Kong/1/2010 etc
bird[inds] 	<- apply( as.matrix(isName[inds]), 1, getEl, ind=2, sep="/")

# if human
place		<- apply(as.matrix(isName), 1, getEl, ind=2, sep="/")
# for non-humans
place[inds] <- apply(as.matrix(isName[inds]), 1, getEl, ind=3, sep="/")

# latin bird order, e.g. gal = galliformes
# obviously wont work for "swine" or "human" etc
birdOrder	<- birdClass(bird)

inds		<- which(birdOrder == "-")
if (length(inds) > 0) {
	birdOrder[inds] <- host[inds]
	print("some bird names not recognised - manually inspection required")
} else {
	print("all bird names are recognised")
}

# get wild and domestic classifications for bird names
inds		<- which(birdOrder != "-")
wildDomestic<- array("-", length(birdOrder))
if (length(inds) > 0) {
	wildDomestic[inds] <- getWildDomestic( bird[inds] )
} else {
	print("no birds to process for wild / domestic")
}

# get continent
continent	<- continentClass(country)
inds		<- which(continent == "-")
if (length(inds) > 0) {
	continent[inds] <- country[inds]
	print("some country names not recognised - manually inspection required")
} else {
	print("all country names are recognised")
}

# now make full table with all details
# do this before any more processing
tbl			<- cbind(taxa, segNumber, accn, isName, subtype, host, birdOrder, wildDomestic, bird, continent, country, place, dateTxt, decDate)
# the column name for taxa doesnt automatically appear
colnames(tbl)[1] 	<- "taxa"

###########################################################
# WRITE INITIAL PROCESSED FIELDS TO FILE

tname <- paste( rootname, "_sequenceNamesTable.csv",sep="")
write.table(tbl, file=tname, col.names=TRUE, row.names=FALSE, sep=",")

############################################################
# SPLIT INTO SEGMENTS AND RENAME
# and rename - dont include the accession number this time

newNames	<- paste(paste(subtype,host,birdOrder,wildDomestic,continent,country,isName,sep="_"), decDate, sep="|")
# suitable for BEAST

newSeqs			  <- seqs
attributes(newSeqs)$names <- newNames

for (s in 1:8) {
	inds	<- which(segNumber == s)

	tname			<- paste(rootname,"_seg",s,"_oldNames_newNames.csv",sep="")
	newTbl		<- cbind(taxa[inds],newNames[inds])
	colnames(newTbl) 	<- c("OriginalNames","NewNames")
	write.table(newTbl, file=tname, sep=",", col.names=TRUE, row.names=FALSE)


	sname			<- paste(rootname,"_seg",s,".fas",sep="")
	newSegSeqs		<- newSeqs[inds]
	write.dna( newSegSeqs, file=sname, format="fasta", nbcol=-1, colsep="")
}








