# script to process Genbank download of African Horse Sickness
# for EPIC
# S. J. Lycett
# 27 Oct 2013
# 2 Dec 2013

library(ape)
source("Rcode//getEl.R")

##########################################################
# helper functions

	getSequencesWithTitle <- function(seqs, title="complete cds", fixed=TRUE, useCase=TRUE) {
		taxa  <- attributes(seqs)$names
		if (useCase) {
			inds  <- grep(title, taxa, fixed=fixed)
		} else {
			inds	<- grep(tolower(title), tolower(taxa), fixed=fixed)
		}
		return( seqs[inds] )
	}

	getRemainingSeqs	   <- function( allSeqs, chosenSeqs ) {
		chosenTaxa	<- attributes(chosenSeqs)$names
		allTaxa	<- attributes(allSeqs)$names
		remainTaxa	<- setdiff(allTaxa, chosenTaxa)
		inds		<- match(remainTaxa, allTaxa)
		return( allSeqs[inds] )
	}

##################################################################################
# Step 1
# Split the downloaded Genbank sequences into vaccine, patent and other sequences

doSplit1 <- FALSE
if (doSplit1) {

	path	<- "D://slycett//epic_data//BACKGROUND_SEQUENCES//"
	path2 <- "AfricanHorseSickness//"
	name  <- "africanHorseSickness_sequences"

	seqs  <- read.dna( paste(path,path2,name,".fasta",sep=""), format="fasta", as.matrix=FALSE)
	taxa  <- attributes(seqs)$names

	# get vaccine sequences
	vacSeqs <- getSequencesWithTitle( seqs, title="vaccine", useCase=FALSE)
	
	# get patent sequences
	patSeqs <- getSequencesWithTitle( seqs, title="patent", useCase=FALSE)

	# get other sequences
	remainSeqs <- getRemainingSeqs( seqs, vacSeqs )
	remainSeqs <- getRemainingSeqs( remainSeqs, patSeqs )
	
	write.dna( vacSeqs, file=paste(path,path2,"africanHorseSickness_vaccine_seqs.fas",sep=""), 
			format="fasta", nbcol=-1, colsep="")

	write.dna( patSeqs, file=paste(path,path2,"africanHorseSickness_patent_seqs.fas",sep=""), 
			format="fasta", nbcol=-1, colsep="")

	# actually did this manually from BioEdit
	write.dna( remainSeqs, file=paste(path,path2,"africanHorseSickness_seqs.fas",sep=""), 
			format="fasta", nbcol=-1, colsep="")

}

##################################################################################
# Step 2
# split the non-vaccine and non-patent sequences from step 1 into segments

doSplit2 <- FALSE
if (doSplit2) {
	path	<- "D://slycett//epic_data//BACKGROUND_SEQUENCES//"
	path2 <- "AfricanHorseSickness//"
	name  <- "africanHorseSickness_seqs"
	seqs	<- read.dna( paste(path,path2,name,".fas",sep=""), format="fasta", as.matrix=FALSE)
	taxa  <- attributes(seqs)$names

	#partialSeqs <- getSequencesWithTitle(seqs, title="partial", useCase=FALSE)
	#completeSeqs<- getSequencesWithTitle(seqs, title="complete", useCase=FALSE)

	segNos <- array(0, length(taxa))

	for (s in 1:10) {
		sinds		  <- grep(paste("segment",s), taxa)
		segNos[sinds] <- s
	}

	# segment 1
	s1inds		<- grep("VP1",taxa)
	s1inds		<- sort(unique(c(s1inds,grep("L1",taxa))))
	segNos[s1inds] 	<- 1

	# segment 2
	s2inds		<- grep("VP2",taxa)
	s2inds		<- sort(unique(c(s2inds,grep("L2",taxa))))
	segNos[s2inds]	<- 2

	# segment 3
	s3inds		<- grep("VP3",taxa)
	s3inds		<- sort(unique(c(s3inds,grep("L3",taxa))))
	segNos[s3inds]	<- 3

	# segment 4
	s4inds		<- grep("VP4",taxa)
	s4inds		<- sort(unique(c(s4inds,grep("M4",taxa))))
	segNos[s4inds]	<- 4

	# segment 5
	s5inds		<- grep("NS1",taxa)
	s5inds		<- sort(unique(c(s5inds,grep("M5",taxa))))
	segNos[s5inds]	<- 5

	# segment 6
	s6inds		<- grep("VP5",taxa)
	s6inds		<- sort(unique(c(s6inds,grep("M6",taxa))))
	segNos[s6inds]	<- 6

	# segment 7
	s7inds		<- grep("VP7",taxa)
	s7inds		<- sort(unique(c(s7inds,grep("S7",taxa))))
	segNos[s7inds]	<- 7

	# segment 8
	s8inds		<- grep("NS2",taxa)
	s8inds		<- sort(unique(c(s8inds,grep("S8",taxa))))
	segNos[s8inds]	<- 8

	# segment 9
	s9inds		<- grep("VP6",taxa)
	s9inds		<- sort(unique(c(s9inds,grep("S9",taxa))))
	segNos[s9inds]	<- 9

	# segment 10
	s10inds		<- grep("NS3", taxa)
	s10inds		<- sort(unique(c(s10inds, grep("S10", taxa))))
	segNos[s10inds] 	<- 10

	for (s in c(1:10,0) ) {
		inds		<- which( segNos==s )
		write.dna( seqs[inds], file=paste(path,path2,"africanHorseSickness_segment",s,".fas",sep=""),
				format="fasta", nbcol=-1, colsep="")
	}
}


##################################################################################
# Step 3
# after manual alignment etc
# count sequences per segment

doCount <- FALSE
if (doCount) {
	path	<- "D://slycett//epic_data//BACKGROUND_SEQUENCES//"
	path2 <- "AfricanHorseSickness//aligned//"
	name  <- "africanHorseSickness_segment"

	for (s in 1:10) {
		sname <- paste(path,path2,name,s,"_al.fas",sep="")
		seqs  <- read.dna( sname, format="fasta", as.matrix=FALSE)
		print( c(s,length(seqs)) )
	}


#[1]  1 7
#[1]  2 55
#[1]  3 22
#[1]  4 21
#[1]  5 63
#[1]  6 32
#[1]  7 56
#[1]  8 44
#[1]  9 7
#[1] 10 205

}

##################################################################################
# Step 4
# extract sequences with named serotypes to separate files

doSerotypes <- FALSE
if (doSerotypes) {
	
	path	<- "D://slycett//epic_data//BACKGROUND_SEQUENCES//"
	path2 <- "AfricanHorseSickness//aligned//"
	path3 <- "AfricanHorseSickness//serotype//"
	path4 <- "AfricanHorseSickness//unal//"
	name  <- "africanHorseSickness_segment"
	fullname <- "africanHorseSickness_sequences_full_info"

	# serotype defined on segment 2
	s <- 2
		sname <- paste(path,path2,name,s,"_al.fas",sep="")
		seqs  <- read.dna( sname, format="fasta", as.matrix=FALSE)
		taxa  <- attributes(seqs)$names
		seqs2 <- getSequencesWithTitle(seqs, title="serotype")
		remain<- getRemainingSeqs(seqs, seqs2)
		
		sname2<- paste(path,path3,name,s,"_serotype.fas",sep="")
		sname3<- paste(path,path3,name,s,"_other.fas",sep="")
		write.dna( seqs2, file=sname2, format="fasta", nbcol=-1, colsep="")
		write.dna( remain, file=sname3, format="fasta", nbcol=-1, colsep="")

	# load genbank file
	gbLines	<- readLines (paste(path,path4,fullname,".gb",sep="") )
	locusLines	<- grep("LOCUS",gbLines)
	accnLines	<- grep("ACCESSION", gbLines)
	endLines	<- grep("//", gbLines)

	s <- 2
		sname <- paste(path,path2,name,s,"_al.fas",sep="")
		seqs  <- read.dna( sname, format="fasta", as.matrix=FALSE)
		taxa  <- attributes(seqs)$names
		gb_id <- apply(as.matrix(taxa), 1, getEl, ind=4, sep="\\|")
		gb_id <- gsub("\\.1", "", gb_id)
		ginds <- unlist(apply(as.matrix(gb_id), 1, grep, gbLines[accnLines]))
		acc_ginds <- accnLines[ginds]
		loc_ginds <- locusLines[ginds]
		end_ginds <- endLines[ginds]

		nrecords  <- length(gb_id)
		for (i in 1:nrecords) {
			recordLines <- gbLines[loc_ginds[i]:end_ginds[i]]
			ft_start	<- grep("FEATURES",recordLines)[1]
			source_start<- grep("source",recordLines)[1]
			seq_start   <- grep("CDS",recordLines)[1]
			strain_sero <- c(grep("serotype",recordLines),grep("strain",recordLines))
			jj		<- which((strain_sero >= source_start) & (strain_sero < seq_start))
			strain_sero <- recordLines[ strain_sero[jj] ][1]
			print( paste(gb_id[i],strain_sero) )
		}

}


