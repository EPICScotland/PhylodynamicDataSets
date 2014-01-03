# script to process Genbank download of African Horse Sickness
# for EPIC
# S. J. Lycett
# 27 Oct 2013

library(ape)

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

doSplit2 <- TRUE
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