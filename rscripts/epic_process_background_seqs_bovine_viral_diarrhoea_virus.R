# script to process GenBank download of Bovine Viral Diarrhoea Virus
# for EPIC
# S. J. Lycett
# 29 Nov 2013
# 17 Dec 2013

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

	calcSeqLen	<- function( seq ) {
		return( length(unlist(seq)) )
	}

	getSequenceLengths <- function( seqs ) {
		return( unlist(lapply(seqs, calcSeqLen)) )
	}


##########

# returns an array of (locus,accn,gi,defs,strain,date) from a GB record

	processGBRecord <- function( lines ) {
		ind   <- grep("LOCUS",lines)[1]
		els   <- strsplit(lines[ind]," ")[[1]]
		eln	<- apply(as.matrix(els), 1, nchar)
		jj	<- which(eln > 0)
		locus <- els[eln[jj[2]]]

		defI	<- grep("DEFINITION", lines)[1]
		accnI	<- grep("ACCESSION",lines)[1]
		versI <- grep("VERSION",lines)[1]

		defs	<- lines[defI:(accnI-1)]
		defs	<- gsub("DEFINITION","",defs)
		if (length(defs) > 1) {
			dd <- defs[1]
			for (k in 2:length(defs)) {
				dd <- paste(dd, defs[k], sep=" ")
			}
			defs <- dd
			defs <- gsub("  ", " ", defs)
		}

		accn	<- gsub("ACCESSION", "", lines[accnI])
		accn  <- gsub(" ", "", accn)

		gi	<- strsplit(lines[versI],"GI\\:")[[1]][2]

		find  <- grep("FEATURES", lines)[1]
		ind	<- grep("[Ss]train", lines)
		if ( length(ind) > 0 ) {
			jj	 <- which(ind >= find)
			if (length(jj) > 0) {
				strain <- lines[ind[jj]]
				strain <- strsplit(strain,"=")[[1]][2]
				strain <- gsub("\"","",strain)
			} else {
				strain <- "-"
			}
		} else {
			strain <- "-"
		}

		ind	<- grep("[Cc]ollection_date", lines)
		if ( length(ind) > 0 ) {
			jj	 <- which(ind >= find)
			date <- lines[ind[jj]]
			date <- strsplit(date,"=")[[1]][2]
			date <- gsub("\"","",date)
		} else {
			date <- "-"
		}

		ind	<- grep("[Cc]ountry", lines)
		if ( length(ind) > 0 ) {
			jj	 <- which(ind >= find)
			if (length(jj) > 0) {
				country <- lines[ind[jj]]
				country <- strsplit(country,"=")[[1]][2]
				country <- gsub("\"","",country)
			} else {
				country <- "-"
			}
		} else {
			country <- "-"
		}


		res <- c(locus,accn,gi,defs,country,strain,date)
		return( res )

	}
	


##################################################################################
# Step 1
# Split the downloaded Genbank sequences into vaccine, patent and other sequences


doSplit1 <- FALSE
if (doSplit1) {

	path	<- "D://slycett//epic_data//BACKGROUND_SEQUENCES//"
	path2 <- "BovineViralDiarrhoeaVirus//"
	virus <- "bvdv"
	name  <- paste(virus,"_sequences",sep="")

	seqs  <- read.dna( paste(path,path2,name,".fasta",sep=""), format="fasta", as.matrix=FALSE)
	taxa  <- attributes(seqs)$names

	# get vaccine sequences
	vacSeqs <- getSequencesWithTitle( seqs, title="vaccine", useCase=FALSE)
	
	# get patent sequences
	patSeqs <- getSequencesWithTitle( seqs, title="patent", useCase=FALSE)

	# get attenuated sequences
	attenSeqs <- getSequencesWithTitle( seqs, title="attenuated", useCase=FALSE)

	# get other sequences
	remainSeqs <- getRemainingSeqs( seqs, vacSeqs )
	remainSeqs <- getRemainingSeqs( remainSeqs, patSeqs )
	remainSeqs <- getRemainingSeqs( remainSeqs, attenSeqs )
	
	write.dna( vacSeqs, file=paste(path,path2,virus,"_vaccine_seqs.fas",sep=""), 
			format="fasta", nbcol=-1, colsep="")

	write.dna( patSeqs, file=paste(path,path2,virus,"_patent_seqs.fas",sep=""), 
			format="fasta", nbcol=-1, colsep="")

	write.dna( attenSeqs, file=paste(path,path2,virus,"_attenuated.fas",sep=""),
			format="fasta", nbcol=-1, colsep="")

	write.dna( remainSeqs, file=paste(path,path2,virus,"_seqs.fas",sep=""), 
			format="fasta", nbcol=-1, colsep="")


	# get sequences listed as complete genome
	completeGenomes <- getSequencesWithTitle( remainSeqs, title="complete", useCase=FALSE)

	write.dna( completeGenomes, file=paste(path,path2,virus,"_completeGenome.fas",sep=""), 
			format="fasta", nbcol=-1, colsep="")

	completeSeqLens <- getSequenceLengths( completeGenomes )
	
	# get partial sequences
	remainSeqs <- getRemainingSeqs( remainSeqs, completeGenomes )

	write.dna( remainSeqs, file=paste(path,path2,virus,"_partialGenome.fas",sep=""), 
			format="fasta", nbcol=-1, colsep="")


	# get sequence lengths of remaining
	seqLens <- getSequenceLengths( remainSeqs )

	# get taxa of remaining
	taxa	  <- attributes(remainSeqs)$names

	# split to long and short sequences
	longSeqs <- which(seqLens >= 10000)	
	shortSeqs <- which(seqLens < 10000)

# long seqs - remove these (have already put complete genomes to a separate file
#gi|2169408|dbj|E01149.1| cDNA to genomic RNA of bovine viral diarrhea virus
#gi|340253344|dbj|HV302650.1| JP 2009291203-A/1: Infectious bovine viral diarrhea virus clone 
#gi|92860856|dbj|DD122170.1| Infectious bovine viral diarrhea virus clone
#gi|38565522|gb|AY442521.1| Bovine viral diarrhea virus 2 clone BV 54 polyprotein gene, partial cds
#gi|551465772|dbj|DI231264.1| KR 1020110091579-A/32: BOVINE VIRAL DIARRHEA VIRUS WITH A MODIFIED ERNS PROTEIN
#gi|551465771|dbj|DI231263.1| KR 1020110091579-A/31: BOVINE VIRAL DIARRHEA VIRUS WITH A MODIFIED ERNS PROTEIN
#gi|416143352|dbj|HV932180.1| JP 2011512877-A/7: Methods for Identifying Cells Suitable for Large-Scale Production of Recombinant Proteins 

	# these are probably not correct BVDV
	write.dna( remainSeqs[longSeqs], file=paste(path,path2,virus,"_partial_long.fas",sep=""), 
			format="fasta", nbcol=-1, colsep="")
                                                      
	# these are the interesting ones                                                                    
	write.dna( remainSeqs[shortSeqs], file=paste(path,path2,virus,"_partial_b.fas",sep=""),
			format="fasta", nbcol=-1, colsep="")

	partials <- remainSeqs[shortSeqs]
	seqLens  <- getSequenceLengths( partials )
	taxa     <- attributes(partials)$names

	# get 5' UTR
	seqs_5utr   <- getSequencesWithTitle( partials, title="UTR", useCase=FALSE)
	seqs_5utr_b <- getSequencesWithTitle( getRemainingSeqs(partials, seqs_5utr), title="untrans", useCase=FALSE)

	# check no dups
	intersect( attributes(seqs_5utr)$names, attributes(seqs_5utr_b)$names )

	write.dna( seqs_5utr, file=paste(path,path2,virus,"_5utr.fas",sep=""), 
			format="fasta", nbcol=-1, colsep="")
	write.dna( seqs_5utr_b, file=paste(path,path2,virus,"_5utr.fas",sep=""), append=TRUE,
			format="fasta", nbcol=-1, colsep="")

}


##################################################################################
# Part 2 - UK sequences

doPart2 <- FALSE
if (doPart2) {

	path	<- "D://slycett//epic_data//BACKGROUND_SEQUENCES//"
	path2 <- "BovineViralDiarrhoeaVirus//"
	virus <- "bvdv_uk"
	name  <- paste(virus,"_sequences",sep="")

	seqs  <- read.dna( paste(path,path2,name,".fasta",sep=""), format="fasta", as.matrix=FALSE)
	taxa  <- attributes(seqs)$names
	nseqs <- length(taxa)

	e2_seqs  <- getSequencesWithTitle(seqs, title="glycoprotein", useCase=FALSE)
	utr_seqs <- getSequencesWithTitle(seqs, title="5'", useCase=FALSE)
	npro_seqs<- getSequencesWithTitle(seqs, title="npro", useCase=FALSE)
	npro_seqs2<- getSequencesWithTitle(seqs, title="N-terminal", useCase=FALSE)

	other_seqs <- getRemainingSeqs(seqs, e2_seqs)
	other_seqs <- getRemainingSeqs(other_seqs, utr_seqs)
	other_seqs <- getRemainingSeqs(other_seqs, npro_seqs)
	other_seqs <- getRemainingSeqs(other_seqs, npro_seqs2)

	write.dna( e2_seqs, file=paste(path,path2,virus,"_e2.fas",sep=""),
			format="fasta", nbcol=-1, colsep="")

	write.dna( utr_seqs, file=paste(path,path2,virus,"_5utr.fas",sep=""),
			format="fasta", nbcol=-1, colsep="")

	write.dna( npro_seqs, file=paste(path,path2,virus,"_npro.fas",sep=""),
			format="fasta", nbcol=-1, colsep="")
	write.dna( npro_seqs2, file=paste(path,path2,virus,"_npro.fas",sep=""),
			format="fasta", nbcol=-1, colsep="", append=TRUE)

	write.dna( other_seqs, file=paste(path,path2,virus,"_other.fas",sep=""),
			format="fasta", nbcol=-1, colsep="")


}


##################################################################################
# 17 Dec 2013
# match strain type and short names to reference sequences

doPart3 <- FALSE
if (doPart3) {
	path		<- "D://slycett//epic_data//BACKGROUND_SEQUENCES//"
	path2 		<- "BovineViralDiarrhoeaVirus//"

	tblName 	<- paste(path,path2,"booth_table3_ref_strains_tbl.txt",sep="")
	refLines	<- readLines( tblName )
	colNames    <- strsplit(refLines[1], "\t")[[1]]
	shortName	<- apply(as.matrix(refLines[2:length(refLines)]), 1, getEl, ind=1, sep="\t")
	strainType	<- apply(as.matrix(refLines[2:length(refLines)]), 1, getEl, ind=2, sep="\t")
	accn_5utr	<- apply(as.matrix(refLines[2:length(refLines)]), 1, getEl, ind=4, sep="\t")
	accn_npro	<- apply(as.matrix(refLines[2:length(refLines)]), 1, getEl, ind=5, sep="\t")

	shortName	<- gsub(" ", "", shortName)

	strainType	<- gsub(" ", "", strainType)
	strainType  <- gsub("\\*", "", strainType)
	strainType  <- gsub("\\^", "", strainType)
	strainType  <- toupper(strainType)

	newNames	<- paste(paste("BVDV",strainType,sep="_"),shortName,sep="/")

	accn_5utr	<- gsub(" ", "", accn_5utr)
	accn_5utr	<- gsub("\\+", "", accn_5utr)
	accn_5utr   <- gsub("SuppliedbyProf.Vilcek", "-", accn_5utr, fixed=TRUE)
	accn_5utr	<- gsub("-", "", accn_5utr, fixed=TRUE)
	inds		<- which(accn_5utr != "")
	accn_5utr	<- accn_5utr[inds]
	names_5utr  <- newNames[inds]


	accn_npro	<- gsub(" ", "", accn_npro)
	accn_npro	<- gsub("\\+", "", accn_npro)
	accn_npro   <- gsub("SuppliedbyProf.Vilcek", "-", accn_npro, fixed=TRUE)
	accn_npro	<- gsub("-", "", accn_npro, fixed=TRUE)
	inds		<- which(accn_npro != "")
	accn_npro	<- accn_npro[inds]
	names_npro  <- newNames[inds]



	# load 5 utr and npro sequences
	name		<- "booth_"
	n2		<- "_refs_trim"
	ext		<- ".fas"
	seqs_5utr   <- read.dna( paste(path,path2,name,"5utr",n2,ext,sep=""), format="fasta", as.matrix=FALSE)
	seqs_npro	<- read.dna( paste(path,path2,name,"npro",n2,ext,sep=""), format="fasta", as.matrix=FALSE)

	# get full taxa names
	taxa_5utr   <- attributes(seqs_5utr)$names
	taxa_npro	<- attributes(seqs_npro)$names

	# match to short taxa names
	inds		<- apply(as.matrix(accn_5utr), 1, grep, taxa_5utr)
	new_taxa_5utr <- names_5utr[inds]

	inds		<- apply(as.matrix(accn_npro), 1, grep, taxa_npro)
	new_taxa_npro <- names_npro[inds]

	# make new sequences with short taxa names
	newSeqs_5utr <- seqs_5utr
	attributes(newSeqs_5utr)$names <- new_taxa_5utr

	newSeqs_npro <- seqs_npro
	attributes(newSeqs_npro)$names <- new_taxa_npro

	# write to file
	write.dna( newSeqs_5utr, file=paste(path,path2,name,"5utr_refs_shortNames.fas",sep=""), format="fasta", 
			colsep="", nbcol=-1)
	write.dna( newSeqs_npro, file=paste(path,path2,name,"npro_refs_shortNames.fas",sep=""), format="fasta", 
			colsep="", nbcol=-1)
}


##################################################################################
# 30 Dec 2013
# attempt to get date information

# load complete info file for all downloaded sequences

doGetDates <- FALSE
if (doGetDates) {
	path	<- "D://slycett//epic_data//BACKGROUND_SEQUENCES//"
	path2 <- "BovineViralDiarrhoeaVirus//"
	virus <- "bvdv"
	gbname<- paste(path,path2,virus,"_sequences_full_info.gb",sep="")
	gbLines <- readLines(gbname)

	recordStarts <- grep("LOCUS",gbLines)
	nrecords	 <- length(recordStarts)
	recordEnds   <- c(recordStarts[2:nrecords]-1,length(gbLines))

	tbl		 <- 	matrix(0, nrecords, 7)
	colnames(tbl)<- c("LOCUS","ACCESSION","GI","DEFINITION","Country","Strain","Collection-date")
	for (j in 1:nrecords) {
		res <- processGBRecord( gbLines[ recordStarts[j] : recordEnds[j] ] )
		tbl[j,] <- res
	}

	write.table( tbl, file = paste(path,path2,virus,"_sequences_extract_info_gb.txt",sep=""),
				sep="\t", col.names=TRUE, row.names=FALSE)

	save(tbl, file=paste(path,path2,virus,"_sequences_extract_info_gb_tbl.Rdata",sep="") )

}


##################################################################################
# 31 Dec 2013
# add dates to background sequences

addDatesToUK <- FALSE
if (addDatesToUK) {
	path	<- "D://slycett//epic_data//BACKGROUND_SEQUENCES//"
	path2 <- "BovineViralDiarrhoeaVirus//"
	virus <- "bvdv"

	# load dates table
	load( paste(path,path2,virus,"_sequences_extract_info_gb_tbl.Rdata",sep="") )

	virus <- "bvdv_uk"
	segs  <- c("5utr","npro","e2")

	# load sequences
	for (i in 1:length(segs)) {
		sname <- paste(path,path2,virus,"_",segs[i],"_al.fas",sep="")
		seqs  <- read.dna( sname, format="fasta", as.matrix=FALSE)
		taxa  <- attributes(seqs)$names
		taxa_gi <- apply(as.matrix(taxa), 1, getEl, ind=2, sep="\\|")

		# match gi numbers to table with info and dates
		inds			<- match(taxa_gi, tbl[,3])
		all( taxa_gi == tbl[inds,3] )
		corresponding_tbl <- tbl[inds,]
		

		# which entries have dates
		inds			<- which(corresponding_tbl[,7] != "-")

	if (length(inds) > 0) {

		dated_taxa		<- taxa[inds]
		dated_tbl		<- corresponding_tbl[inds,]
		
		# these taxa seem to have dates but no strain info
		dates			<- dated_tbl[,7]
		strains		<- dated_tbl[,6]
		strains[ which(strains=="-") ] <- "BVDV"

		country		<- dated_tbl[,5]
		country		<- gsub(" ","_",country)

		# extract isolate number from DEFINITION
		isNumber		<- apply(as.matrix(dated_tbl[,4]), 1, getEl, ind=2, sep="isolate")
		isNumber		<- apply(as.matrix(isNumber), 1, getEl, ind=2, sep=" ")
		isolate		<- paste("isolate-",isNumber,sep="")

		# extract year from dates
		# date format seems to be Mon-Year
		days			<- apply(as.matrix(dates), 1, getEl, ind=3, sep="-", fromEnd=TRUE)
		months		<- apply(as.matrix(dates), 1, getEl, ind=2, sep="-", fromEnd=TRUE)
		years			<- apply(as.matrix(dates), 1, getEl, ind=1, sep="-", fromEnd=TRUE)		

		mtxt			<- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
		mm			<- match(months, mtxt)
		yy			<- as.integer(years)
		dd			<- array(15, length(dates))	# days default to 15th of month
		if (length(days) == length(dates) ) {
			jj		<- which( apply(as.matrix(days), 1, nchar) > 0 )
			dd[jj]	<- as.integer( days[jj] )
		} 
		decDates		<- array(0, length(dates))
		for (k in 1:length(dates)) {
			decDates[k] <- calcDecimalDate(dd[k],mm[k],yy[k])
		}

		newIsolateNames	<- paste(strains,country,isolate,years,sep="/")
		newNames		<- paste(newIsolateNames,dates,decDates,sep="|")

		newTbl		<- cbind(dated_taxa,newNames,newIsolateNames,dated_tbl,decDates)
		colnames(newTbl)[1] <- "OriginalSequenceNames"

		sname2 			<- paste(path,path2,virus,"_",segs[i],"_dated_only.fas",sep="")
		inds	 			<- match(dated_taxa, taxa)
		all( dated_taxa == taxa[inds] )
		seqs2  			<- seqs[inds]
		attributes(seqs2)$names <- newNames
		write.dna( seqs2, file=sname2, format="fasta", nbcol=-1, colsep="")
		
		tname2			<- paste(path,path2,virus,"_",segs[i],"_dated_only_info_tbl.txt",sep="")
		write.table( newTbl, file=tname2, sep="\t", col.names=TRUE, row.names=FALSE)

	  } else {
		print(paste("No dates for",segs[i]))
	  }
	}	

}

#############################################################################################
# 2 Jan 2014
# split dated UK sequences to BVDV type

doSplitUKtoType <- FALSE
if (doSplitUKtoType) {

	source("Rcode//selected_sequences_from_figTree.R")

	path	<- "D://slycett//epic_data//BACKGROUND_SEQUENCES//"
	path2 <- "BovineViralDiarrhoeaVirus//"	
	name 	<- "bvdv_uk_5utr_dated_only_and_new_archive"
	trExt <- "_tn93_nj_colByType"
	figTreeName <- paste(path,path2,name,trExt,".figTree",sep="")
	seqsName <- paste(path,path2,name,".fas",sep="")

	selected_sequences_from_figTree( seqsName, figTreeName, mega=TRUE )

#[1] "Sequences & FigTree taxa match"
#[1] "Number of sequences 425 D://slycett//epic_data//BACKGROUND_SEQUENCES//BovineViralDiarrhoeaVirus//bvdv_uk_5utr_dated_only_and_new_archive_color65536.fas"
#[1] "Number of sequences 20 D://slycett//epic_data//BACKGROUND_SEQUENCES//BovineViralDiarrhoeaVirus//bvdv_uk_5utr_dated_only_and_new_archive_color65332.fas"
#[1] "Number of sequences 63 D://slycett//epic_data//BACKGROUND_SEQUENCES//BovineViralDiarrhoeaVirus//bvdv_uk_5utr_dated_only_and_new_archive_color16776961.fas"
 
#425 => BVDV_1A
#20 => BDVD_1I
#63 => BVDV_1B

# reload sequences and tabulate by year etc
	types <- c("BVDV_1A","BVDV_1B","BVDV_1I")
	ext 	<- "_uk_5utr_dated_only_and_new_archive.fas"

	all_st 	<- c()
	all_years 	<- c()

	for (i in 1:length(types)) {
		sname <- paste(path,path2,types[i],ext,sep="")
		seqs  <- read.dna( sname, format="fasta", as.matrix=FALSE)
		taxa  <- as.matrix(attributes(seqs)$names)
		st    <- apply(taxa, 1, getEl, ind=1, sep="/")
		decDate <- apply(taxa, 1, getEl, ind=1, fromEnd=TRUE, sep="\\|")
		dateTxt <- apply(taxa, 1, getEl, ind=2, fromEnd=TRUE, sep="\\|")
		year    <- floor(as.numeric(decDate))
		month	  <- apply(as.matrix(dateTxt), 1, getEl, ind=2, fromEnd=TRUE, sep="-")

		jj	  <- which(st == "BVDV")
		kk	  <- which(st != "BVDV")
		st[jj]  <- paste(types[i],"(genBank)",sep="")
		st[kk]  <- paste(types[i],"(new)",sep="")

		all_st  <- c(all_st, st)
		all_years <- c(all_years, year)
	}
	
	write.table( table(all_st,all_years), col.names=TRUE, row.names=TRUE,
				sep="\t", file=paste(path,"subtypes_table",ext,"_table.txt",sep="") )

}


#############################################################################################
# 7 Jan 2014
# extract Booth et al sequences which have locations

doBoothLocations <- FALSE
if (doBoothLocations) {

	path	<- "D://slycett//epic_data//BACKGROUND_SEQUENCES//"
	path2 <- "BovineViralDiarrhoeaVirus//"

	tbl	<- read.table( paste(path,path2,"booth_table_1.txt",sep=""),
					sep="\t", header=TRUE) 

	virus <- "bvdv_uk"
	segs  <- c("5utr","npro","e2")

	for (s in 1:length(segs)) {
		seqs <- read.dna( paste(path,path2,virus,"_",segs[s],"_al.fas",sep=""),
				format="fasta", as.matrix=FALSE)
		taxa <- attributes(seqs)$names
		isName <- apply(as.matrix(taxa),1,getEl,ind=2,sep="_isolate")
		isName <- apply(as.matrix(isName), 1, getEl, ind=2, sep="_")
		farmNo <- apply(as.matrix(isName), 1, getEl, ind=1, sep="-")
		isNo	 <- apply(as.matrix(isName), 1, getEl, ind=2, sep="-")
		farmNo <- as.integer(farmNo)
		isNo	 <- as.integer(isNo)

		inds	 <- match(farmNo, tbl[,1])
		jj	 <- which(is.finite(inds))
		all( farmNo[jj] == tbl[inds[jj],1])

		if (length(jj) > 0) {

		selInds <- jj
		selTaxa <- taxa[jj]
		selIsName <- isName[jj]
		selSeqs <- seqs[jj]
		selInfo <- tbl[inds[jj],]

		host	  <- paste(selInfo[,3])
		host[ which(host=="Unknown") ] <- "Cattle"

		location<- paste(selInfo[,2])
		location<- gsub(" ","",location)

		selNewName <- paste("BVDV/",host,"/UK/",location,"/",selIsName,sep="")
		selInfo <- cbind(selTaxa,selNewName,selInfo)

		attributes(selSeqs)$names <- selNewName

		rname	  <- paste(path,path2,virus,"_booth_loc_",segs[s],sep="")
		sname	  <- paste(rname,".fas",sep="")
		tname	  <- paste(rname,"_sequenceInfo.txt",sep="")

		write.dna( selSeqs, file=sname, format="fasta", nbcol=-1, colsep="")
		write.table( selInfo, file=tname, sep="\t", col.names=TRUE, row.names=FALSE)

		} else {
			print( paste("No sequences for",segs[s]) )
		}

	}

	# now concatenate where possible
	# Booth has only 5utr and npro

	seqs_5utr <- read.dna( paste(path,path2,virus,"_booth_loc_5utr.fas",sep=""),
				format="fasta", as.matrix=FALSE)
	seqs_npro <- read.dna( paste(path,path2,virus,"_booth_loc_npro.fas",sep=""),
				format="fasta", as.matrix=FALSE)

	t1	    <- attributes(seqs_5utr)$names
	t2	    <- attributes(seqs_npro)$names

	joint	    <- intersect(t1,t2)

	inds	    <- match(joint,t1)
	write.dna( seqs_5utr[inds], file=paste(path,path2,virus,"_booth_loc_common_5utr.fas",sep=""),
				format="fasta", nbcol=-1, colsep="")

	inds	    <- match(joint,t2)
	write.dna( seqs_npro[inds], file=paste(path,path2,virus,"_booth_loc_common_npro.fas",sep=""),
				format="fasta", nbcol=-1, colsep="")

	library(seqinr)

	jname <- paste(path,path2,virus,"_booth_loc_joint_5utr_npro.fas",sep="")
	for (i in 1:length(joint)) {
		i1 <- which(t1 == joint[i])
		i2 <- which(t2 == joint[i])
		s1 <- c2s( unlist(as.character(seqs_5utr[i1])) )
		s2 <- c2s( unlist(as.character(seqs_npro[i2])) )
		ss <- paste(s1,"---------",s2,sep="")
		write( paste(">",joint[i],sep=""), file=jname, append=(i>1))
		write( ss, file=jname, append=TRUE)
	}

	# table(location)
#location
#          Ayr BuryStEdmunds      Carlisle    Carmarthen    Chelmsford    Galashiels 
#            2            12             3             2            13             5 
#    Inverness     Inverurie     Newcastle       Reigate       Taunton 
#            9             7             8             1            36


	# manually looked up locations
	# load locations info

	loc_info <- read.table( paste(path,path2,"bvdv_uk_booth_locations.txt",sep=""), header=TRUE, sep="\t")

	# load unique sequences (5utr & npro)
	seqs <- read.dna( paste(path,path2,virus,"_booth_loc_joint_5utr_npro_unique.fas",sep=""),
			format="fasta", as.matrix=FALSE)
	taxa <- attributes(seqs)$names
	locs <- apply(as.matrix(taxa), 1, getEl, ind=4, sep="/")
	host <- apply(as.matrix(taxa), 1, getEl, ind=2, sep="/")

	minds 	<- match(locs, loc_info[,1])
	all( locs == loc_info[minds,1])
	country 	<- loc_info[minds,2]
	region	<- loc_info[minds,3]
	lat		<- loc_info[minds,4]
	lon		<- loc_info[minds,5]	

#> table(country)
#country
# England Scotland    Wales 
#      44       13        1

	country2 <- paste(country)
	country2[ which(country != "Scotland") ] <- "RUK"

#> table(country2)
#country2
#     RUK Scotland 
#      45       13 


#> table(region)
#region
# NorthEngland NorthScotland  SouthEngland SouthScotland          West 
#            4             6            17             7            24 

	region2 <- paste(region)
	region2[ grep("North", region) ]    <- "North"
	region2[ grep("Scotland", region) ] <- "North"
	region2	<- gsub("SouthEngland","South",region2)

#> table(region2)
#region2
#North South  West 
#   17    17    24 

	traits_tbl <- cbind(taxa,paste(country),country2,paste(region),region2,lat,lon)
	colnames(traits_tbl) <- c("traits","cont","country","reg","region","lat","lon")

	tname <-  paste(path,path2,virus,"_booth_loc_joint_5utr_npro_unique_traits.txt",sep="")
	write.table(traits_tbl,file=tname,sep="\t",col.names=TRUE,row.names=FALSE)
}


#############################################################################################
# 8 Jan 2014
# ran BEAST discrete traits with Booth locations
# to process rate matrix

doBoothRates <- FALSE
if (doBoothRates) {

	path1	<- "D://slycett//epic_data//BACKGROUND_SEQUENCES//"
	path2 <- "BovineViralDiarrhoeaVirus//"
	path3	<- "beast//no_dates_discrete_traits//"

	path <- paste(path1,path2,path3,sep="")

	name <- "booth_joint_asym_reg5_1"
	trait <- "reg"

	rate_hpds <- get_BEAST_180_relative_rates(path=path,name=name,trait=trait,burnin=1000,step=1)
	rel_net   <- get_rel_network(rate_hpds)
	write.table(rel_net, file=paste(path,name,".",trait,".relRates.net",sep=""),
				sep="\t", col.names=FALSE, row.names=FALSE)



	name <- "booth_joint_asym_regions_2"
	trait <- "region"

	rate_hpds <- get_BEAST_180_relative_rates(path=path,name=name,trait=trait,burnin=1000,step=1)
	rel_net   <- get_rel_network(rate_hpds)
	write.table(rel_net, file=paste(path,name,".",trait,".relRates.net",sep=""),
				sep="\t", col.names=FALSE, row.names=FALSE)


	name <- "booth_joint_asym_country_2"
	trait <- "country"

	rate_hpds <- get_BEAST_180_relative_rates(path=path,name=name,trait=trait)
	rel_net   <- get_rel_network(rate_hpds)
	write.table(rel_net, file=paste(path,name,".",trait,".relRates.net",sep=""),
				sep="\t", col.names=FALSE, row.names=FALSE)



	name  <- "booth_joint_asym_place_1"
	trait <- "place"

	rate_hpds <- get_BEAST_180_relative_rates(path=path,name=name,trait=trait,burnin=1000,step=1)
	rel_net   <- get_rel_network(rate_hpds)
	write.table(rel_net, file=paste(path,name,".",trait,".relRates.net",sep=""),
				sep="\t", col.names=FALSE, row.names=FALSE)

	name  <- "booth_joint_sym_place_1"
	trait <- "place"

	rate_hpds <- get_BEAST_180_relative_rates(path=path,name=name,trait=trait,burnin=1000,step=1)
	rel_net   <- get_rel_network(rate_hpds)
	write.table(rel_net, file=paste(path,name,".",trait,".relRates.net",sep=""),
				sep="\t", col.names=FALSE, row.names=FALSE)


	name  <- "booth_joint_asym_bssvs_place_1"
	trait <- "place"

	rate_hpds <- get_BEAST_BSSVS_Rates_174(path=path,name=name,burnin=1000,step=1)
	rel_net   <- get_network(rate_hpds, thres=0.5)
	write.table(rel_net, file=paste(path,name,".",trait,".relRates.net",sep=""),
				sep="\t", col.names=FALSE, row.names=FALSE)

	rel_net   <- get_network(rate_hpds, thres=0.0)
	write.table(rel_net, file=paste(path,name,".",trait,".relRates_thres0.net",sep=""),
				sep="\t", col.names=FALSE, row.names=FALSE)


# reminder of traits
	traits_tbl <- read.table( paste(path,"bvdv_uk_booth_loc_joint_5utr_npro_unique_traits.txt",sep=""),
			sep="\t", header=TRUE)

	table(traits_tbl$country)

#table(traits_tbl$country)
#
#     RUK Scotland 
#      45       13


#table(traits_tbl$reg)
#
#NorthEngland NorthScotland  SouthEngland SouthScotland          West 
#            4             6            17             7            24 

 place <- apply(as.matrix(traits_tbl[,1]), 1, getEl, ind=4, sep="/")

#> table(place)
#place
#          Ayr BuryStEdmunds      Carlisle    Carmarthen    Chelmsford    Galashiels     Inverness     Inverurie 
#            2             8             3             1             8             5             2             4 
#    Newcastle       Reigate       Taunton 
#            1             1            23 

}

