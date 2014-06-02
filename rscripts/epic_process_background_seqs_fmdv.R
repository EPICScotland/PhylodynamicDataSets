# script to process Genbank download of Foot and Mouth Disease Viruses
# for EPIC
# S. J. Lycett
# 28 May 2014
# 1 June 2014

# search term = "Foot-and-mouth disease virus"[porgn:__txid12110] 

library(ape)
source("Rcode//getEl.R")
source("Rcode//calcDecimalDate.R")

#########################################################################################################################
# supporting functions

getFeaturesTable <- function( rec ) {
	ind1 <- grep("FEATURES", rec )
	ind2 <- grep("ORIGIN", rec )
	ft   <- rec[ind1:(ind2-1)]
	return( ft )
}

numberCDS <- function( ft ) {
	ind <- grep("CDS", ft)
	return( length( ind ) )
}



getCDS <- function( ft ) {
	ind <- grep("CDS", ft)
	cds <- ft[ind:length(ft)]
	return( cds )
}

getCDSProduct <- function( ft, product="p72" ) {
	inds 		<- grep("CDS", ft)
	prods		<- grep("product=",ft)
	jj   		<- grep(product, ft[prods] )
	kk   		<- prods[jj]
	temp 		<- which(inds < kk )
	if ( length(temp) > 0 ) {
		cds_start 	<- inds[ temp[ length(temp) ] ]
		cds_end	<- inds[ temp[ length(temp) ] + 1 ] -1
		if ((cds_end == cds_start) | (cds_end > length(ft))) {
			cds_end <- length(ft)
		}
	} else {
		cds_start <- inds[1]
		cds_end   <- length(ft)
	}
	return( ft[cds_start:cds_end] )
}

getFirstCDSProduct <- function( ft ) {
	inds 		<- grep("CDS", ft)
	if (length(inds) > 1) {
		cds_start <- inds[1]
		cds_end   <- inds[2]-1
	} else {
		cds_start <- inds[1]
		cds_end   <- length(ft)
	}
	return( ft[cds_start:cds_end] )
}

getInfoFromFeaturesTable	<- function( name="collection_date", ft=ft ) {
	ind <- grep(name, ft)
	if (length(ind) > 0) {
		ll  <- ft[ind]
		els <- strsplit(ll, "=")[[1]]
		ans <- els[2]
		ans <- gsub("\"", "", ans)
		return( ans )
	} else {
		return ( -1 )
	}
}

getLocus		<- function( rec ) {
	els		<- strsplit(rec[1]," ")[[1]]
	nn    	<- apply(as.matrix(els), 1, nchar)
	jj		<- which(nn > 0)
	els		<- els[jj]
	locus 	<- els[2]
	bp		<- els[3]
	recDate 	<- els[ length(els) ]
	res		<- c(locus,bp,recDate)
	return (res)
}

getGenBankInfo	<- function( rec, names1 =c("host","country","collection_date","serotype","isolate","strain","segment"),
						names2=c("gene","product","codon_start","note"), product="x" ) {
	loc	<- getLocus( rec )
	ft 	<- getFeaturesTable( rec )
	ncds  <- numberCDS( ft )
	
  if (ncds == 0) {

	ans			<- apply( as.matrix(names1), 1, getInfoFromFeaturesTable, ft=ft )
	ans2			<- array("-1", length(names2))
	ans			<- c(loc, ans, ans2, "-1", "-1")
	ans			<- as.matrix(ans, length(ans), 1)
	rownames(ans) 	<- c("LOCUS","LENGTH","RECORD_DATE",names1,names2,"COMPLEMENT","position")

  } else {

	if (ncds > 1) {
		if (product != "x") {
			cds <- getCDSProduct(ft, product=product)
		} else {
			cds <- getFirstCDSProduct(ft)
		}
	} else {
		cds   <- getCDS( ft )
	}
	cinds <- grep("complement",cds)
	if ( length(cinds) > 0 ) {
		complement <- TRUE
	} else {
		complement <- FALSE
	}

	cels  <- strsplit(cds[1]," ")[[1]]
	temp  <- apply(as.matrix(cels), 1, nchar)
	jj    <- which(temp > 0)
	cels  <- cels[jj]
	pos   <- cels[ length(cels) ]

	ans	<- apply( as.matrix(names1), 1, getInfoFromFeaturesTable, ft=ft )
	ans2	<- apply( as.matrix(names2), 1, getInfoFromFeaturesTable, ft=cds )

	ans	<- c(loc, ans, ans2, complement, pos)
	ans	<- as.matrix(ans, length(ans), 1)
	rownames(ans) <- c("LOCUS","LENGTH","RECORD_DATE",names1,names2,"COMPLEMENT","position")

   }

	return ( ans )
}

getNuclsFromGenBank <- function( rec ) {
	loc	     <- getLocus( rec )[1]

	ii	     <- grep("ORIGIN",rec)[1]
	nuclsLines <- rec[(ii+1):(length(rec)-1)]
	nuclsLines <- gsub("[0-9]", "", nuclsLines)
	nuclsLines <- gsub(" ", "", nuclsLines)
	nuclsLines <- gsub("//", "", nuclsLines)

	temp	     <- c()
	for (j in 1:length(nuclsLines)) {
		temp <- paste(temp, nuclsLines[j],sep="")
	}
	nuclsLines <- toupper(temp)

	ans	     <- c(loc,nuclsLines)
	ans	     <- as.matrix(ans, length(ans), 1)
	rownames(ans) <- c("LOCUS","NUCLS")
	return( ans )
}

#########################################################################################################################
# process the data (global data)

	path	<- "D://slycett//epic_data//BACKGROUND_SEQUENCES//"
	path2 <- "FootAndMouthDiseaseVirus//2014_may_28//"
	name  <- "fmdv_sequence_genBank_full"
	
	gblines 		<- readLines( paste(path,path2,name,".gb",sep="") )
	recordStarts 	<- grep("LOCUS",gblines)
	nrec			<- length(recordStarts)
	recordEnds		<- c(recordStarts[2:nrec]-1,length(gblines))

	recs			<- vector("list",nrec)
	for (i in 1:nrec) {
		recs[[i]] <- gblines[ recordStarts[i]:recordEnds[i] ]
	}

	nCDS			<- unlist( lapply(recs, numberCDS) )
	jj			<- which(nCDS > 1)

	infos			<- lapply(recs, getGenBankInfo)
	nels			<- length(infos[[1]])
	header		<- rownames(infos[[1]])
	infos 		<- matrix( unlist(infos), nels, length(infos))
	infos			<- t(infos)
	colnames(infos)	<- header

	# write all info to file
	write.table( infos, file=paste(path,path2,name,"_process_info.txt",sep=""), sep="\t", col.names=TRUE, row.names=FALSE)

	all_locii		<- infos[,1]

	# write info to file which has at least 1 cds
	inds 			<- which(nCDS >= 1)
	infos			<- infos[inds,]

	write.table( infos, file=paste(path,path2,name,"_process_info_with_cds.txt",sep=""), sep="\t", col.names=TRUE, row.names=FALSE)
	
	###############
	# re-load table
	infos	 <- read.table( paste(path,path2,name,"_process_info_with_cds.txt",sep=""), sep="\t", header=TRUE )
	header <- colnames(infos)
	infos  <- as.matrix(infos)
	

	# find the individual segments / genes
	gene_col		<- which(header=="gene")
	prod_col		<- which(header=="product")
	seg_col		<- which(header=="segment")
	note_col		<- which(header=="note")

	table(infos[,seg_col])
	table(infos[,gene_col])
	table(infos[,prod_col])
	table(infos[,note_col])

	protNames<- c("[Vv][Pp]1","polyprotein")
	

	infos[,gene_col] <- tolower(infos[,gene_col])
	regAssign	     <- array(0, length( infos[,gene_col] ) )

	for (s in 1:length(protNames)) {

		inds1 <- grep( protNames[s], infos[,gene_col] )
		inds2 <- grep( protNames[s], infos[,prod_col] )
		inds3 <- grep( protNames[s], infos[,note_col] )

		sinds  <- sort(unique( c(inds1, inds2, inds3) ) )
		regAssign[sinds] <- s
	}


	table(infos[which(regAssign==0),gene_col])
	table(infos[which(regAssign==0),prod_col])

	unknown_inds <- which(regAssign == 0)
	unknown_locii <- infos[unknown_inds,1]
	minds		 <- match(unknown_locii, all_locii)
	print(minds)

	#write unprocessed records to separate file
	unprocName <- paste(path,path2,name,"_unprocessed.gb",sep="")
	for (j in 1:length(minds)) {
		write( recs[[minds[j]]], file=unprocName, append=(j>1))
	}



	infos <- cbind(infos, regAssign)
	write.table( infos, file=paste(path,path2,name,"_process_info_with_cds_and_regAssign.txt",sep=""), sep="\t", col.names=TRUE, row.names=FALSE)
	


	# split to regions
	reg_col     <- length(infos[1,])
	regions	<- infos[,reg_col]
	ureg		<- sort(unique(regions))
	pn       <- c("unassign","VP1","polyprotein")
	for (s in 1:length(ureg)) {
		sinds <- which(infos[,reg_col]==ureg[s])
		print( c(ureg[s], length(sinds) ) )

		sname <- paste(path,path2,"fmdv_",pn[s],".gb",sep="")
		sname2<- paste(path,path2,"fmdv_",pn[s],"_nucls.fas",sep="")
		s_loc <- infos[sinds,1]
		minds <- match(s_loc, all_locii)
		all( s_loc == all_locii[minds] )
		
		sname0<- paste(path,path2,"fmdv_",pn[s],"_accn_list.txt",sep="")
		write(s_loc, file=sname0)

		for (j in 1:length(minds) ) {
			write( recs[[minds[j]]], file=sname, append=(j>1) )

			nucls <- getNuclsFromGenBank( recs[[minds[j]]] )
			write( paste(">",nucls[1],sep=""), file=sname2, append=(j>1) )
			write( nucls[2], file=sname2, append=TRUE)
		}
	
		sname3 <- paste(path,path2,"fmdv_",pn[s],"_process_info.txt",sep="")
		write.table( infos[sinds,], file=sname3, sep="\t", col.names=TRUE, row.names=FALSE)
		
	}



# number of sequences per region
#unassign		0 613
#VP1    		1 3344
#polyprotein    	2 1491

###################################################################################
# remove some of the short sequences

	minLen <- c(500,500,500)

############
# VP1 ONLY #
############
	s <- 2
		reg_info <- as.matrix( read.table( paste(path,path2,"fmdv_",pn[s],"_process_info.txt",sep=""), sep="\t", header=TRUE) )
		pos_col  <- grep("position",colnames(reg_info))
		codon_col<- grep("codon_start",colnames(reg_info))
		comp_col <- grep("COMPLEMENT",colnames(reg_info))
		pos      <- reg_info[,pos_col]
		cstart   <- reg_info[,codon_col]
		complement<- reg_info[,comp_col]=="TRUE"

		len1	   <- apply(as.matrix(pos), 1, getEl, ind=1, sep="\\.\\.")
		len1	   <- gsub("<", "", len1)
		len1	   <- gsub("complement\\(", "", len1)
		len1	   <- as.integer(len1)
		len2	   <- apply(as.matrix(pos), 1, getEl, ind=2, sep="\\.\\.")
		len2	   <- gsub(">", "", len2)
		len2     <- gsub(")", "", len2)
		len2	   <- as.integer(len2)

		len      <- len2-len1

		inds     <- which(len >= minLen[s])
		reg_info2<- reg_info[inds,]
		write.table( reg_info2, paste(path,path2,"fmdv_",pn[s],"_process_info_long.txt",sep=""), sep="\t", col.names=TRUE, row.names=FALSE)

		seqs	   <- read.dna( paste(path,path2,"fmdv_",pn[s],"_nucls.fas",sep=""), format="fasta", as.matrix=FALSE)
		taxa	   <- attributes(seqs)$names

		sel_locii <- reg_info2[,1]
		minds	   <- match(sel_locii, taxa)
		all( sel_locii == taxa[minds])

		write.dna( seqs[minds], file=paste(path,path2,"fmdv_",pn[s],"_accns_al0.fas",sep=""), format="fasta", nbcol=-1, colsep="")

		# re-read in sequences, have been partially aligned
		seqs	   <- read.dna( paste(path,path2,"fmdv_",pn[s],"_accns_al2.fas",sep=""), format="fasta", as.matrix=FALSE)
		taxa	   <- attributes(seqs)$names

		# by serotype
		ser_col <- grep("serotype",colnames(reg_info2))
		seros   <- reg_info2[,ser_col]
		seros   <- gsub(" ","-",seros)
		seros   <- gsub("FMD-O","O",seros)
		seros   <- gsub("T-","T",seros)
		seros   <- gsub("a-","a",seros)
		table(seros)

		useros  <- sort(unique(seros))

		for (st in 1:length(useros)) {
			inds 		<- which(seros == useros[st])
			sel_locii 	<- reg_info2[inds,1]
			minds		<- match(sel_locii, taxa)
			all (sel_locii == taxa[minds] )

			sname    <- paste(path,path2,"fmdv_",pn[s],"_sero",useros[st],"_accns_al2.fas",sep="")
			write.dna( seqs[minds], file=sname, format="fasta", nbcol=-1, colsep="")
		}
	}


###################################################################################
# reload aligned sequences (some short ones are now removed) - ** VP1 ONLY **

	s	<- 2
	infos <- read.table( paste(path,path2,"fmdv_",pn[s],"_process_info_long.txt",sep=""), sep="\t", header=TRUE)
	header<- colnames(infos)

	host_col	<- which(header=="host")
	host		<- infos[,host_col]
	table(host)
	
	country_col <- which(header=="country")
	country	<- infos[,country_col]
	country	<- apply(as.matrix(country), 1, getEl, ind=1, sep=":")
	country	<- gsub(" and ","-and-",country)
	country	<- gsub(" ","", country)
	country	<- gsub("'", "", country)

	date_col    <- which(header=="collection_date")
	dateTxt	<- as.matrix(infos[,date_col])
	year		<- apply(as.matrix(dateTxt), 1, getEl, ind=1, sep="-", fromEnd=TRUE)
	year[which(year==1)] <- 0

	str1		<- infos[,which(header=="strain")]
	str2		<- infos[,which(header=="isolate")]
	strain	<- paste(str2)
	inds		<- which(str2 == "-1")
	strain[inds]<- str1[inds]

	inds		<- intersect( which(year == 0), grep("/", strain) )
	inds		<- setdiff(inds, grep("XJ1",strain))
	tempY		<- apply( as.matrix(strain[inds]), 1, getEl, ind=1, fromEnd=TRUE, sep="/")
	tempY		<- gsub("c", "", tempY)
	tempY		<- gsub(")", "", tempY)
	tempY		<- as.integer(tempY)
	inds2		<- which(tempY < 15)
	tempY[inds2]<- 2000+tempY[inds2]
	inds2		<- which(tempY < 100)
	tempY[inds2]<- 1900+tempY[inds2]
	inds2		<- which(!is.finite(tempY))
	tempY[inds2]<- 0
	
	year[inds]	<- tempY
	dateTxt[inds]<- paste(tempY)

	ii		<- which( apply(as.matrix(dateTxt), 1, nchar) > 4  )
	month		<- array(0, length(dateTxt))
	month[ii]	<- apply(as.matrix(dateTxt[ii]), 1, getEl, ind=2, sep="-", fromEnd=TRUE)

	ii		<- which( apply(as.matrix(dateTxt), 1, nchar) > 8  )
	day		<- array(0, length(dateTxt))
	day[ii]	<- apply(as.matrix(dateTxt[ii]), 1, getEl, ind=3, sep="-", fromEnd=TRUE)

	decDates	<- array(0, length(dateTxt))
	monthNames  <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
	mm		<- match(month,monthNames)
	jj		<- which(!is.finite(mm))
	mm[jj]	<- 0
	for (i in 1:length(dateTxt)) {
		if (year[i] > 0) {
			decDates[i] <- calcDecimalDate(as.integer(day[i]), as.integer(mm[i]), as.integer(year[i]))
		}
	}
	decDatesTxt <- gsub(" ","",format(decDates, digits=8))
	decDatesTxt <- gsub("0.0000","0",decDatesTxt,fixed=TRUE)

	

	newNames 	<- paste(host, country, year, strain, dateTxt, decDatesTxt, sep="|")

	newInfos    <- infos
	newInfos    <- cbind(newNames, newInfos)
	colnames(newInfos) <- c("SequenceName",header)

	write.table( newInfos, paste(path,path2,"fmdv_",pn[s],"_process_info_long_with_newNames.txt",sep=""), 
		sep="\t", col.names=TRUE, row.names=FALSE)

	newInfos 	<- as.matrix(newInfos)

	useros 	<- c("A","Asia1","C","O","SAT1","SAT2","SAT3")

	# re-load sequences for renaming (these have been aligned to each other)
		for (s in 1:length(useros)) {
		
		fname <- paste(path,path2,"fmdv_VP1_sero",useros[s],"_accns_al.fas",sep="")
		
		seqs  <- read.dna( fname, format="fasta", as.matrix=FALSE)
		

		taxa  <- attributes(seqs)$names
		minds <- match(taxa, newInfos[,2])
		all( taxa == newInfos[minds,2])

		newSeqs <- seqs
		attributes(newSeqs)$names <- newInfos[minds,1]
		fname2  <- paste(path,path2,"fmdv_VP1_sero",useros[s],"_names.fas",sep="")
		write.dna( newSeqs, file=fname2, format="fasta", nbcol=-1, colsep="")

		print(fname2)

		write.table( newInfos[minds,], file=paste(path,path2,"fmdv_VP1_sero",useros[s],"_names_sequenceInformation.txt",sep=""),
					sep="\t", col.names=TRUE, row.names=FALSE)

		nu	<- length(unique(newInfos[minds,1]))
		if ((nu > 1) & (nu == length(newInfos[minds,1]))) {

			tname	<- paste(path,path2,"fmdv_VP1_sero",useros[s],"_accns_al3_tn93_nj.nwk",sep="")
			tr    <- read.tree(tname)
			taxa  <- attributes(seqs)$names
			minds <- match(taxa, newInfos[,2])
			all( taxa == newInfos[minds,2] )

			newTr <- tr
			newTr$tip.label <- newInfos[minds,1]
			tname2	    <- paste(path,path2,"fmdv_VP1_sero",useros[s],"_names_tn93_nj.nwk",sep="")
			write.tree( newTr, file=tname2 )
			print( tname2 )

		} else {
			print( paste(useros[s],"not unique names") )
		}

	}

	######
	# extract just sequences with date, host and country information
	
	host	   <- gsub("not known","-1",host)
	host	   <- gsub("BHK-21 cells","-1",host)
	host	   <- gsub("host unknown","-1",host)



	table(host[sel_inds])

#		Aepyceros melampus (impala)       African buffalo (Syncerus caffer)                            BHK-21 cells 
#                                      2                                      10                                       2 
#                              blue bull                     Bos grunniens (yak)                                 Bos sp. 
#                                      1                                       5                                       1 
#                             Bos taurus                                  bovine         Bubalus bubalis (water buffalo) 
#                                     32                                     192                                       2 
#                                buffalo                            buffalo calf                                 caprine 
#                                     63                                       4                                       1 
#                                 cattle                                 gazelle                                    goat 
#                                   1201                                       1                                       7 
#                                  goats                            host unknown                                     pig 
#                                      1                                       6                                      42 
#                                   pigs                                 porcine                                   sheep 
#                                      1                                      35                                       6 
#                             Sus scrofa                                   swine Tragelaphus strepsiceros (greater kudu) 
#                                     72                                       1                                       1 
#                          water buffalo                              wildebeest                                     yak 
#                                     10                                       1                                       2

	host <- gsub("Aepyceros melampus (impala)","impala",host,fixed=TRUE)
	host <- gsub("African buffalo (Syncerus caffer)","africanbuffalo",host,fixed=TRUE)
	host <- gsub("blue bull","cattle",host)
	host <- gsub("Bos taurus","cattle",host)
	host <- gsub("buffalo calf","buffalo",host)
	host <- gsub("pigs","porcine",host)
	host <- gsub("swine","porcine",host)
	host <- gsub("Sus scrofa","wildboar",host)
	host <- gsub("Tragelaphus strepsiceros (greater kudu)","antelope_kudu",host,fixed=TRUE)
	host <- gsub("Bubalus bubalis (water buffalo)","waterbuffalo",host,fixed=TRUE)
	host <- gsub("Bos sp.","bovine",host,fixed=TRUE)
	host <- gsub("Bos grunniens (yak)","yak",host,fixed=TRUE)
	host <- gsub("unvaccinated cattle","cattle",host)
	
	newNames 	<- paste(host, country, year, strain, dateTxt, decDatesTxt, sep="|")
	newInfos[,1]<- newNames
	
	sel_inds <- which( (host != "-1") & (decDates > 0) & (country != "-1") )
	write.table( newInfos[sel_inds,], file=paste(path,path2,"fmdv_VP1_selected_sequenceInformation.txt",sep=""), 
				sep="\t", col.names=TRUE, row.names=FALSE)
	



	sel_accns <-  newInfos[sel_inds,2]
	for (s in 1:length(useros)) {
		fname <- paste(path,path2,"fmdv_VP1_sero",useros[s],"_accns_al.fas",sep="")
		seqs  <- read.dna( fname, format="fasta", as.matrix=FALSE)
		taxa  <- attributes(seqs)$names

		
		sel_taxa 	<- intersect(sel_accns, taxa)
		minds1	<- match(sel_taxa, taxa)
		all( sel_taxa == taxa[minds1] )
		newSeqs	<- seqs[minds1]
		
		minds2 	<- match(sel_taxa, newInfos[,2])
		all( sel_taxa == newInfos[minds2,2])
		all( attributes(newSeqs)$names == newInfos[minds2,2] )
		attributes(newSeqs)$names <- newInfos[minds2,1]

		fname2  <- paste(path,path2,"fmdv_VP1_sero",useros[s],"_selected_names.fas",sep="")
		#write.dna( newSeqs, file=fname2, format="fasta", nbcol=-1, colsep="")

		#print( paste(fname2,length(minds2)) )

		#write.table( newInfos[minds2,], file=paste(path,path2,"fmdv_VP1_sero",useros[s],"_selected_names_sequenceInformation.txt",sep=""),
		#			sep="\t", col.names=TRUE, row.names=FALSE)

		print( paste(useros[s],length(taxa),length(minds2)) )
	}

#[1] "A 677 494"
#[1] "Asia1 226 166"
#[1] "C 29 2"
#[1] "O 1203 779"
#[1] "SAT1 38 38"
#[1] "SAT2 225 167"
#[1] "SAT3 4 4"

#[1] "D://slycett//epic_data//BACKGROUND_SEQUENCES//FootAndMouthDiseaseVirus//2014_may_28//fmdv_VP1_seroA_selected_names.fas 494"
#[1] "D://slycett//epic_data//BACKGROUND_SEQUENCES//FootAndMouthDiseaseVirus//2014_may_28//fmdv_VP1_seroAsia1_selected_names.fas 166"
#[1] "D://slycett//epic_data//BACKGROUND_SEQUENCES//FootAndMouthDiseaseVirus//2014_may_28//fmdv_VP1_seroC_selected_names.fas 2"
#[1] "D://slycett//epic_data//BACKGROUND_SEQUENCES//FootAndMouthDiseaseVirus//2014_may_28//fmdv_VP1_seroO_selected_names.fas 779"
#[1] "D://slycett//epic_data//BACKGROUND_SEQUENCES//FootAndMouthDiseaseVirus//2014_may_28//fmdv_VP1_seroSAT1_selected_names.fas 38"
#[1] "D://slycett//epic_data//BACKGROUND_SEQUENCES//FootAndMouthDiseaseVirus//2014_may_28//fmdv_VP1_seroSAT2_selected_names.fas 167"
#[1] "D://slycett//epic_data//BACKGROUND_SEQUENCES//FootAndMouthDiseaseVirus//2014_may_28//fmdv_VP1_seroSAT3_selected_names.fas 4"

	source("Rcode//fourWaySampler.R")

	# re-load these and select one per host per country per year using fourWaySampler
	recomb_name		<- paste(path,path2,"fmdv_VP1_subsampled_host_country_year.fas",sep="")
	
	for (s in 1:length(useros)) {
		seqs 				<- read.dna( paste(path,path2,"fmdv_VP1_sero",useros[s],"_selected_names.fas",sep=""), format="fasta", as.matrix=FALSE)
		taxa 				<- as.matrix( attributes(seqs)$names )
		sel_host 			<- apply(taxa, 1, getEl, ind=1, sep="\\|")
		sel_country 		<- apply(taxa, 1, getEl, ind=2, sep="\\|")
		sel_year			<- apply(taxa, 1, getEl, ind=3, sep="\\|")
		taxa				<- paste(useros[s],taxa,sep="|")
		attributes(seqs)$names 	<- taxa
		chosen_inds <- fourWaySampler( sel_host, sel_country, sel_year, sel_year, maxPer=1 )
		print( paste(useros[s],length(taxa),length(chosen_inds)) )
		
		write.dna(seqs[chosen_inds], file=recomb_name, append=(s>1), format="fasta", nbcol=-1, colsep="")
	}

#	[1] "A 494 153"
#	[1] "Asia1 166 52"
#	[1] "C 2 1"
#	[1] "O 779 190"
#	[1] "SAT1 38 19"
#	[1] "SAT2 167 85"
#	[1] "SAT3 4 2"

###################################################################################
###################################################################################
# Get the Carrillo full genome sequences in the genBank download from the PopSet
# http://www.ncbi.nlm.nih.gov/popset?DbFrom=nuccore&Cmd=Link&LinkName=nuccore_popset&IdsFromResult=46810794
	

	path		   <- "D://slycett//epic_data//BACKGROUND_SEQUENCES//FootAndMouthDiseaseVirus//"
	path2		   <- "carrillo_2005//"
	
	name		   	<- "carrillo_sequences_genBank_full"
	gblines 		<- readLines( paste(path,path2,name,".gb",sep="") )
	recordStarts 	<- grep("LOCUS",gblines)
	nrec			<- length(recordStarts)
	recordEnds		<- c(recordStarts[2:nrec]-1,length(gblines))

	recs			<- vector("list",nrec)
	for (i in 1:nrec) {
		recs[[i]] <- gblines[ recordStarts[i]:recordEnds[i] ]
	}

	nCDS			<- unlist( lapply(recs, numberCDS) )
	jj			<- which(nCDS > 1)

	infos			<- lapply(recs, getGenBankInfo)
	nels			<- length(infos[[1]])
	header		<- rownames(infos[[1]])
	infos 		<- matrix( unlist(infos), nels, length(infos))
	infos			<- t(infos)
	colnames(infos)	<- header
	accns			<- infos[,1]
	seros			<- infos[,which(header=="serotype")]
	country		<- infos[,which(header=="country")]
	country		<- apply(as.matrix(country), 1, getEl, ind=1, sep=":")
	country		<- gsub(" ","", country)
	strain		<- infos[,9]
	
	details		<- read.table( paste(path,path2,"carrillo_table1_details.csv",sep=""), sep=",", header=TRUE)
	

	sname			<- paste(path,path2,"carrillo_2005_nucls.fas",sep="")
	year			<- array(0, length(recs))
	seqName		<- array(0, length(recs))
	for (i in 1:length(recs)) {
		dinds		<- match( accns[i], details[,4])
		res		<- getNuclsFromGenBank(recs[[i]])
		if (accns[i] == res[1]) {
			if (country[i] == -1) {
				country[i] <- details[dinds,2]
			}
			year[i]	<- details[dinds,3]
			seqName[i]	<- paste(accns[i],seros[i],country[i],year[i],strain[i],sep="|")
		
			write( paste(">",seqName[i],sep=""), file=sname, append=(i>1))
			write(res[2], file=sname, append=TRUE)
		} else {
			print(paste("Warning",i,accns[i],res[1],"not match"))
		}
	}

	
###################################################################################
###################################################################################
# Get the Cottam 2006 data (full genomes)
# GenBank Search term = fmdv ukg 2001 complete genome, 44 sequences

	path		   <- "D://slycett//epic_data//BACKGROUND_SEQUENCES//FootAndMouthDiseaseVirus//"
	path2		   <- "cottam_2006//"

	# read GenBank records
	name		   	<- "fmdv_ukg_2001_complete_genome_sequences_genBank_full"
	gblines 		<- readLines( paste(path,path2,name,".gb",sep="") )
	recordStarts 	<- grep("LOCUS",gblines)
	nrec			<- length(recordStarts)
	recordEnds		<- c(recordStarts[2:nrec]-1,length(gblines))

	recs			<- vector("list",nrec)
	for (i in 1:nrec) {
		recs[[i]] <- gblines[ recordStarts[i]:recordEnds[i] ]
	}

	nCDS			<- unlist( lapply(recs, numberCDS) )
	jj			<- which(nCDS > 1)

	infos			<- lapply(recs, getGenBankInfo)
	nels			<- length(infos[[1]])
	header		<- rownames(infos[[1]])
	infos 		<- matrix( unlist(infos), nels, length(infos))
	infos			<- t(infos)
	colnames(infos)	<- header
	accns			<- infos[,1]
	seros			<- infos[,which(header=="serotype")]
	country		<- infos[,which(header=="country")]
	country		<- apply(as.matrix(country), 1, getEl, ind=1, sep=":")
	country		<- gsub(" ","", country)
	strain		<- infos[,9]


	# read details from Cottam 2006 Table 1
	details	   <- read.table( paste(path,path2,"cottam_2006_table_1.csv",sep=""), sep=",", header=TRUE)
	details_accns  <- gsub(" ","",paste(details[,1]))
	details_year   <- as.integer(apply(as.matrix(details[,5]), 1, getEl, ind=3, sep="/"))
	details_month  <- as.integer(apply(as.matrix(details[,5]), 1, getEl, ind=2, sep="/"))
	details_day	   <- as.integer(apply(as.matrix(details[,5]), 1, getEl, ind=1, sep="/"))
	decDates	   <- array(0, length(details_year))
	for (i in 1:length(details_year)) {
		decDates[i] <- calcDecimalDate(details_day[i], details_month[i], details_year[i])
	}
	decDatesTxt	   <- format(decDates, digits=8)

	# make new sequence names - e.g. >DQ404180|Pig|Essex|Farm1|UKG/11/2001|19/02/2001|2001.1342
	details_newNames <- paste(details_accns, details[,4], details[,6], 
						paste("Farm",details[,2],sep=""), gsub("G ","G/",details[,3]), 
						details[,5], decDatesTxt, sep="|")
	details_newNames <- gsub(" ","", details_newNames)

	# find these sequences in the download
	minds		   <- match(details_accns, infos[,1])
	all( details_accns == infos[minds,1] )
	
	# write to fasta file with new names
	fname		   <- paste(path,path2,"cottam_2006_nucls.fas",sep="")
	for (i in 1:length(details_accns)) {
		res <- getNuclsFromGenBank( recs[[minds[i]]] )
		if (res[1] == details_accns[i]) {
			write( paste(">",details_newNames[i],sep=""), file=fname, append=(i>1) )
			write( res[2], file=fname, append=TRUE)
		} else {
			print(paste("Problem with",i,details_accns[i],res[1]) )
		}
	}

	##################
	# read details from Cottam 2008 Table 1
	path2		   <- "cottam_2008//"
	details	   <- read.table( paste(path,path2,"cottam_2008_table_1.csv",sep=""), sep=",", header=TRUE)
	details_accns  <- gsub(" ","",paste(details[,2]))
	details_year   <- as.integer(apply(as.matrix(details[,4]), 1, getEl, ind=3, sep="/"))
	details_month  <- as.integer(apply(as.matrix(details[,4]), 1, getEl, ind=2, sep="/"))
	details_day	   <- as.integer(apply(as.matrix(details[,4]), 1, getEl, ind=1, sep="/"))
	decDates	   <- array(0, length(details_year))
	for (i in 1:length(details_year)) {
		decDates[i] <- calcDecimalDate(details_day[i], details_month[i], details_year[i])
	}
	decDatesTxt	   <- format(decDates, digits=8)

	

	# make new sequence names
	details_newNames1 <- paste(details_accns, details[,3], "UK", 
						paste("Sample",details[,1],sep=""), sep="|")
	details_newNames2 <- paste(details[,4],decDatesTxt,sep="|")

	# find these sequences in the download
	minds		   <- match(details_accns, infos[,1])
	all( details_accns == infos[minds,1] )

	details_newNames  <- paste(details_newNames1, infos[minds,9], details_newNames2, sep="|")
	
	# write to fasta file with new names
	fname		   <- paste(path,path2,"cottam_2008_nucls.fas",sep="")
	for (i in 1:length(details_accns)) {
		res <- getNuclsFromGenBank( recs[[minds[i]]] )
		if (res[1] == details_accns[i]) {
			write( paste(">",details_newNames[i],sep=""), file=fname, append=(i>1) )
			write( res[2], file=fname, append=TRUE)
		} else {
			print(paste("Problem with",i,details_accns[i],res[1]) )
		}
	}


###################################################################################
###################################################################################
# Cottam 2008, the 2007 FMDV outbreak
# genbank search term (UKG 2007) AND "Foot-and-mouth disease virus - type O"[porgn:__txid12118] 

	path		   <- "D://slycett//epic_data//BACKGROUND_SEQUENCES//FootAndMouthDiseaseVirus//"
	path2		   <- "cottam_fmdv_2007_outbreak//"

	# read GenBank records
	name		   	<- "fmdv_2007_sequences_genBank_full"
	gblines 		<- readLines( paste(path,path2,name,".gb",sep="") )
	recordStarts 	<- grep("LOCUS",gblines)
	nrec			<- length(recordStarts)
	recordEnds		<- c(recordStarts[2:nrec]-1,length(gblines))

	recs			<- vector("list",nrec)
	for (i in 1:nrec) {
		recs[[i]] <- gblines[ recordStarts[i]:recordEnds[i] ]
	}

	nCDS			<- unlist( lapply(recs, numberCDS) )
	jj			<- which(nCDS > 1)

	infos			<- lapply(recs, getGenBankInfo)
	nels			<- length(infos[[1]])
	header		<- rownames(infos[[1]])
	infos 		<- matrix( unlist(infos), nels, length(infos))
	infos			<- t(infos)
	colnames(infos)	<- header
	accns			<- infos[,1]
	seros			<- infos[,which(header=="serotype")]
	country		<- infos[,which(header=="country")]
	country		<- apply(as.matrix(country), 1, getEl, ind=1, sep=":")
	country		<- gsub(" ","", country)
	strain		<- infos[,8]


	# read details from Cottam 2008 Table 1
	details	   <- read.table( paste(path,path2,"cottam_2008_table_1.txt",sep=""), sep=",", header=TRUE)
	details_accns  <- gsub(" ","",paste(details[,4]))
	details_year   <- as.integer(apply(as.matrix(details[,2]), 1, getEl, ind=3, sep="/"))
	details_month  <- as.integer(apply(as.matrix(details[,2]), 1, getEl, ind=2, sep="/"))
	details_day	   <- as.integer(apply(as.matrix(details[,2]), 1, getEl, ind=1, sep="/"))
	decDates	   <- array(0, length(details_year))
	for (i in 1:length(details_year)) {
		decDates[i] <- calcDecimalDate(details_day[i], details_month[i], details_year[i])
	}
	decDatesTxt	   <- format(decDates, digits=8)

	minds		  <- match(details_accns, infos[,1])
	host		  <- paste(details[,5])
	host		  <- gsub("BOVINE","Cattle",host)
	host		  <- gsub("OVINE","Sheep",host)

	# make new sequence names
	details_newNames <- paste(details_accns, host, "UK", 
						paste("Holding-",details[,1],sep=""), details[,3], details[,2],decDatesTxt,sep="|")


	# write to fasta file with new names
	fname		   <- paste(path,path2,"fmdv_2007_nucls.fas",sep="")
	for (i in 1:length(details_accns)) {
		res <- getNuclsFromGenBank( recs[[minds[i]]] )
		if (res[1] == details_accns[i]) {
			write( paste(">",details_newNames[i],sep=""), file=fname, append=(i>1) )
			write( res[2], file=fname, append=TRUE)
		} else {
			print(paste("Problem with",i,details_accns[i],res[1]) )
		}
	}

	
