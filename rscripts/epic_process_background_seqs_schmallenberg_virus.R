# script to process Genbank download of Schmallenberg Virus
# for EPIC
# S. J. Lycett
# 28 May 2014

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
# process the data

	path	<- "D://slycett//epic_data//BACKGROUND_SEQUENCES//"
	path2 <- "SchmallenbergVirus//2014_may_28//"
	name  <- "schmallenberg_sequence_genBank_full"
	
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

	# names of the proteins / segments
	# S, M, L
	# N = nucleocapsid is on S segment
	# RdRp is on L segment

	segNames <- c("S","M","L")
	protNames<- c("N","M","RdRp")

	infos[,gene_col] <- toupper(infos[,gene_col])

	segAssign	     <- array(0, length( infos[,gene_col] ) )
	for (s in 1:length(segNames)) {
		inds1	<- grep( segNames[s], infos[,gene_col] )
		inds2 <- grep( segNames[s], infos[,prod_col] )

		inds3 <- grep( protNames[s], infos[,gene_col] )
		inds4 <- grep( protNames[s], infos[,prod_col] )

		inds5 <- grep( segNames[s],  infos[,seg_col] )
		inds6 <- grep( protNames[s], infos[,seg_col] )

		inds7	<- grep( segNames[s], infos[,note_col] )
		inds8 <- grep( protNames[s], infos[,note_col] )

		inds9	<- which( infos[,seg_col] == s )

		sinds <- sort( unique( c(inds1, inds2, inds3, inds4, inds5, inds6, inds7, inds8, inds9) ) )
		segAssign[sinds] <- s
	}


	infos[which(segAssign==0),gene_col]
	infos[which(segAssign==0),prod_col]

	unknown_inds <- which(segAssign == 0)
	unknown_locii <- infos[unknown_inds,1]
	minds		 <- match(unknown_locii, all_locii)
	print(minds)

	# write unprocessed records to separate file
	#unprocName <- paste(path,path2,name,"_unprocessed.gb",sep="")
	#for (j in 1:length(minds)) {
	#	write( recs[[minds[j]]], file=unprocName, append=(j>1))
	#}



	infos <- cbind(infos, segAssign)
	write.table( infos, file=paste(path,path2,name,"_process_info_with_cds_and_segAssign.txt",sep=""), sep="\t", col.names=TRUE, row.names=FALSE)
	


	# split to segment
	seg_col     <- length(infos[1,])
	for (s in 1:length(segNames)) {
		sinds <- which(infos[,seg_col]==s)
		print( c(s, length(sinds) ) )

		sname <- paste(path,path2,"schmallenberg_seg",segNames[s],".gb",sep="")
		sname2<- paste(path,path2,"schmallenberg_seg",segNames[s],"_nucls.fas",sep="")
		s_loc <- infos[sinds,1]
		minds <- match(s_loc, all_locii)
		for (j in 1:length(minds) ) {
			write( recs[[minds[j]]], file=sname, append=(j>1) )

			nucls <- getNuclsFromGenBank( recs[[minds[j]]] )
			write( paste(">",nucls[1],sep=""), file=sname2, append=(j>1) )
			write( nucls[2], file=sname2, append=TRUE)
		}
	
		sname3 <- paste(path,path2,"schmallenberg_seg",segNames[s],"_process_info.txt",sep="")
		write.table( infos[sinds,], file=sname3, sep="\t", col.names=TRUE, row.names=FALSE)
		
	}

# number of sequences per segment
# S  1 35
# M  2 39
# L  3 15


###################################################################################
# reload aligned sequences (some short ones are now removed)

	infos <- read.table( paste(path,path2,name,"_process_info_with_cds_and_segAssign.txt",sep=""), sep="\t", header=TRUE)
	header<- colnames(infos)

	

	host_col	<- which(header=="host")
	host		<- infos[,host_col]

	#> table(host)
	#host
	#                   -1            Bos taurus                  goat            Ovis aries Ovis orientalis aries                 sheep 
	#                    3                    17                     2                    11                    36                    20

	host		<- gsub("Bos taurus","Cattle",host)
	host		<- gsub("goat","Goat",host)
	host		<- gsub("Ovis aries","Sheep",host)
	host		<- gsub("Ovis orientalis aries","Sheep",host)
	host		<- gsub("sheep","Sheep",host)
	host[ which(host==-1) ] <- "Lab"
	table(host)
	
	country_col <- which(header=="country")
	country	<- infos[,country_col]
	country	<- apply(as.matrix(country), 1, getEl, ind=1, sep=":")
	country	<- gsub(" and ","-and-",country)
	country	<- gsub(" ","", country)

	date_col    <- which(header=="collection_date")
	dateTxt	<- infos[,date_col]
	year		<- apply(as.matrix(dateTxt), 1, getEl, ind=1, sep="-", fromEnd=TRUE)
	year[which(year==1)] <- 0

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

	str1		<- infos[,which(header=="strain")]
	str2		<- infos[,which(header=="isolate")]
	strain	<- paste(str2)
	inds		<- which(str2 == "-1")
	strain[inds]<- str1[inds]
	newNames 	<- paste(host, country, year, strain, dateTxt, decDatesTxt, sep="|")

	newInfos    <- infos
	newInfos    <- cbind(newNames, newInfos)
	colnames(newInfos) <- c("SequenceName",header)

	for (s in 1:length(segNames)) {
		seqs <- read.dna( paste(path,path2,"schmallenberg_",segNames[s],"_accns.fas",sep=""), format="fasta", as.matrix=FALSE)
		taxa <- attributes(seqs)$names
		minds<- match(taxa, infos[,1])
		all( taxa == infos[minds,1] )

		newSeqs 			  <- seqs
		attributes(newSeqs)$names <- newNames[minds]

		write.table( newInfos[minds,], file=paste(path,path2,"schmallenberg_",segNames[s],"_names_info.txt",sep=""), 
					sep="\t", col.names=TRUE, row.names=FALSE)

		write.dna( newSeqs, file=paste(path,path2,"schmallenberg_",segNames[s],"_names.fas",sep=""), format="fasta", nbcol=-1, colsep="")
	}

	# load in trees and rename also
	for (s in 1:length(segNames)) {
		tr 	<- read.tree( paste(path,path2,"schmallenberg_",segNames[s],"_accns_tn93_nj.nwk",sep="") )
		tr	<- ladderize(tr)
		minds <- match(tr$tip.label, newInfos[,2])
		all( tr$tip.label == newInfos[minds,2] )
		tr$tip.label <- paste(newInfos[minds,1])
		write.tree( tr, file=paste(path,path2,"schmallenberg_",segNames[s],"_names_tn93_nj.nwk",sep="") )
		rm(tr)

		tr 	<- read.tree( paste(path,path2,"schmallenberg_",segNames[s],"_accns_gtr_ml_bs100.nwk",sep="") )
		tr	<- ladderize(tr)
		minds <- match(tr$tip.label, newInfos[,2])
		all( tr$tip.label == newInfos[minds,2] )
		tr$tip.label <- paste(newInfos[minds,1])
		write.tree( tr, file=paste(path,path2,"schmallenberg_",segNames[s],"_names_gtr_ml_bs100.nwk",sep="") )
		rm(tr)
	}
	