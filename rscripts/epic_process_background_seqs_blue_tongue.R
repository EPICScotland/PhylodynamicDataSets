# script to process Genbank download of Blue Tongue Virus
# for EPIC
# S. J. Lycett
# 19 May 2014

library(ape)
source("Rcode//getEl.R")
source("Rcode//calcDecimalDate.R")

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
	path2 <- "BlueTongueVirus//2014_may_19//"
	name  <- "blueTongue_sequence_genBank_full"
	
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
	

	# names of the proteins / segments (from AfricanHorseSickness)
	# L1	3500	VP1	core	segment 1
	# L2	3229	VP2	outer shell	segment 2
	# L3	2792	VP3	core	segment 3
	# M4	1978	VP4	core	segment 4
	# M5	1751	NS1	non-structural	segment 5
	# M6	1566	VP5	outer shell	segment 6
	# S7	1179	VP7	core	segment 7
	# S8	116	NS2	non-structural	segment 8
	# S9	1100	VP6	core	segment 9
	# S10	756	NS3/3A	non-structural	segment 10

	segNames <- c("L1","L2","L3","M4","M5","M6","S7","S8","S9","S10")
	protNames<- c("VP1","VP2","VP3","VP4","NS1","VP5","VP7","NS2","VP6","NS3")

	infos[,prod_col] <- gsub("viral protein","VP",infos[,prod_col])
	infos[,prod_col] <- gsub("nonstructural protein","NS",infos[,prod_col])
	infos[,prod_col] <- gsub("non-structural protein","NS",infos[,prod_col], fixed=TRUE)
	infos[,prod_col] <- gsub("NS ","NS", infos[,prod_col])

	infos[,gene_col] <- toupper(infos[,gene_col])

	jj		     <- which(infos[,gene_col] == "S1")
	if (length(jj) > 0) {
		infos[jj,gene_col] <- "S10"
	}

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
	unprocName <- paste(path,path2,name,"_unprocessed.gb",sep="")
	for (j in 1:length(minds)) {
		write( recs[[minds[j]]], file=unprocName, append=(j>1))
	}



	infos <- cbind(infos, segAssign)
	write.table( infos, file=paste(path,path2,name,"_process_info_with_cds_and_segAssign.txt",sep=""), sep="\t", col.names=TRUE, row.names=FALSE)
	


	# split to segment
	seg_col     <- length(infos[1,])
	for (s in 1:10) {
		sinds <- which(infos[,seg_col]==s)
		print( c(s, length(sinds) ) )

		sname <- paste(path,path2,"blueTongue_seg",s,".gb",sep="")
		sname2<- paste(path,path2,"blueTongue_seg",s,"_nucls.fas",sep="")
		s_loc <- infos[sinds,1]
		minds <- match(s_loc, all_locii)
		for (j in 1:length(minds) ) {
			write( recs[[minds[j]]], file=sname, append=(j>1) )

			nucls <- getNuclsFromGenBank( recs[[minds[j]]] )
			write( paste(">",nucls[1],sep=""), file=sname2, append=(j>1) )
			write( nucls[2], file=sname2, append=TRUE)
		}
	
		sname3 <- paste(path,path2,"blueTongue_seg",s,"_process_info.txt",sep="")
		write.table( infos[sinds,], file=sname3, sep="\t", col.names=TRUE, row.names=FALSE)
		
	}





	
