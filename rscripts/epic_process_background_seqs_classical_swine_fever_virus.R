# script to process Genbank download of Classical Swine Fever Virus
# for EPIC
# S. J. Lycett
# 12 May 2014

library(ape)
source("Rcode//getEl.R")

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
		cds_start[inds[1]]
		cds_end <- length(ft)
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

getGenBankInfo	<- function( rec, names1 =c("host","country","collection_date","isolate","strain","note"),
						names2=c("gene","product","codon_start"), product="polyprotein" ) {
	loc	<- getLocus( rec )
	ft 	<- getFeaturesTable( rec )
	ncds  <- numberCDS( ft )
	if (ncds > 1) {
		cds <- getCDSProduct(ft, product=product)
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
	return ( ans )
}

#########################################################################################################################

# process full genomes with full genbank info
	path	<- "D://slycett//epic_data//BACKGROUND_SEQUENCES//"
	path2 <- "ClassicalSwineFever//2014_may_12//"
	name  <- "csfv_complete_genome_genBank_full_61"

	gblines <- readLines( paste(path,path2,name,".gb",sep="") )

	recordStarts 	<- grep("LOCUS",gblines)
	nrec			<- length(recordStarts)
	recordEnds		<- c(recordStarts[2:nrec]-1,length(gblines))

	recs			<- vector("list",nrec)
	for (i in 1:nrec) {
		recs[[i]] <- gblines[ recordStarts[i]:recordEnds[i] ]
	}

	infos			<- lapply(recs, getGenBankInfo)
	nels			<- length(infos[[1]])
	header		<- rownames(infos[[1]])
	infos 		<- matrix( unlist(infos), nels, length(infos))
	infos			<- t(infos)
	colnames(infos)	<- header

	# write all info to file
	write.table( infos, file=paste(path,path2,name,"_process_info.txt",sep=""), sep="\t", col.names=TRUE, row.names=FALSE)

	

