# script to process Genbank download of African Swine Fever Virus
# for EPIC
# S. J. Lycett
# 6 May 2014

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

getGenBankInfo	<- function( rec, names1 =c("host","country","collection_date","isolate","strain"),
						names2=c("gene","product","codon_start"), product="[Pp]72" ) {
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
# process the P72 data

	path	<- "D://slycett//epic_data//BACKGROUND_SEQUENCES//"
	path2 <- "AfricanSwineFeverVirus//2014_may_6//"
	name  <- "asfv_p72_genBank_full"
	
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

	# select the P72
	temp1		<- grep("B646L",infos[,9])
	temp2		<- grep("B646L",infos[,10])
	temp3		<- grep("[Pp]72", infos[,9])
	temp4		<- grep("[Pp]72", infos[,10])
	selInds	<- sort(unique(c(temp1,temp2,temp3,temp4)))
	selinfo	<- infos[selInds,]

	# write selected info to file
	write.table( selinfo, file=paste(path,path2,name,"_selected_p72_process_info.txt",sep=""), sep="\t", col.names=TRUE, row.names=FALSE)

	# split to short, medium, long for processing
	seqlens	<- as.integer(selinfo[,2])
	short_inds  <- which(seqlens < 1000)
	med_inds	<- which((seqlens >= 1000) & (seqlens < 10000))
	long_inds	<- which(seqlens >= 10000)

	write.table( selinfo[short_inds,], file=paste(path,path2,name,"_selected_p72_process_info_short.txt",sep=""), sep="\t", col.names=TRUE, row.names=FALSE)
	write.table( selinfo[med_inds,],  file=paste(path,path2,name,"_selected_p72_process_info_med.txt",sep=""), sep="\t", col.names=TRUE, row.names=FALSE)
	write.table( selinfo[long_inds,], file=paste(path,path2,name,"_selected_p72_process_info_long.txt",sep=""), sep="\t", col.names=TRUE, row.names=FALSE)

	short_gb_inds <- match(selinfo[short_inds,1], infos[,1])
	outname	  <- paste(path,path2,name,"_selected_p72_short.gb",sep="")
	for (i in 1:length(short_gb_inds)) {
		write(recs[[short_gb_inds[i]]], file=outname, append=(i>1))
	}

	med_gb_inds <- match(selinfo[med_inds,1], infos[,1])
	outname	  <- paste(path,path2,name,"_selected_p72_med.gb",sep="")
	for (i in 1:length(med_gb_inds)) {
		write(recs[[med_gb_inds[i]]], file=outname, append=(i>1))
	}

	long_gb_inds <- match(selinfo[long_inds,1], infos[,1])
	outname	  <- paste(path,path2,name,"_selected_p72_long.gb",sep="")
	for (i in 1:length(long_gb_inds)) {
		write(recs[[long_gb_inds[i]]], file=outname, append=(i>1))
	}

#########################################################################################################################
# find information about studies from papers

source("Rcode/getEl.R")
source("Rcode/calcDecimalDate.R")

	short_info <- read.table( paste(path,path2,name,"_selected_p72_process_info_short.txt",sep=""),
						header=TRUE)

	cn	   		<- colnames(short_info)
	host_col 		<- which(cn=="host")
	country_col 	<- which(cn=="country")
	date_col		<- which(cn=="collection_date")
	is_col		<- which(cn=="isolate")
	str_col		<- which(cn=="strain")

	short_info		<- as.matrix(short_info)
	colnames(short_info) <- cn

	jj					<- grep("Congo",short_info[,country_col])
	short_info[jj,country_col]	<- "DemocraticRepublicCongo"

	jj					<- grep("Ivoire", short_info[,country_col])
	short_info[jj,country_col]	<- "CotedIvoire"

	els					<- apply(as.matrix(short_info[,country_col]), 1, getEl, ind=1, sep=":")
	els					<- gsub(" ","", els)
	short_info[,country_col]	<- els

	
# Lubisi
	lubisi_data <- readLines( paste(path,path2,"lubisi_table1.txt",sep="") )
	lubisi_accn <- apply(as.matrix(lubisi_data), 1, getEl, ind=1, fromEnd=TRUE, sep=" ")
	lubisi_host <- apply(as.matrix(lubisi_data), 1, getEl, ind=2, fromEnd=TRUE, sep=" ")
	lubisi_year <- apply(as.matrix(lubisi_data), 1, getEl, ind=3, fromEnd=TRUE, sep=" ")
	lubisi_strain <- apply(as.matrix(lubisi_data), 1, getEl, ind=1, sep=" ")
	lubisi_country<- apply(as.matrix(lubisi_data), 1, getEl, ind=2, sep=" ")

	# correct host names
	jj 			<- which(lubisi_host=="Phaecochoerus")
	lubisi_host[jj]	<- "Phacochoerus aethiopicus"
	jj			<- which(lubisi_host=="SusScrofa")
	lubisi_host[jj]	<- "Sus scrofa"

	
	minds	<- match(lubisi_accn, short_info[,1])
	all( short_info[minds,1]==lubisi_accn )

	# update host column
	all( short_info[minds,host_col] == -1 )		
	short_info[minds,host_col] <- lubisi_host

	# update year column
	all( short_info[minds,date_col] == -1 )
	short_info[minds,date_col] <- lubisi_year
	
	# update strain column, note that isolate is populated from genbank
	all( short_info[minds,str_col] == -1 )
	short_info[minds,str_col] <- lubisi_strain

	# correct country
	jj <- which( short_info[minds,country_col] != lubisi_country )
	short_info[minds[jj],]
	lubisi_country[jj]
	short_info[minds[jj],country_col]
	short_info[minds[jj],country_col] <- "DemocraticRepublicCongo"

# Boshov
	boshov_data <- readLines(paste(path,path2,"boshov_table1.csv",sep="") )
	boshov_data <- boshov_data[4:length(boshov_data)]
	boshov_strain <- apply(as.matrix(boshov_data), 1, getEl, ind=1, sep=",")
	boshov_country<- apply(as.matrix(boshov_data), 1, getEl, ind=2, sep=",")
	boshov_country<- apply(as.matrix(boshov_country), 1, getEl, ind=1, fromEnd=TRUE, sep="<comma> ")
	boshov_country<- gsub(" ", "", boshov_country)
	boshov_year   <- apply(as.matrix(boshov_data), 1, getEl, ind=3, sep=",")
	boshov_accn	  <- apply(as.matrix(boshov_data), 1, getEl, ind=4, sep=",")

	binds		  <- match(boshov_accn, short_info[,1])

	# didnt find one entry
	jj		  <- which(is.finite(binds))

	# also remove Lubisi entries, because they are done already
	jj		  <- setdiff(jj, grep("Lubi",boshov_data) )
	all( boshov_accn[jj] == short_info[binds[jj],1] )

	# update date
	all( short_info[binds[jj],date_col] == -1 )
	short_info[binds[jj],date_col] <- boshov_year[jj]

	# countries are OK
	all( short_info[binds[jj],country_col] == boshov_country[jj] )

	# update strain
	all( short_info[binds[jj],str_col] == -1 )
	short_info[binds[jj],str_col] <- boshov_strain[jj]

# write updated information to file
	write.table(short_info, paste(path,path2,name,"_selected_p72_process_info_short_updated.txt",sep=""),
				col.names=TRUE, row.names=FALSE)

# some manual corrections to dates in file
	short_info <- read.table( paste(path,path2,name,"_selected_p72_process_info_short_updated2.txt",sep=""),
				header=TRUE)
	cn	     <- colnames(short_info)
	short_info <- as.matrix(short_info)


	short_is 		<- short_info[,is_col]
	short_str		<- short_info[,str_col]
	ii	   		<- which((short_is != "-1") & (short_str == "-1"))
	short_str[ii] 	<- short_is[ii]
	short_str		<- gsub("\\.","/",short_str)
	short_str		<- gsub(" ","_",short_str)

# these Nigerian sequences from Owolodun 2010 only have hapoltypes not strain names
# and it is not possible to exactly match them with the strain, host or date date from the paper
	ii			<- which(short_str == -1)
	find_accns		<- short_info[ii,1]
	for ( i in 1:length(find_accns) ) {
		j 	<- grep(find_accns[i], gblines[recordStarts])
		rec 	<- gblines[recordStarts[j]:recordEnds[j]]
		ft	<- getFeaturesTable(rec)
		r_str <- getInfoFromFeaturesTable(ft, name="haplotype")
		short_str[ii[i]] <- r_str
	}

	short_info[,str_col] <- short_str

	ii	     		<- which(short_info[,date_col] == -1)
	ii			<- intersect(ii, grep("/",short_info[,str_col]))
	temp_str		<- short_info[ii,str_col]
	temp_last		<- as.integer(apply(as.matrix(temp_str), 1, getEl, ind=1, sep="/", fromEnd=TRUE))
	jj			<- which(temp_last < 15)
	temp_last[jj]	<- 2000 + temp_last[jj]
	jj			<- which(temp_last < 100)
	temp_last[jj]	<- 1900 + temp_last[jj]
	jj			<- which(temp_last < 1000)
	temp_last[jj]	<- NA
	jj			<- which(!is.finite(temp_last))
	temp_last[jj]	<- -1
	short_info[ii,date_col] <- temp_last
	short_info[ii,]

# corrections
	jj				<- which(short_info[,is_col]=="NUR 90/1")
	short_info[jj,date_col] <- 1990

	jj				<- which(short_info[,is_col]=="BEN 97/4")
	short_info[jj,date_col] <- 1997

	jj				<- which(short_info[,is_col]=="KAB 94/1")
	short_info[jj,date_col] <- 1994

	write.table(short_info, paste(path,path2,name,"_selected_p72_process_info_short_updated3.txt",sep=""),
				col.names=TRUE, row.names=FALSE)

###
# more corrections, using Bastos 2003

	short_info <- read.table( paste(path,path2,name,"_selected_p72_process_info_short_updated3.txt",sep=""),
				header=TRUE)
	short_info <- as.matrix(short_info)

	print( short_info[ which(short_info[,date_col]==-1), str_col] )

# "Ghana"        "NIG-2"        "IC/576"       "VIR_Uganda"   "NH/P68"       
# "Georgia_2007" 
# "CVR Tet-36"   "CVR Tet-29"   "CVR Tet-27"  
# "CVR Tet-21"   "CVR Tet-20"   
# "Togo"         "Hinde_II"     "Dezda"        "ZAR85"        "Dom_Rep"      "TAN/Kwh12"    "Tengani" 

	jj				<- which(short_info[,str_col]=="Ghana")
	short_info[jj,date_col] <- 2000

	jj				<- which(short_info[,str_col]=="NIG-2")
	short_info[jj,date_col] <- 1998

	jj				<- which(short_info[,str_col]=="Togo")
	short_info[jj,str_col]  <- "Togo/98"
	short_info[jj,date_col]	<- 1998

	jj				<- which(short_info[,str_col]=="IC/576")
	short_info[jj,date_col] <- 1996

	jj				<- which(short_info[,str_col]=="Hinde_II")
	short_info[jj,date_col]	<- 1959	

	jj				<- which(short_info[,str_col]=="Dezda")
	short_info[jj,date_col] <- 1986

	jj				<- which(short_info[,str_col]=="ZAR85")
	short_info[jj,date_col] <- 1985

	jj				<- which(short_info[,str_col]=="Dom_Rep")
	short_info[jj,str_col]  <- "DomRep/79"
	short_info[jj,date_col] <- 1979

	jj				<- which(short_info[,str_col]=="TAN/Kwh12")
	short_info[jj,date_col] <- 1968
	short_info[jj,host_col] <- "warthog Phacochoerus spp."


	jj				<- which(short_info[,str_col]=="Tengani")
	short_info[jj,str_col]  <- "Tengani/60"
	short_info[jj,date_col] <- 1960
	short_info[jj,host_col] <- "warthog Phacochoerus spp."


	jj				<- which(short_info[,str_col]=="Georgia_2007")
	short_info[jj,date_col] <- 2007

	jj				<- which(short_info[,str_col]=="RSA/1/99/W")
	short_info[jj,country_col] <- "SouthAfrica"

# try to get more date details
	ii				<- which(short_info[,date_col] == -1)
	find_accns			<- short_info[ii,1]
	for (i in 1:length(find_accns)) {
		j 	<- grep(find_accns[i], gblines[recordStarts])
		rec 	<- gblines[recordStarts[j]:recordEnds[j]]
		ft	<- getFeaturesTable(rec)
		#print(ft)
	}


# try to get more host details
	ii				<- which(short_info[,host_col] == -1)
	find_accns			<- short_info[ii,1]
	for (i in 1:length(find_accns)) {
		j 	<- grep(find_accns[i], gblines[recordStarts])
		rec 	<- gblines[recordStarts[j]:recordEnds[j]]
		#ft	<- getFeaturesTable(rec)
	}
	print( short_info[ii,str_col] )

	jj				<- grep("RSA/1/99", short_info[,str_col])
	short_info[jj,host_col] <- "warthog Phacochoerus spp."

	jj				<- grep("RSA/W/1/99", short_info[,str_col])
	short_info[jj,host_col] <- "warthog Phacochoerus spp."

	jj				<- which(short_info[,str_col]=="Ourt/88/1")
	short_info[jj,host_col] <- "Ornithodoros ticks"

	jj				<- which(short_info[,str_col]=="KAV89/1")
	short_info[jj,host_col] <- "Ornithodoros ticks"

	jj				<- grep("VICT_90/11", short_info[,str_col])
	short_info[jj,date_col] <- 1990
	short_info[jj,host_col] <- "Ornithodoros ticks"

	jj				<- grep("VICT/90/1", short_info[,str_col])
	short_info[jj,host_col] <- "Ornithodoros ticks"


	ii				<- which(short_info[,host_col] == -1)
	print( short_info[ii,str_col] )


	write.table(short_info, paste(path,path2,name,"_selected_p72_process_info_short_updated4.txt",sep=""),
				col.names=TRUE, row.names=FALSE)

	short_date		<- short_info[,date_col]

	short_year		<- apply(as.matrix(short_date), 1, getEl,ind=1, fromEnd=TRUE, sep="-")
	short_year[which(short_year == "1")] <- 0
	short_year[which(short_year == "NK")] <- 0
	short_year[grep("<",short_year)] <- 1959 
	short_year		<- as.integer(short_year)

	short_month		<- array(0, length(short_date))
	ii			<- setdiff(grep("-",short_date),which(short_year==0))
	short_month[ii]	<- apply(as.matrix(short_date[ii]), 1, getEl, ind=2, fromEnd=TRUE, sep="-")

	mm			<- match( short_month[ii], c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec") )
	short_month[ii]   <- mm
	short_month		<- as.integer(short_month)

	ii			<- setdiff(grep("-",short_date),which(short_month==0))
	short_day		<- array(0, length(short_date))
	short_day[ii]	<- apply(as.matrix(short_date[ii]), 1, getEl, ind=1, sep="-")
	ii			<- which(apply(as.matrix(short_day), 1, nchar) != 2)
	short_day[ii]	<- 0
	short_day		<- as.integer(short_day)

	short_dd		<- array(0,length(short_date))
	for (i in 1:length(short_date)) {
		short_dd[i] <- calcDecimalDate(short_day[i], short_month[i], short_year[i])
	}

	new_shortNames	<- paste(short_info[,1], short_info[,country_col], short_year, short_info[,str_col], short_date, short_dd, sep="|")

	write.table(new_shortNames, paste(path,path2,name,"_selected_p72_short_newNames.txt",sep=""),
				col.names=FALSE, row.names=FALSE)


# medium length sequences
	med_info <- read.table( paste(path,path2,name,"_selected_p72_process_info_med.txt",sep=""), header=TRUE)
	med_info <- as.matrix(med_info)

	jj					<- grep("Congo",med_info[,country_col])
	med_info[jj,country_col]	<- "DemocraticRepublicCongo"

	jj					<- grep("Ivoire", med_info[,country_col])
	med_info[jj,country_col]	<- "CotedIvoire"

	els					<- apply(as.matrix(med_info[,country_col]), 1, getEl, ind=1, sep=":")
	els					<- gsub(" ","", els)
	med_info[,country_col]		<- els


	all(med_info[,str_col] == -1)

	med_is 		<- med_info[,is_col]
	med_str		<- med_info[,str_col]
	ii	   		<- which((med_is != "-1") & (med_str == "-1"))
	med_str[ii] 	<- med_is[ii]
	med_str		<- gsub("\\.","/",med_str)
	med_str		<- gsub(" ","_",med_str)
	med_info[,str_col]<- med_str

	ii	   <- which(med_info[,date_col] == -1)
	find_accns<- med_info[ii,1]
	for ( i in 1:length(find_accns) ) {
		j 	<- grep(find_accns[i], gblines[recordStarts])
		rec 	<- gblines[recordStarts[j]:recordEnds[j]]
		ft	<- getFeaturesTable(rec)

		if (med_info[ii[i], host_col] == -1) {
			r_host<- getInfoFromFeaturesTable(ft, name="isolation_source")
			med_info[ii[i], host_col] <- r_host
			print(r_host)
		}

		if (med_info[ii[i],date_col] == -1) {
			r_nn <- getInfoFromFeaturesTable(ft, name="note")
			r_y  <- gsub("isolated in ","",r_nn)
			med_info[ii[i],date_col] <- r_y
			print(r_nn)
		}
	}
	
#	correction
	ii	<- which(med_info[,is_col]=="BEN/1/97")
	med_info[ii,date_col] <- 1997

	write.table( med_info, file=paste(path,path2,name,"_selected_p72_process_info_med_updated3.txt",sep=""),
				col.names=TRUE, row.names=FALSE)

	med_date		<- med_info[,date_col]

	med_year		<- apply(as.matrix(med_date), 1, getEl,ind=1, fromEnd=TRUE, sep="-")
	med_year[which(med_year == "1")] <- 0
	med_year[which(med_year == "NK")] <- 0
	med_year[grep("<",med_year)] <- 1959 
	med_year		<- as.integer(med_year)

	med_month		<- array(0, length(med_date))
	ii			<- setdiff(grep("-",med_date),which(med_year==0))
	med_month[ii]	<- apply(as.matrix(med_date[ii]), 1, getEl, ind=2, fromEnd=TRUE, sep="-")

	mm			<- match( med_month[ii], c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec") )
	med_month[ii]   <- mm
	med_month		<- as.integer(med_month)

	ii			<- setdiff(grep("-",med_date),which(med_month==0))
	med_day		<- array(0, length(med_date))
	med_day[ii]	<- apply(as.matrix(med_date[ii]), 1, getEl, ind=1, sep="-")
	ii			<- which(apply(as.matrix(med_day), 1, nchar) != 2)
	med_day[ii]	<- 0
	med_day		<- as.integer(med_day)

	med_dd		<- array(0,length(med_date))
	for (i in 1:length(med_date)) {
		med_dd[i] <- calcDecimalDate(med_day[i], med_month[i], med_year[i])
	}

	new_medNames	<- paste(med_info[,1], med_info[,country_col], med_year, med_info[,str_col], med_date, med_dd, sep="|")

	write.table(new_medNames, paste(path,path2,name,"_selected_p72_med_newNames.txt",sep=""),
				col.names=FALSE, row.names=FALSE)


# long length sequences

	long_info <- read.table( paste(path,path2,name,"_selected_p72_process_info_long_updated.txt",sep=""), header=TRUE)
	long_date		<- long_info[,date_col]
	long_year		<- as.integer(long_date)
	long_dd		<- array(0,length(long_date))
	for (i in 1:length(long_date)) {
		long_dd[i] <- calcDecimalDate(0, 0, long_year[i])
	}
	new_longNames	<- paste(long_info[,1], long_info[,country_col], long_year, long_info[,str_col], long_date, long_dd, sep="|")

	write.table(new_longNames, paste(path,path2,name,"_selected_p72_long_newNames.txt",sep=""),
				col.names=FALSE, row.names=FALSE)



################################################################

source("Rcode/matrixTranslate.R")
source("Rcode/fourWaySampler.R")

# concatenate all short, med and long information
	final_infos <- read.table( paste(path,path2,"asfv_p72_genBank_full_selected_p72_process_info_all_final.txt",sep=""), header=TRUE)
	final_names <- readLines(  paste(path,path2,"asfv_p72_genBank_full_selected_p72_final_newNames.txt",sep="") )
	final_accns <- final_infos[,1]

# load in partial alignment
	seqs <- read.dna( paste(path,path2,"asfv_p72_newNames_partial_al.fas",sep=""), format="fasta", as.matrix=FALSE)
	taxa <- attributes(seqs)$names
	accns<- apply(as.matrix(taxa), 1, getEl, ind=1, sep="\\|")

	minds<- match(accns, final_accns)
	sinfo<- final_infos[minds,]
	sinfo<- cbind(sinfo,taxa)
	write.table(sinfo, file=paste(path,path2,"asfv_p72_newNames_partial_al_info.txt",sep=""), col.names=TRUE, row.names=FALSE)

	mtData <- matrixTranslate(seqs)
	ngapsAmbs <- countNuclsTypes( mtData )

# drop sequences with no date or country information
	ex_inds	<- which( (sinfo[,date_col] == -1) | (sinfo[,country_col] == -1) )
	inc_accns	<- setdiff(accns, accns[ex_inds])
	inc_inds	<- match(inc_accns, accns)

	write.dna( seqs[inc_inds], file=paste(path,path2,"asfv_p72_newNames_partial_incs.fas",sep=""), format="fasta", nbcol=-1, colsep="")
	write.table( sinfo[inc_inds,], file=paste(path,path2,"asfv_p72_newNames_partial_incs_info,txt",sep=""), col.names=TRUE, row.names=FALSE)

	seqs <- seqs[inc_inds]
	sinfo<- sinfo[inc_inds,]

# use sub-sampling
	taxa    <- attributes(seqs)$names
	country <- apply(as.matrix(taxa), 1, getEl, ind=2, sep="\\|")
	year	  <- apply(as.matrix(taxa), 1, getEl, ind=3, sep="\\|")
	traits1 <- country
	traits2 <- country
	traits3 <- year
	traits4 <- year
	chosen_inds <- fourWaySampler(traits1,traits2,traits3,traits4,maxPer=1)

	write.dna( seqs[chosen_inds], file=paste(path,path2,"asfv_p72_partial_onePerCountryPerYear.fas",sep=""), format="fasta", nbcol=1, colsep="")
	write.table( sinfo[chosen_inds,], file=paste(path,path2,"asfv_p72_partial_onePerCountryPerYear_info.txt",sep=""), col.names=TRUE, row.names=FALSE)
	
