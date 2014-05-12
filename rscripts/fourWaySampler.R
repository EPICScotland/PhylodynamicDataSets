# four way subsampler
# for H7N9 work; 19 Sept 2013
# S J Lycett
# 19 Sept 2013


	# helper function
	getIndexList <- function( ulistNames, listNames ) {
		inds	 <- vector("list", length(ulistNames))
		for (i in 1:length(ulistNames)) {
			inds[[i]] <- which(listNames == ulistNames[i])
		}
		return( inds )
	}

# MAIN function
fourWaySampler	<- function( traits1, traits2, traits3, traits4, maxPer=5 ) {

	utraits1	<- sort(unique(traits1))
	utraits2	<- sort(unique(traits2))
	utraits3	<- sort(unique(traits3))
	utraits4	<- sort(unique(traits4))

	t1_inds	<- getIndexList(utraits1, traits1)
	t2_inds	<- getIndexList(utraits2, traits2)
	t3_inds	<- getIndexList(utraits3, traits3)
	t4_inds	<- getIndexList(utraits4, traits4)

	chosen_inds	<- c()
	for (i in 1:length(utraits1)) {
	  inds1		<- t1_inds[[i]]

	  for (j in 1:length(utraits2)) {
		inds2	<- intersect(inds1, t2_inds[[j]] )

		if (length(inds2) > 0) {
		  for (k in 1:length(utraits3)) {
			inds3	<- intersect(inds2, t3_inds[[k]] )

			if (length(inds3) > 0) {
			  for (m in 1:length(utraits4)) {
				inds4 <- intersect(inds3, t4_inds[[m]])

				if (length(inds4) > 0) {
					if (length(inds4) <= maxPer) {
						chosen_inds <- c(chosen_inds,inds4)
					} else {
						chosen_inds <- c(chosen_inds,
									sample(inds4, maxPer) )
					}

				}
			  }
			}

		  }				
		}
	  }
	}

	return( chosen_inds )
}

	