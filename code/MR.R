suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(viper)))

mrDiscovery = function(ExS, regulon, gs, nsample = 50){

	ref = apply(ExS, 1, mean)
	dM  = ExS - ref
	score  = rep(NA, ncol(dM))
	
	# sample stratification
	sizeCheck = any(sapply(gs, length) < 10)
	if(sizeCheck){
		if(length(gs) == 1)
			score = apply(dM[rownames(dM) %in% gs[[1]], ], 2, mean)
		else if(length(gs) == 2)
			score = apply(dM[rownames(dM) %in% gs[[1]], ], 2, mean) - apply(dM[rownames(dM) %in% gs[[2]], ], 2, mean)
	} else {
		for(i in 1:ncol(dM)){
			temp = SimpleRankTest(dM[,i], gs)
			if(length(gs) == 1)
				score[i] = temp[1]
			else if(length(gs) == 2)
				score[i] = temp[1] - temp[2]
	}
	names(score) = colnames(dM)
	}

	# Master regulator discovery
	nsample = min(as.integer(ncol(dM)/2), 50)
	ExH = ExS[, names(sort(score, decreasing = TRUE)[1:nsample])]
	ExL = ExS[, names(sort(score, decreasing = FALSE)[1:nsample])]
	signature = rowTtest(ExH, ExL)
	signature = (qnorm(signature$p.value / 2, lower.tail = FALSE) * sign(signature$statistic))[, 1]
	nullmodel = ttestNull(ExH, ExL, per = 500, repos = TRUE, verbose = FALSE)
	mrs 	  = msviper(signature, regulon, nullmodel, verbose = FALSE)
	res = summary(mrs, 1000)
	return(res)
}

mrMulti = function(EXList, GRList, gs, nsample = 50){
	if(length(EXList) != length(GRList))
		stop("Size not matched.\n")
	for(i in 1:length(EXList)){
		cat("Running data set", names(EXList)[i], "\n")
		res = mrDiscovery(EXList[[i]], GRList[[i]], gs, nsample = nsample)
		res$dataset = names(EXList)[i]
		if(i == 1)
			mergedT = res
		else
			mergedT = rbind(mergedT, res)
	}
	return(mergedT)
}

SimpleRankTest <- function(Pvalue, gsList, size = 4) {
# The Pvalue is pvalue vector with Entrez ID as names
# gsList is the gene set list

	GeneID <- names(Pvalue)
	M <- length(gsList)
	zscore <- rep(NaN,M)
	
	N <- length(Pvalue)
	RankP <- rank(Pvalue,na.last = TRUE)
	
	for(i in 1:M) {

		Index <- GeneID %in% gsList[[i]] 
		n1 <- sum( Index )
		
		if(n1 < size)
			next
		n2 <- N - n1
		R1 <- sum(RankP[Index])
		K1 <- R1 - n1*(n1+1)/2
		Umean <- n1*n2/2
		Uvar  <- sqrt(n1*n2*(n1+n2+1)/12)
		zscore[i] <- (K1-Umean)/Uvar
	}

	# return values
	names(zscore) <- names(gsList)
	return(zscore)
}

mrSummary = function(res, pvalue = 0.01){
	uTable = filter(res, NES > 0, p.value < pvalue)
	uTable$group = "Activated"
	dTable = filter(res, NES < 0, p.value < pvalue)
	dTable$group = "Repressed"
	res = rbind(uTable, dTable)
	sumx = table(res$Regulon, res$group)
	return(sumx)
}

updateTable = function(cTable, space = 0:7, N00 = TRUE){
	cM = as.matrix(cTable)
	newM = matrix(0, length(space), length(space))
	dimnames(newM) = list(as.character(space), as.character(space))
	newM[rownames(cM), colnames(cM)] = cM
	if(N00)
		newM["0", "0"] = 0
	return(newM)
}
