#As in LumpTest.R except a change in optimizeRateGroup to allow single rate matrix
##########################################################################
LumpabilityTest <- function(fileDir, seqFile, parMatrixFile, NewickTreeFile, initialPi, grpSmat, grpPiVector, beta
			   , InvProbSame, pxInv, verbose, isReversible, Grouping, numItr=20){

	ResultNames <- paste("LumpResults", Grouping, ".txt", sep="")
	ResultNames <- c("OrigResultsN.txt", ResultNames)
	Grouping    <- c("N", Grouping)

	Comparisons <- cbind(ResultNames, Grouping)
	
	uniqueRs    <- paste(grpPiVector[-1], grpSmat)
	uniqueRs    <- unique(uniqueRs)
	dfreedom    <- 2*length(uniqueRs)

	Likelihoods <- apply(Comparisons, 1, AnalyzeGroup
		, fileDir = fileDir
		, seqFile = seqFile
		, parMatrixFile = parMatrixFile
		, NewickTreeFile = NewickTreeFile
		, initialPi = initialPi
		, grpSmat = grpSmat
		, grpPiVector = grpPiVector
		, beta = beta
		, InvProbSame = InvProbSame
		, pxInv = pxInv
		, verbose = verbose
		, isReversible = isReversible
		, numItr= numItr)


	UncnLikelihood  <- Likelihoods[1]
	ConsLikelihoods <- Likelihoods[-1]

	LRs         <- 2*(UncnLikelihood-ConsLikelihoods)
	Pvalues     <- 1-pchisq(LRs, dfreedom)
	StatResults <- cbind(LRs, Pvalues)
	row.names(StatResults) <- Grouping[-1]
	colnames(StatResults) <- c("LR", "P-value")

	inLumpResult <- file.path(fileDir, "Test_for_Lumpability.txt")
		write(paste("Lumping", "LR", "P-value", sep="\t"), inLumpResult)
		write.table(StatResults, inLumpResult, append=TRUE, row.names=TRUE, col.names=FALSE, sep="\t")



	StatResults

}
		
AnalyzeGroup <- function(GroupInfo, fileDir, seqFile, parMatrixFile, NewickTreeFile, initialPi, grpSmat, grpPiVector, beta
			   , InvProbSame, pxInv, verbose, isReversible, numItr){

		print(paste("Calculating loglikelihood for", GroupInfo[1]))
		print(GroupInfo[2])
		IndLikelihood <- estimateMleTest(fileDir, seqFile, parMatrixFile, NewickTreeFile, GroupInfo[1], initialPi, grpSmat
			  	, grpPiVector, beta, InvProbSame, pxInv, verbose, isReversible, GroupInfo[2], numItr)

		IndLikelihood
		
}


#########################################################################
#									#
# Codes to compute the statistics under unconstrained and constrained	#
#  conditions, likelihood ratio tests and PP-plots 			#
#									#
#########################################################################

#########################################################################################
#											#
# Get lumpability PP-plots and parameter estimate statistics				#
#											#
#########################################################################################

computeLumpStats <- function(mainDir, outFolder, lumpFolder, LR_File, isReversible){

	fileDir	    <- file.path(mainDir, outFolder)
	lumpDir	    <- file.path(mainDir, lumpFolder)
	Directories <- c(fileDir, lumpDir)

	OutNames    <- sapply(Directories, dir, pattern=".out")
	OutNames    <- matrix(OutNames, ncol=2)

	colnames(OutNames) <- NULL
	nReps	           <- dim(OutNames)[1]

	numchars    <- as.numeric(names(table(nchar(OutNames[,1]))))
	sortNames   <- NULL
	for(i in 1:length(numchars)){
		sortNames <- rbind(sortNames, OutNames[nchar(OutNames[, 1])==numchars[i], ])
	}

	BasicInfo  <- file.path(fileDir, sortNames[1,1])
	BasicInfo  <- scan(BasicInfo, what="character", sep="\r")
	kTaxa     <<- as.numeric(BasicInfo[2]) + 1
	cStat	  <<- length(unlist(strsplit(BasicInfo[3*kTaxa+9], split="\t")))
	numPar    <<- length(unlist(strsplit(BasicInfo[kTaxa+10], split="\t")))
	edges	  <<- 2*kTaxa-2
	if(isReversible) {freePar <<- numPar-(cStat+1)
		    }else{freePar <<- numPar-1}

	trueParameters <- GetTrueParamInfo(sortNames[1, 1], fileDir, isReversible)
	trueParameters <- round(trueParameters, 5)
	trueParMatrix  <- trueParameters[1:(edges*numPar)]
	trueParMatrix  <- matrix(trueParMatrix, edges, numPar, byrow=TRUE)

	if(isReversible){RateMats <- t(apply(trueParMatrix, 1, NormalizeRevPars))
		 } else {RateMats <- t(apply(trueParMatrix, 1, NormalizeNoRevPars))}

	Rtype      <- apply(RateMats[, -numPar], 1, paste, collapse="")
	degFreedom <- 2*length(unique(Rtype))

	UncnInfos <- sapply(sortNames[, 1], GetOneOutInfo, fileDir=fileDir, isReversible)
	LumpInfos <- sapply(sortNames[, 2], GetOneOutInfo, fileDir=lumpDir, isReversible)

	# Likelihood-ratio-related information and PP-plot
	##################################################
	LR_index  <- edges*numPar+2*cStat+4
	LR_stat   <- 2*(UncnInfos[LR_index, ]-LumpInfos[LR_index, ])

	P_values  <- 1-pchisq(LR_stat, degFreedom)

	PP_File   <- file.path(mainDir, "PP_plot.pdf")
		pdf(height=5, width=5, file=PP_File)
		Pvalues <- sort(P_values)
		UnifDis <- 1:nReps/(nReps+1)
		plot(UnifDis, Pvalues, ylab="Observed p values", xlab="Expected p values", main=" ", xlim=c(0,1), ylim=c(0,1), pch=19, col="red", cex=0.5)
		abline(0, 1, h=0.05, lwd=1)
		abline(0.05, 0, h=0.05, lwd=1)
		dev.off()

	LR_matrix <- cbind(LR_stat, P_values) 
	inLR_File <- file.path(mainDir, LR_File)
		write("Likelihood Ratio Test", inLR_File)
		write(paste("Number of degrees of freedom =", degFreedom), inLR_File, append = TRUE)
		write(" ", inLR_File, append = TRUE)

		write(paste("LR statistics and P-values"), inLR_File, append = TRUE)
		write(t(round(LR_matrix, 5)), inLR_File, ncolumns=2, append=TRUE, sep="\t")
	

	# Delta-related information
	###########################
	UncnDeltas  = UncnInfos[edges*numPar+2*cStat+2, ]
	LumpDeltas  = LumpInfos[edges*numPar+2*cStat+2, ]

		DeltaInfo  <- cbind(UncnDeltas, LumpDeltas)
		row.names(DeltaInfo) <- sortNames[, 1]

		inDeltaFile <- file.path(mainDir, "deltaLog Info.txt")
			write("Delta LogLs for unconstrained and lumpability conditions", inDeltaFile)
			write(paste("", "Unconstrained", "Lumpability", sep="\t"), inDeltaFile, append=TRUE)
			write.table(DeltaInfo, inDeltaFile, append=TRUE, row.names=TRUE, col.names=FALSE, sep="\t")

	# Iteration related information
	###############################
	UncnIters   = UncnInfos[edges*numPar+2*cStat+3, ]			
	ConnIters   = LumpInfos[edges*numPar+2*cStat+3, ]

		Iterations   = cbind(UncnIters, ConnIters)

		allItersFreq = table(c(UncnIters, ConnIters))
		uncItersFreq = table(UncnIters)
		conItersFreq = table(ConnIters)
		
		freqSummary  = matrix(0, length(allItersFreq), 2)
		row.names(freqSummary) = names(allItersFreq)

		freqSummary[row.names(uncItersFreq), 1] = table(UncnIters)
		freqSummary[row.names(conItersFreq), 2] = table(ConnIters)
		freqSummary = cbind(as.numeric(row.names(freqSummary)), freqSummary)
		
	
		inIterFile <- file.path(mainDir, "Iterations.txt")
			write("Number of iterations for unconstrained and lumpability files", inIterFile)
			write.table(Iterations, inIterFile, append=TRUE, row.names=TRUE, col.names=FALSE, sep="\t")
			write(" ", inIterFile, append=TRUE)
			write("Number of iterations frequency table", inIterFile, append=TRUE)
			write(paste("Iterations", "Unconstrained", "Lumpability", sep="\t"), inIterFile, append=TRUE)
			write.table(freqSummary, inIterFile, append=TRUE, row.names=FALSE, col.names=FALSE, sep="\t")

	# Replicate analysis related information
	########################################
	UncnMeans       <- apply(UncnInfos, 1, mean)
	UncnSds         <- apply(UncnInfos, 1, sd)
	UncnStats       <- paste(trueParameters, round(UncnMeans, 5), round(UncnSds/sqrt(nReps), 5), round(UncnSds/UncnMeans, 5), sep="|")

		parMatUncnStats <- UncnStats[1:(edges*numPar)]
		parMatUncnStats <- matrix(parMatUncnStats, edges, numPar, byrow=TRUE)
		rootUncnStats   <- UncnStats[(edges*numPar+1):(edges*numPar+cStat)]
		betaUncnStats   <- UncnStats[(edges*numPar+cStat+1)]
		invarUncnStats  <- UncnStats[(edges*numPar+cStat+2):(edges*numPar+2*cStat+1)]
	
		inOutFile <- file.path(mainDir, "STATS_uncnPar.txt")
	
			write("Unconstrained parameter estimation statistics", inOutFile)		
			write(" ", inOutFile, append=TRUE)
			write("Information per parameter: True parameter, Sample mean, Standard error, Coefficient of variation", inOutFile, append=TRUE)
			write(" ", inOutFile, append=TRUE)
	
			if(isReversible){
				write("S-matrix statistics", inOutFile, append=TRUE)
				write(t(parMatUncnStats[, 1:freePar]), inOutFile, ncolumns=freePar, append=TRUE, sep="\t")
				write(" ", inOutFile, append=TRUE)
				write("PI-matrix statistics", inOutFile, append=TRUE)
				write(t(parMatUncnStats[, (freePar+1):(freePar+cStat)]), inOutFile, ncolumns=cStat, append=TRUE, sep="\t")
			}else{
				write("R-matrix statistics", inOutFile, append=)
				write(t(parMatUncnStats[, 1:freePar]), inOutFile, ncolumns=freePar, append=TRUE, sep="\t")
			}
	
			write(" ", inOutFile, append=TRUE)
			write("Edge length statistics", inOutFile, append = TRUE)
			write(t(parMatUncnStats[, numPar]), inOutFile, ncolumns=1, append=TRUE, sep="\t")

			write(" ", inOutFile, append=TRUE)
			write("Root-vector statistcs", inOutFile, append = TRUE)
			write(rootUncnStats, inOutFile, ncolumns=cStat, append=TRUE, sep="\t")
	
			write(" ", inOutFile, append=TRUE)
			write("Beta statistics", inOutFile, append=TRUE)
			write(betaUncnStats, inOutFile, append=TRUE)
	
			write(" ", inOutFile, append=TRUE)
			write("Prob of the elements of the state space among the invariant sites", inOutFile, append=TRUE)
			write(invarUncnStats, inOutFile, ncolumns=cStat, append=TRUE, sep="\t")
	
			write(" ", inOutFile, append=TRUE)
			write("Total number of replicates", inOutFile, append=TRUE)
			write(nReps, inOutFile, append=TRUE)


	LumpMeans <- apply(LumpInfos, 1, mean)
	LumpSds   <- apply(LumpInfos, 1, sd)
	LumpStats <- paste(trueParameters, round(LumpMeans, 5), round(LumpSds/sqrt(nReps), 5), round(LumpSds/LumpMeans, 5), sep="|")

		parMatLumpStats <- LumpStats[1:(edges*numPar)]
		parMatLumpStats <- matrix(parMatLumpStats, edges, numPar, byrow=TRUE)
		rootLumpStats   <- LumpStats[(edges*numPar+1):(edges*numPar+cStat)]
		betaLumpStats   <- LumpStats[(edges*numPar+cStat+1)]
		invarLumpStats  <- LumpStats[(edges*numPar+cStat+2):(edges*numPar+2*cStat+1)]

		inLumpFile <- file.path(mainDir, "STATS lumpPar Stats.txt")
	
			write("Lumpability-constrained parameter estimation statistics", inLumpFile)		
			write(" ", inLumpFile, append=TRUE)
			write("Information per parameter: True parameter, Sample mean, Standard error, Coefficient of variation", inLumpFile, append=TRUE)
			write(" ", inLumpFile, append=TRUE)
	
			if(isReversible){
				write("S-matrix statistics", inLumpFile, append=TRUE)
				write(t(parMatLumpStats[, 1:freePar]), inLumpFile, ncolumns=freePar, append=TRUE, sep="\t")
				write(" ", inLumpFile, append=TRUE)
				write("PI-matrix statistics", inLumpFile, append=TRUE)
				write(t(parMatLumpStats[, (freePar+1):(freePar+cStat)]), inLumpFile, ncolumns=cStat, append=TRUE, sep="\t")
			}else{
				write("R-matrix statistics", inLumpFile, append=)
				write(t(parMatLumpStats[, 1:freePar]), inLumpFile, ncolumns=freePar, append=TRUE, sep="\t")
			}
	
			write(" ", inLumpFile, append=TRUE)
			write("Edge length statistcs", inLumpFile, append = TRUE)
			write(t(parMatLumpStats[, numPar]), inLumpFile, ncolumns=1, append=TRUE, sep="\t")

			write(" ", inLumpFile, append=TRUE)	
			write("Root-vector statistics", inLumpFile, append = TRUE)
			write(rootLumpStats, inLumpFile, ncolumns=cStat, append=TRUE, sep="\t")
	
			write(" ", inLumpFile, append=TRUE)
			write("Beta statistics", inLumpFile, append=TRUE)
			write(betaLumpStats, inLumpFile, append=TRUE)
	
			write(" ", inLumpFile, append=TRUE)
			write("Prob of the elements of the state space among the invariant sites", inLumpFile, append=TRUE)
			write(invarLumpStats, inLumpFile, ncolumns=cStat, append=TRUE, sep="\t")
	
			write(" ", inLumpFile, append=TRUE)
			write("Total number of replicates", inLumpFile, append=TRUE)
			write(nReps, inLumpFile, append=TRUE)

	LR_matrix
}



#####################################################################
# Codes to compute the statistics from a single folder		    #
#####################################################################

computeStatistics <- function(fileDir, isReversible, meanStdFile, deltaLogFile, itrFile){

	OutNames  <- dir(fileDir, pattern=".out")
	nReps	  <- length(OutNames)
	numchars  <- as.numeric(names(table(nchar(OutNames))))
	sortNames <- NULL
	for(i in 1:length(numchars)){
		sortNames <- c(sortNames, OutNames[nchar(OutNames)==numchars[i]])
	}
	sortNames

	BasicInfo  <- file.path(fileDir, sortNames[1])
	BasicInfo  <- scan(BasicInfo, what="character", sep="\r")
	kTaxa     <<- as.numeric(BasicInfo[2]) + 1
	cStat	  <<- length(unlist(strsplit(BasicInfo[3*kTaxa+9], split="\t")))
	numPar    <<- length(unlist(strsplit(BasicInfo[kTaxa+10], split="\t")))
	edges	  <<- 2*kTaxa-2
	if(isReversible) {freePar <<- numPar-(cStat+1)
		    }else{freePar <<- numPar-1}

	truePars  <- GetTrueParamInfo(OutNames[1], fileDir, isReversible)
	allInfos  <- sapply(sortNames, GetOneOutInfo, fileDir=fileDir, isReversible)

	# Delta-related information
	###########################
	Deltas        <- allInfos[edges*numPar+2*cStat+2, ]
	names(Deltas) <- sortNames

	inDeltaFile <- file.path(fileDir, deltaLogFile)
		write("Delta LogL", inDeltaFile)
		write.table(Deltas, inDeltaFile, append=TRUE, row.names=TRUE, col.names=FALSE, sep="\t")

	# Iteration related information
	###############################
	Iterations = allInfos[dim(allInfos)[1], ]
	itrSummary = table(Iterations)

	numIters   = names(itrSummary)
	itrSummary = cbind(as.numeric(numIters), itrSummary)

	inIterFile <- file.path(fileDir, itrFile)
		write("Number of iterations per file", inIterFile)
		write.table(Iterations, inIterFile, append=TRUE, row.names=TRUE, col.names=FALSE)
		write(" ", inIterFile, append=TRUE)
		write("Number of iterations frequency table", inIterFile, append=TRUE)
		write(paste("Iterations", "Frequency", sep="\t"), inIterFile, append=TRUE)
		write.table(itrSummary, inIterFile, append = TRUE, row.names=FALSE, col.names=FALSE, sep="\t")

	# Replicate analysis related information
	########################################
	Means    <- apply(allInfos, 1, mean)
	Sds      <- apply(allInfos, 1, sd)
	Stats    <- paste(round(truePars, 5), round(Means, 5), round(Sds/sqrt(nReps), 5), round(Sds/Means, 5)
	                , (abs(Means-truePars)<Sds/sqrt(nReps)), (abs(Means-truePars) < 2*Sds/sqrt(nReps))
			, (abs(Means-truePars)<3*Sds/sqrt(nReps)), sep="|")

	parMatStats <- Stats[1:(edges*numPar)]
	parMatStats <- matrix(parMatStats, edges, numPar, byrow=TRUE)
	rootStats   <- Stats[(edges*numPar+1):(edges*numPar+cStat)]
	betaStats   <- Stats[(edges*numPar+cStat+1)]
	invarStats  <- Stats[(edges*numPar+cStat+2):(edges*numPar+2*cStat+1)]
	
	inOutFile <- file.path(fileDir, meanStdFile)
	
		write("Parameter estimation statistics", inOutFile)		
		write(" ", inOutFile, append=TRUE)
		write("Information per parameter: True parameter, Sample mean, Standard error, Coefficient of variation, does the estimate fall within one SE, 2SE, 3SE", inOutFile, append=TRUE)
		write(" ", inOutFile, append=TRUE)
	
		if(isReversible){
			write("S-matrix statistics", inOutFile, append=TRUE)
			write(t(parMatStats[, 1:freePar]), inOutFile, ncolumns=freePar, append=TRUE, sep="\t")
			write(" ", inOutFile, append=TRUE)
			write("PI-matrix statistics", inOutFile, append=TRUE)
			write(t(parMatStats[, (freePar+1):(freePar+cStat)]), inOutFile, ncolumns=cStat, append=TRUE, sep="\t")
		}else{
			write("R-matrix statistics", inOutFile, append=)
			write(t(parMatStats[, 1:freePar]), inOutFile, ncolumns=freePar, append=TRUE, sep="\t")
			}
	
		write(" ", inOutFile, append=TRUE)
		write("Edge length statistics", inOutFile, append = TRUE)
		write(t(parMatStats[, numPar]), inOutFile, ncolumns=1, append=TRUE, sep="\t")

		write(" ", inOutFile, append=TRUE)
		write("Root-vector statistcs", inOutFile, append = TRUE)
		write(rootStats, inOutFile, ncolumns=cStat, append=TRUE, sep="\t")
	
		write(" ", inOutFile, append=TRUE)
		write("Beta statistics", inOutFile, append=TRUE)
		write(betaStats, inOutFile, append=TRUE)
	
		write(" ", inOutFile, append=TRUE)
		write("Prob of the elements of the state space among the invariant sites", inOutFile, append=TRUE)
		write(invarStats, inOutFile, ncolumns=cStat, append=TRUE, sep="\t")
	
		write(" ", inOutFile, append=TRUE)
		write("Total number of replicates", inOutFile, append=TRUE)
		write(nReps, inOutFile, append=TRUE)
}

#########################################################################
# Extracting the true parameter information from one replicate		#
# Input global parameters						#
#	kTaxa								#
#	edges
#	numPar
#	freePar
#########################################################################

GetTrueParamInfo <- function(OutFile, fileDir, isReversible){

	BasicInfo = file.path(fileDir, OutFile)
	BasicInfo = scan(BasicInfo, what="character", sep="\r")
	ParMat    = Re(as.numeric(unlist(strsplit(BasicInfo[(kTaxa+10):(3*kTaxa+7)], split="\t"))))
	ParMat    = matrix(ParMat, edges, numPar, byrow=TRUE)

	if(isReversible){ParMat = apply(ParMat, 1, NormalizeRevPars)
		   }else{ParMat = apply(ParMat, 1, NormalizeNoRevPars) }

	Root      = as.numeric(unlist(strsplit(BasicInfo[3*kTaxa+9], split="\t")))
	Beta      = as.numeric(BasicInfo[3*kTaxa+11])
	Invariant = as.numeric(unlist(strsplit(BasicInfo[3*kTaxa+13], split="\t")))
	AllInfo   = c(ParMat, Root, Beta, Invariant, 0, 0)
	AllInfo
}


#################################################################
#								#
# Extracting the basic information from a single file		#
#								#
# Input global parameters					#
#	kTaxa							#
#################################################################

GetOneOutInfo <- function(OutFile, fileDir, isReversible){

	BasicInfo = file.path(fileDir, OutFile)
	BasicInfo = scan(BasicInfo, what="character", sep="\r")

	Deltas     = as.numeric(unlist(strsplit(BasicInfo[3*kTaxa+19], split="\t")))
	DeltaDiff  = Deltas[2]-Deltas[1]
	Likelihood = Deltas[1]

	ItrNumber = as.numeric(BasicInfo[3*kTaxa+15])

	ParMat    = unlist(strsplit(BasicInfo[(3*kTaxa+21):(5*kTaxa+18)], split="\t"))
	if(!isReversible) ParMat = gsub("\\+0i", "", ParMat); #print(ParMat)
	ParMat	  = as.numeric(ParMat)

	Root      = as.numeric(unlist(strsplit(BasicInfo[5*kTaxa+20], split="\t")))
	Beta      = as.numeric(BasicInfo[5*kTaxa+22])
	Invariant = as.numeric(unlist(strsplit(BasicInfo[5*kTaxa+24], split="\t")))
	AllInfo   = c(ParMat, Root, Beta, Invariant, DeltaDiff, ItrNumber, Likelihood)
	AllInfo
}


#########################################################################################
# 						    					#
# Codes to check lumpability			    					#
# V matrix must be an input parameter		    					#
# If all matrices are lumpable, the output should be a matrix with only			#
# zeros											#
# 						    					#
#########################################################################################

####################################################################################
# Calculates the frobenius norm of the matrix VUPV-PV for all the output.out files #
# in a folder.									   #
#										   #
# If the data was generated under lumpable conditions, the output should be a	   #
# matrix with only zeros							   #
####################################################################################

CheckLumpability <- function(lumpDir, isReversible, Vmat){

	lumpNames  <- dir(lumpDir, pattern=".out")
	nReps	  <- length(lumpNames)
	numchars  <- as.numeric(names(table(nchar(lumpNames))))
	SortNames <- NULL
	for(i in 1:length(numchars)){
		SortNames <- c(SortNames, lumpNames[nchar(lumpNames)==numchars[i]])
	}
	SortNames

	BasicInfo  <- file.path(lumpDir, SortNames[1])
	BasicInfo  <- scan(BasicInfo, what="character", sep="\r")
	kTaxa     <<- as.numeric(BasicInfo[2]) + 1
	cStat	  <<- length(unlist(strsplit(BasicInfo[3*kTaxa+9], split="\t")))
	numPar    <<- length(unlist(strsplit(BasicInfo[kTaxa+10], split="\t")))
	edges	  <<- 2*kTaxa-2
	if(isReversible) {freePar <<- numPar-(cStat+1)
		    }else{freePar <<- numPar-1}
	
	V      <<- Vmat
	sapply(SortNames, GetLumpMat, lumpDir=lumpDir, isReversible)

}

#########################################################################################
# It calculates the frobenius norm for each parameter set of a parameter matrix		#
#########################################################################################
GetLumpMat <- function(lumpFile, lumpDir, isReversible){
	BasicInfo = file.path(lumpDir, lumpFile)
	BasicInfo = scan(BasicInfo, what="character", sep="\r")
	#ParMat    = Re(as.numeric(unlist(strsplit(BasicInfo[(3*kTaxa+21):(5*kTaxa+18)], split="\t"))))

	ParMat    = unlist(strsplit(BasicInfo[(3*kTaxa+21):(5*kTaxa+18)], split="\t"))
	if(!isReversible) ParMat = gsub("\\+0i", "", ParMat)
	ParMat	  = as.numeric(ParMat)

	ParMat    = matrix(ParMat, ncol=numPar, byrow=TRUE)
	Mindex    = apply(ParMat, 1, getM, isReversible=isReversible, V=V)
	Mindex[Mindex <= 1e-6 ] = 0
	Mindex
}

#########################################################################################
# It calculates the frobenius norm for a single set of parameters			#
#########################################################################################
getM = function(pars, isReversible, V){
	if(isReversible){
		PI = diag(pars[(freePar+1):(freePar+cStat)])
		P  = RevPmat(pars)
	} else {
		PI = diag(GetStatDist(pars))
		P  = NoRevPmat(pars)
	}
	U  = solve(t(V)%*%PI%*%V)%*%t(V)%*%PI
	M  = V%*%U%*%P%*%V - P%*%V
	sum(sqrt(M^2))
}


####################################################
#						   #
# Codes to analyze bootstrapped replicates	   #
#						   #
####################################################

analyzeTest <- function(fileDir, templateFile, BTFile, isReversible, InvProbSame, verbose, Group=NULL, numItr){

	TempInfo	<- scanTemplate(fileDir, templateFile)
	#
	# From previous command I got this info
	#	kTaxa
	#	TaxaNames
	#	mergeMatrix
	#	mergeVector
	#	rMerg
	#	parMatrix
	#	Beta
	#	ProbInv
	#	Pi0
	#	cStat
	#	edges

	if(is.null(Group)) Group="N" else Group=Group

	Group		<<- Group
	InitParMatrix	<<- parMatrix
	Beta0		<<- Beta
	alpha		<<- 1-Beta
	rootVector	<<- Pi0		

	TempFile	 <- file.path(fileDir, templateFile)
	MoreInfo	 <- scan(TempFile, what="character", sep="\r")
	grpSmat		<<- as.numeric(unlist(strsplit(MoreInfo[kTaxa+6], split="\t")))
	 sMatRows	<<- split(1:edges, grpSmat)
	 sMatGroups	<<- 1:length(sMatRows)

	grpPiVector	<<- as.numeric(unlist(strsplit(MoreInfo[kTaxa+8], split="\t")))
	 piVecRows	<<- split(1:edges, grpPiVector[-1])
	 piVecGroups	<<- 1:length(piVecRows)
	 rootVecRow	<<- grpPiVector[1]


	# These indexes will be used to extract info from the unique pattern Info matrix (Not yet generated).
	pConst		<<- kTaxa+1
	Countings	<<- kTaxa+2

	BTAccess	 <- file.path(fileDir, BTFile)
	BTInfo		<<- scan(BTAccess, what="character", sep="\r")
	nSites		<<- as.numeric(unlist(strsplit(BTInfo[1], split="\t"))[2])
	nReps		<<- length(BTInfo)/(kTaxa+1)

	firstBT		 <- BTInfo[2:(kTaxa+1)]
	firstBT		 <- substr(firstBT, 11, nSites+10)
	stateSpace	<<- names(table(unlist(strsplit(firstBT, split=""))))
	numericS	<<- 1:cStat

	#seqIndex <- seq(1, 90, kTaxa+1)	
	#BTRepl	 <- cbind(1:10, seqIndex)

	#seqIndex <- seq(91, 180, kTaxa+1)	
	#BTRepl	 <- cbind(11:20, seqIndex)

	#seqIndex <- seq(181, 270, kTaxa+1)	
	#BTRepl	 <- cbind(21:30, seqIndex)

	#seqIndex <- seq(271, 360, kTaxa+1)	
	#BTRepl	 <- cbind(31:40, seqIndex)

	#seqIndex <- seq(361, length(BTInfo), kTaxa+1)	
	#BTRepl	 <- cbind(41:nReps, seqIndex)

	seqIndex <- seq(1, length(BTInfo), kTaxa+1)
	BTRepl	 <- cbind(1:nReps, seqIndex)

	if(isReversible) { freePar <<- cStat*(cStat-1)/2; numPar <<- freePar+cStat+1
		    }else{ freePar <<- cStat*(cStat-1)  ; numPar <<- freePar+1 }


	if(Group=="N") {
		dir.create(file.path(fileDir, "Output Replicates"), showWarnings = FALSE)
		BTDir	<- file.path(fileDir, "Output Replicates")
	}else{
		FolderName <-  paste("Lumped Replicates", Group)
		dir.create(file.path(fileDir, FolderName), showWarnings = FALSE)
		BTDir      <- file.path(fileDir, FolderName)
	}

	BTDir <<- BTDir

	apply(BTRepl, 1, getOutputFile
		, isReversible	= isReversible
		, InvProbSame	= InvProbSame
		, displayInfo	= verbose
		, numItr	= numItr)

}

##################################################################
# Data analysis for each individual sequence in the BT file
#
# Input global parameters
#	kTaxa
#	BTInfo
#	nSites
#	cStat
#	stateSpace
#	numericS
#	ProbInv
#
#	sMatRows   = split(1:edges, grpSmat)
#	sMatGroups = 1:length(sMatRows)
#
#	piVecRows   = split(1:edges, grpPiVector[-1])
#	piVecGroups = 1:length(piVecRows)
#	rootVecRow  = grpPiVector[1]
#
#	Beta0
#	InitParMatrix
#	BTDir
#
#	Beta0
#	TaxaNames
#	Group	Group for lumping the state space
#
###################################################################

getOutputFile <- function(replicateInfo, isReversible, InvProbSame, displayInfo, numItr){
	repNumber <- replicateInfo[1]
	repIndex  <- replicateInfo[2]

	seqInfo   <- BTInfo[(repIndex+1):(repIndex+kTaxa)]
	seqMatrix <- justGetSequences(seqInfo)

	taxaNames <- row.names(seqMatrix)

	seqMatrix     <- ChangeTaxaOrder(seqMatrix, TaxaNames)	

	patternInfo   <- getUniquePatternsL(seqMatrix, ProbInv)
	patternMat   <<- patternInfo$patMat
	InvOrder     <<- patternInfo$InvOrder
	rearrPxInv   <<- (1:cStat)[order(InvOrder)]

	if(displayInfo) print(paste("Replicate number", repNumber))
	Loglikelihood <- mleLump(parMatrix    = InitParMatrix
			        , isReversible = isReversible
				, InvProbSame  = InvProbSame
				, numItr       = numItr
				, displayInfo = displayInfo)

	OutFileName <- paste("Output", repNumber, ".out", sep="")

	write.info(fileDir	= BTDir
		, OutFile	= OutFileName
		, taxaNames	= TaxaNames
		, stateSpace    = stateSpace
		, grpSmat	= grpSmat
		, grpPiVector	= grpPiVector
		, parMatrix	= InitParMatrix
		, initialPi	= Pi0
		, pxInv		= ProbInv
		, numItr	= Loglikelihood$totalItr
		, outLogL	= Loglikelihood$currentLogL
		, outUncLogL    = Loglikelihood$unconsLogL
		, outparMatrix  = round(Loglikelihood$ParMat, 6)
		, outRootVector = Loglikelihood$rootVector
		, outBeta	= Loglikelihood$Beta
		, outpxInv	= Loglikelihood$pxInv)

}

#########################################################################
# Getting the sequence matrix with numeric alphabet from an alignment	#
#
# Input global parameters
#	nSites
#	kTaxa
#	cStat
#	stateSpace
#	numericS
#
# Definde global parameter
# 	
#	Global parameter taxaNames is defined here
#########################################################################
justGetSequences <- function(seqInfo){

	taxaNames  <- substr(seqInfo, 1, 10)
	taxaNames  <- gsub(" ", "", taxaNames)

	seqInfo    <- substr(seqInfo, 11, nSites+10)

	seqMat     <- NULL
	for(i in 1:kTaxa){
		seqMat.i = unlist(strsplit(seqInfo[i], NULL))
		seqMat   = rbind(seqMat, seqMat.i)
	}

	seqMat		  <- apply(seqMat, 2, DirectRecode, stateSpace, numericS, cStat)
	row.names(seqMat) <- taxaNames

	seqMat

}

#################################################
#						#
# LUMPABILITY-RELATED CODES			#
#						#
#################################################
#################################################
# Gy.to.S					#
#						#
# Restriction matrix for reversible cases	#
#						#
# Input global parameters			#
#	freePar	s				#
#	cStat	c				#
#################################################

Gy.to.S=function(parameters, Group){

	pix = parameters[(freePar+1):(freePar+cStat)]
	pix = pix/sum(pix)
	ss  = parameters[1:freePar]

	if(Group=="B") {A=matrix(c(  1,  -1,   0,   0,   0,  0, 0.5, 0.5,  -1,   0,   0,  0,   1,   1,   1,   0,   0,  0,   0,   0,   0,   1,   0,  0,   0,   0,   0,   0,   1,  0,   0,   0,   0,   0,   0,  1), 6, 6, byrow=TRUE)}
	if(Group=="D") {A=matrix(c(  1,   0,   0,  -1,   0,  0, 0.5,   0,   0, 0.5,  -1,  0,   1,   0,   0,   1,   1,  0,   0,   1,   0,   0,   0,  0,   0,   0,   1,   0,   0,  0,   0,   0,   0,   0,   0,  1), 6, 6, byrow=TRUE)}
	if(Group=="H") {A=matrix(c(  0,   1,   0,  -1,   0,  0,   0, 0.5,   0, 0.5,   0, -1,   0,   1,   0,   1,   0,  1,   1,   0,   0,   0,   0,  0,   0,   0,   1,   0,   0,  0,   0,   0,   0,   0,   1,  0), 6, 6, byrow=TRUE)}
	if(Group=="V") {A=matrix(c(  0,   0,   1,   0,  -1,  0,   0,   0, 0.5,   0, 0.5, -1,   0,   0,   1,   0,   1,  1,   1,   0,   0,   0,   0,  0,   0,   1,   0,   0,   0,  0,   0,   0,   0,   1,   0,  0), 6, 6, byrow=TRUE)}
	if(Group=="R") {A=matrix(c(  1,   0,   0,  -1,   0,  0,   0,   0,   1,   0,   0, -1,   1,   0,   1,   1,   0,  1,   1,   0,  -1,   1,   0, -1,   0,   1,   0,   0,   0,  0,   0,   0,   0,   0,   1,  0), 6, 6, byrow=TRUE)}
	if(Group=="Y") {A=matrix(c(  1,   0,  -1,   0,   0,  0,   0,   0,   0,   1,   0, -1,   1,   0,   1,   1,   0,  1,   1,   0,   1,  -1,   0, -1,   0,   1,   0,   0,   0,  0,   0,   0,   0,   0,   1,  0), 6, 6, byrow=TRUE)}
	if(Group=="W") {A=matrix(c(  1,   0,   0,   0,  -1,  0,   0,   1,   0,   0,   0, -1,   1,   1,   0,   0,   1,  1,   1,  -1,   0,   0,   1, -1,   0,   0,   1,   0,   0,  0,   0,   0,   0,   1,   0,  0), 6, 6, byrow=TRUE)}
	if(Group=="S") {A=matrix(c(  1,  -1,   0,   0,   0,  0,   0,   0,   0,   0,   1, -1,   1,   1,   0,   0,   1,  1,   1,   1,   0,   0,  -1, -1,   0,   0,   1,   0,   0,  0,   0,   0,   0,   1,   0,  0), 6, 6, byrow=TRUE)}
	if(Group=="K") {A=matrix(c(  0,   1,  -1,   0,   0,  0,   0,   0,   0,   1,  -1,  0,   0,   1,   1,   1,   1,  0,   0,   1,   1,  -1,  -1,  0,   1,   0,   0,   0,   0,  0,   0,   0,   0,   0,   0,  1), 6, 6, byrow=TRUE)}
	if(Group=="M") {A=matrix(c(  0,   1,   0,  -1,   0,  0,   0,   0,   1,   0,  -1,  0,   0,   1,   1,   1,   1,  0,   0,   1,  -1,   1,  -1,  0,   1,   0,   0,   0,   0,  0,   0,   0,   0,   0,   0,  1), 6, 6, byrow=TRUE)}
	if(Group=="RY"){A=matrix(c(pix[2],     0,  pix[4], -pix[2],       0, -pix[4], pix[1],      0,  -pix[1], pix[3],       0, -pix[3], 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1,  0, -1, -1,  0, 1), 6, 6, byrow=TRUE)}
	if(Group=="SW"){A=matrix(c(pix[2], pix[3],     0,       0,  -pix[2], -pix[3], pix[1], -pix[1],       0,      0,  pix[4], -pix[4], 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, -1,  0,  0, -1, 1), 6, 6, byrow=TRUE)}
	if(Group=="KM"){A=matrix(c(     0, pix[3], pix[4], -pix[3], -pix[4],       0,      0,  pix[1], -pix[1], pix[2], -pix[2],       0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0,  1, -1, -1,  1, 0), 6, 6, byrow=TRUE)}

	A     <- gram.schmidt(A)
	y     <- c(0, 0, (A %*% ss)[-(1:2), 1])
	s.hat <- c(solve(A)%*%y)
	s.hat
}	

#################################################
# Gy.to.R					#
#						#
# Restriction matrix for non-reversible cases	#
#						#
# Input global parameters			#
#	freePar	r				#
#	cStat	c				#
#################################################	

Gy.to.R=function(parameters, Group){

	if(Group=="B") {A=matrix(c(1,  -1,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0.5, 0.5, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0
				,  1,   1,  1, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0,   0,  0, 1, 0, 0, 0, 0, 0, 0, 0, 0	
				,  0,   0,  0, 0, 1, 0, 0, 0, 0, 0, 0, 0,  0,   0,  0, 0, 0, 1, 0, 0, 0, 0, 0, 0
				,  0,   0,  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0,   0,  0, 0, 0, 0, 0, 1, 0, 0, 0, 0
				,  0,   0,  0, 0, 0, 0, 0, 0, 1, 0, 0, 0,  0,   0,  0, 0, 0, 0, 0, 0, 0, 1, 0, 0
				,  0,   0,  0, 0, 0, 0, 0, 0, 0, 0, 1, 0,  0,   0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 1), 12, 12, byrow=TRUE)}
	if(Group=="D") {A=matrix(c(0, 0, 0,   1,  -1,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.5, 0.5, -1, 0, 0, 0, 0, 0, 0
				 , 0, 0, 0,   1,   1,  1, 0, 0, 0, 0, 0, 0, 1, 0, 0,   0,   0,  0, 0, 0, 0, 0, 0, 0
				 , 0, 1, 0,   0,   0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 1,   0,   0,  0, 0, 0, 0, 0, 0, 0
				 , 0, 0, 0,   0,   0,  0, 1, 0, 0, 0, 0, 0, 0, 0, 0,   0,   0,  0, 0, 1, 0, 0, 0, 0
				 , 0, 0, 0,   0,   0,  0, 0, 0, 1, 0, 0, 0, 0, 0, 0,   0,   0,  0, 0, 0, 0, 1, 0, 0
				 , 0, 0, 0,   0,   0,  0, 0, 0, 0, 0, 1, 0, 0, 0, 0,   0,   0,  0, 0, 0, 0, 0, 0, 1), 12, 12, byrow=TRUE)}
	if(Group=="H") {A=matrix(c(0, 0, 0, 0, 0, 0,   1,  -1,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.5, 0.5, -1, 0, 0, 0
				 , 0, 0, 0, 0, 0, 0,   1,   1,  1, 0, 0, 0, 1, 0, 0, 0, 0, 0,   0,   0,  0, 0, 0, 0
				 , 0, 1, 0, 0, 0, 0,   0,   0,  0, 0, 0, 0, 0, 0, 1, 0, 0, 0,   0,   0,  0, 0, 0, 0
				 , 0, 0, 0, 1, 0, 0,   0,   0,  0, 0, 0, 0, 0, 0, 0, 0, 1, 0,   0,   0,  0, 0, 0, 0
				 , 0, 0, 0, 0, 0, 1,   0,   0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0,   0,  0, 1, 0, 0
				 , 0, 0, 0, 0, 0, 0,   0,   0,  0, 0, 1, 0, 0, 0, 0, 0, 0, 0,   0,   0,  0, 0, 0, 1), 12, 12, byrow=TRUE)}
	if(Group=="V") {A=matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 0,   1,  -1,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.5, 0.5, -1
				 , 0, 0, 0, 0, 0, 0, 0, 0, 0,   1,   1,  1, 1, 0, 0, 0, 0, 0, 0, 0, 0,   0,   0,  0
				 , 0, 1, 0, 0, 0, 0, 0, 0, 0,   0,   0,  0, 0, 0, 1, 0, 0, 0, 0, 0, 0,   0,   0,  0
				 , 0, 0, 0, 1, 0, 0, 0, 0, 0,   0,   0,  0, 0, 0, 0, 0, 1, 0, 0, 0, 0,   0,   0,  0
				 , 0, 0, 0, 0, 0, 1, 0, 0, 0,   0,   0,  0, 0, 0, 0, 0, 0, 0, 1, 0, 0,   0,   0,  0
				 , 0, 0, 0, 0, 0, 0, 0, 1, 0,   0,   0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 1,   0,   0,  0), 12, 12, byrow=TRUE)}
	if(Group=="R") {A=matrix(c(0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 1, 0, -1
				 , 0, 0, 0, 1,  1, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 1, 0,  1
				 , 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0, 1, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0
				 , 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0,  0, 1, 0, 0, 0, 0, 0,  0
				 , 0, 0, 0, 0,  0, 0, 1, 0, 0, 0, 0,  0, 0, 0, 0, 0,  0, 0, 0, 1, 0, 0, 0,  0
				 , 0, 0, 0, 0,  0, 0, 0, 0, 1, 0, 0,  0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 1,  0), 12, 12, byrow=TRUE)}
	if(Group=="Y") {A=matrix(c(1, 0, -1, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 1, -1, 0, 0, 0
				 , 1, 0,  1, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 1,  1, 0, 0, 0
				 , 0, 1,  0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0,  0, 1, 0, 0, 0, 0,  0, 0, 0, 0
				 , 0, 0,  0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0,  0, 0, 0, 1, 0, 0,  0, 0, 0, 0
				 , 0, 0,  0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0,  0, 1, 0, 0
				 , 0, 0,  0, 0, 0, 0, 0, 0,  0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0,  0, 0, 0, 1), 12, 12, byrow=TRUE)}
	if(Group=="W") {A=matrix(c(0, 0, 0, 1, 0, -1, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 1, 0, -1, 0, 0, 0
				 , 0, 0, 0, 1, 0,  1, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 1, 0,  1, 0, 0, 0
				 , 1, 0, 0, 0, 0,  0, 0, 0,  0, 0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,  0, 0, 0, 0
				 , 0, 0, 1, 0, 0,  0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 1,  0, 0, 0,  0, 0, 0, 0
				 , 0, 0, 0, 0, 0,  0, 0, 1,  0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,  0, 1, 0, 0
				 , 0, 0, 0, 0, 0,  0, 0, 0,  0, 0, 1, 0, 0, 0, 0, 0, 0,  0, 0, 0,  0, 0, 0, 1), 12, 12, byrow=TRUE)}
	if(Group=="S") {A=matrix(c(1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 1,  1, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0
				 , 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 1,  1
				 , 0,  0, 1, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,  0, 0, 1, 0, 0, 0, 0, 0, 0, 0,  0
				 , 0,  0, 0, 0, 1, 0, 0, 0, 0, 0, 0,  0, 0,  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0
				 , 0,  0, 0, 0, 0, 0, 1, 0, 0, 0, 0,  0, 0,  0, 0, 0, 0, 0, 0, 1, 0, 0, 0,  0
				 , 0,  0, 0, 0, 0, 0, 0, 0, 1, 0, 0,  0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 1, 0,  0), 12, 12, byrow=TRUE)}
	if(Group=="K") {A=matrix(c(0, 1, -1, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 1, -1, 0, 0, 0, 0, 0, 0
				 , 0, 1,  1, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 1,  1, 0, 0, 0, 0, 0, 0
				 , 1, 0,  0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 1, 0,  0, 0, 0, 0, 0, 0, 0
				 , 0, 0,  0, 0, 0,  0, 1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,  0, 0, 1, 0, 0, 0, 0
				 , 0, 0,  0, 0, 0,  0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0, 0,  0, 0, 0, 0, 1, 0, 0
				 , 0, 0,  0, 0, 0,  0, 0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,  0, 0, 0, 0, 0, 0, 1), 12, 12, byrow=TRUE)}
	if(Group=="M") {A=matrix(c(0, 0, 0, 0, 0, 0, 1, -1, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 1, -1, 0
				 , 0, 0, 0, 0, 0, 0, 1,  1, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 1,  1, 0
				 , 1, 0, 0, 0, 0, 0, 0,  0, 0, 0,  0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0, 0,  0, 0
				 , 0, 0, 1, 0, 0, 0, 0,  0, 0, 0,  0, 0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,  0, 0
				 , 0, 0, 0, 0, 1, 0, 0,  0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 1, 0,  0, 0, 0,  0, 0
				 , 0, 0, 0, 0, 0, 0, 0,  0, 1, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,  0, 1), 12, 12, byrow=TRUE)}
	if(Group=="RY"){A=matrix(c( 1, 0, -1, 0,  0, 0,  0,  1, -1,  0, 0,  0, 0, 0,  0, 1, -1, 0,  0,  0,  0,  1, 0, -1
                                  , 1, 0,  1, 0,  0, 0,  0, -1, -1,  0, 0,  0, 1, 0, -1, 0,  0, 0,  0, -1,  1,  0, 0,  0
				  , 1, 0,  1, 0,  0, 0,  0,  1,  1,  0, 0,  0, 0, 0,  0, 1,  1, 0,  0,  0,  0, -1, 0, -1
				  , 0, 0,  0, 1, -1, 0,  0,  0,  0, -1, 0,  1, 0, 0,  0, 1,  1, 0,  0,  0,  0,  1, 0,  1
				  , 0, 1,  0, 0,  0, 0,  0,  0,  0,  0, 0,  0, 0, 0,  0, 0,  0, 1,  0,  0,  0,  0, 0,  0
				  , 0, 0,  0, 0,  0, 0,  1,  0,  0,  0, 0,  0, 0, 0,  0, 0,  0, 0,  0,  0,  0,  0, 1,  0), 12, 12, byrow=TRUE)}
	if(Group=="SW"){A=matrix(c(1, -1, 0, 0, 0,  0,  0, 0,  0, 0,  1, -1, 0,  0, 0, 1, 0, -1,  1, 0, -1, 0,  0,  0
				 , 1,  1, 0, 0, 0,  0,  0, 0,  0, 0, -1, -1, 1, -1, 0, 0, 0,  0,  0, 0,  0, 0, -1,  1
				 , 1,  1, 0, 0, 0,  0,  0, 0,  0, 0,  1,  1, 0,  0, 0, 1, 0,  1, -1, 0, -1, 0,  0,  0
				 , 0,  0, 0, 1, 0, -1, -1, 0,  1, 0,  0,  0, 0,  0, 0, 1, 0,  1,  1, 0,  1, 0,  0,  0
				 , 0,  0, 1, 0, 0,  0,  0, 0,  0, 0,  0,  0, 0,  0, 0, 0, 1,  0,  0, 0,  0, 0,  0,  0
				 , 0,  0, 0, 0, 0,  0,  0, 1,  0, 0,  0,  0, 0,  0, 0, 0, 0,  0,  0, 0,  0, 1,  0,  0), 12, 12, byrow=TRUE)}
	if(Group=="KM"){A=matrix(c(0, 1, -1, 0,  1, -1, 0,  0, 0,  0,  0, 0, 0, 0,  0, 0,  0,  0, 1, -1, 0,  1, -1, 0
				 , 0, 1,  1, 0, -1, -1, 0,  0, 0,  0,  0, 0, 0, 1, -1, 0, -1,  1, 0,  0, 0,  0,  0, 0
				 , 0, 1,  1, 0,  1,  1, 0,  0, 0,  0,  0, 0, 0, 0,  0, 0,  0,  0, 1,  1, 0, -1, -1, 0
				 , 0, 0,  0, 0,  0,  0, 1, -1, 0, -1,  1, 0, 0, 0,  0, 0,  0,  0, 1,  1, 0,  1,  1, 0
				 , 1, 0,  0, 0,  0,  0, 0,  0, 0,  0,  0, 0, 0, 0,  0, 1,  0,  0, 0,  0, 0,  0,  0, 0
				 , 0, 0,  0, 0,  0,  0, 0,  0, 1,  0,  0, 0, 0, 0,  0, 0,  0,  0, 0,  0, 0,  0,  0, 1), 12, 12, byrow=TRUE)}
	rx    <- parameters[1:freePar]
	y     <- c(0, 0, (A %*% rx)[-(1:2), 1])
	r.hat <- c(solve(A)%*%y)
	r.hat
}

############################################
# Gram.Schmidt Matrix
# This one works for any number of vectors
############################################
gram.schmidt <- function(v){
	r     = dim(v)[1]
	u     = matrix(0, dim(v)[1], dim(v)[2])
	u[1,] = v[1, ]
	for(k in 2:r){
		proj = 0
		vk   = v[k, ]
		for(j in 1:(k-1)){
			proj.j = c((u[j,]%*%vk)/(u[j,]%*%u[j, ]))*u[j,]
			proj   = proj+proj.j
		}
		u[k, ] = vk - proj
	}
	u			  
}


####################################################
#						   #
# Conditional transition matrices		   #
#						   #
####################################################

#################################################
# Reversible conditional matrices	    	#
#						#
# Input Global parameters			#
# cStat		c				#
# freePar	s				#
# numPar	l				#
# Group		Group				#
#################################################
RevPmat = function(parvec){

	if(Group=="N") {
		svec = parvec[1:freePar]
	} else {
		svec  = Gy.to.S(parvec, Group)
	}

	pivec = parvec[(freePar+1):(freePar+cStat)]
	pivec = pivec/sum(pivec)
	tval  = parvec[numPar]

	PI    = diag(pivec)
	Smat  = matrix(0, cStat, cStat)
	Smat[lower.tri(Smat)] = svec
	Smat  = Smat + t(Smat)

	Rmat  = Smat%*%PI
	Rmat  = Rmat - diag(apply(Rmat, 1, sum))
	eig    = eigen( sqrt(PI)%*%Rmat%*%solve(sqrt(PI)) )
	lambda = eig$values
	L      = eig$vectors
	
	solve(sqrt(PI))%*%L%*%diag(exp(lambda*tval))%*%t(L)%*%sqrt(PI)
}

#################################################
# Non-Reversible conditional matrices	    	#
#						#
#						#
# Input Global parameters			#
# cStat		c				#
# freePar	r				#
# numPar	l				#
#################################################

NoRevPmat = function(parvec){

	if(Group=="N"){
		rvec  = parvec[1:freePar]
	} else {
		rvec  = Gy.to.R(parvec, Group)
	}

	tval  = parvec[numPar]

	Rmat  = rbind(0, matrix(rvec, cStat, cStat-1))
	Rmat  = matrix(c(Rmat, 0), cStat, cStat)

	Rmat  = Rmat - diag(apply(Rmat, 1, sum))

	eig    = eigen(Rmat)
	lambda = eig$values
	M      = eig$vectors
	
	Re(M%*%diag(exp(lambda*tval))%*%solve(M))
}

#########################################################################
# Getting the stationary distribution for non symmetric PiR matrices	#
#
# Input Global parameters
# cStat		c
# freePar	r
#########################################################################

GetStatDist <- function(parvec){
	rvec  = parvec[1:freePar]

	Rmat  = rbind(0, matrix(rvec, cStat, cStat-1))
	Rmat  = matrix(c(Rmat, 0), cStat, cStat)
	Rmat  = Rmat - diag(apply(Rmat, 1, sum))

	Eig   = eigen(Rmat)
	Val0  = match(0, round(Eig$values, 5))
	pix   = Re(solve(Eig$vectors)[Val0, ])

	pix/sum(pix)
}

#########################################################################
#						   			#
# Conditional transition matrix derivatives with respect to the 	#
#	edge lengths				   			#
#						   			#
#########################################################################

#####################################################
# Reversible conditional matrix derivatives	    #
#
#
# Input Global parameters
# cStat		c
# freePar	s
# numPar	l
#####################################################
RevDeriv = function(parvec){

	if(Group=="N") {
		svec = parvec[1:freePar]
	} else {
		svec  = Gy.to.S(parvec, Group)
	}

	pivec = parvec[(freePar+1):(freePar+cStat)]
	pivec = pivec/sum(pivec)
	tval  = parvec[numPar]

	PI    = diag(pivec)
	Smat  = matrix(0, cStat, cStat)
	Smat[lower.tri(Smat)] = svec
	Smat  = Smat + t(Smat)
	Rmat  = Smat%*%PI
	Rmat  = Rmat - diag(apply(Rmat, 1, sum))

	eig    = eigen( sqrt(PI)%*%Rmat%*%solve(sqrt(PI)) )
	lambda = eig$values
	L      = eig$vectors

	solve(sqrt(PI))%*%L%*%diag(lambda)%*%diag(exp(lambda*tval))%*%t(L)%*%sqrt(PI)
}

#####################################################
# Non-Reversible conditional matrix derivatives	    #
#
# Input Global parameters
# cStat		c
# freePar	r
# numPar	l
#####################################################

NoRevDeriv = function(parvec){

	if(Group=="N"){
		rvec  = parvec[1:freePar]
	} else {
		rvec  = Gy.to.R(parvec, Group)
	}

	tval  = parvec[numPar]

	Rmat  = rbind(0, matrix(rvec, cStat, cStat-1))
	Rmat  = matrix(c(Rmat, 0), cStat, cStat)
	Rmat  = Rmat - diag(apply(Rmat, 1, sum))

	eig    = eigen(Rmat)
	lambda = eig$values
	M      = eig$vectors

	Re(M%*%diag(lambda)%*%diag(exp(lambda*tval))%*%solve(M))
}


#########################################################
# Rearrange transition matrix vectors	   		#
#							#
# P-matrix vectors come in the form			#
# {p11, p21, ..., pc1, p12, p22, ..., pcc}		#
#							#
# This comand rearranges the P-matrix vectors		#
# to the form						#
# {p11, p12, ..., p1c, p21, p22, ..., pcc}		#
#
# Input Global parameters
# cStat		c
#########################################################

rearrangePvector <- function(Pvec){
	Pmat <- matrix(Pvec, cStat, cStat)
	Pvec <- c(t(Pmat))
	Pvec
}

#########################################################
#							#
# Normalizing the parameters				#
#							#
#########################################################

#################################################
# Symmetric Pi * R matrices			#
#
# Input Global parameters
# cStat		c
# freePar	s
#################################################
NormalizeRevPars <- function(parvec){
	svec  = parvec[1:freePar]
	pivec = parvec[(freePar+1):(freePar+cStat)]
	pivec = pivec/sum(pivec)
	tval  = parvec[numPar]

	PI    = diag(pivec)
	Smat  = matrix(0, cStat, cStat)
	Smat[lower.tri(Smat)] = svec
	Smat  = Smat + t(Smat)

	Rmat  = Smat%*%PI
	Rmat  = Rmat - diag(apply(Rmat, 1, sum))

	PIR	   = PI%*%Rmat
	RateSubst  = -sum(diag(PIR))

	svec	   = svec/RateSubst
	tval	   = tval*RateSubst
	
	c(svec, pivec, tval)
}

#################################################
# Non-symmetric Pi * R matrices			#
#
# Input Global parameters
# cStat		c
# freePar	r
# numPar	l
#################################################

NormalizeNoRevPars <- function(parvec){
	rvec  = parvec[1:freePar]
	tval  = parvec[numPar]

	Rmat  = rbind(0, matrix(rvec, cStat, cStat-1))
	Rmat  = matrix(c(Rmat, 0), cStat, cStat)
	Rmat  = Rmat - diag(apply(Rmat, 1, sum))

	Eig   = eigen(Rmat)
	Val0  = match(0, round(Eig$values, 5))
	pix   = solve(Eig$vectors)[Val0, ]

	pix   = pix/sum(pix)
	PI    = diag(pix)

	PIR	   = PI%*%Rmat
	RateSubst  = -sum(diag(PIR))

	rvec	   = rvec/RateSubst
	tval	   = tval*RateSubst
	
	c(rvec, tval)
}


#########################################################
#							#
# Internal codes of the Maximum Likelihood estimates	#
#							#
#########################################################

#################################################
# Calculating the maximum likelihood parameters	#
#
# Input Global parameters
#################################################
estimateMleTest <- function(fileDir, seqFile, parMatrixFile, NewickTreeFile, OutFile, initialPi, grpSmat, grpPiVector, beta
			   , InvProbSame, pxInv, verbose, isReversible, Group=NULL, numItr=20){


	if(is.null(Group)) Group="N" else Group=Group
	
	sequenceFile <- file.path(fileDir, seqFile)
	sequenceInfo <- getSequencesL(sequenceFile)

	# Global parameters are defined here
	####################################
	Group        <<- Group
	rootVector   <<- initialPi
	nSites       <<- sequenceInfo$n
	cStat        <<- sequenceInfo$c
	kTaxa        <<- sequenceInfo$k
	rMerg        <<- kTaxa-1

	stateSpace   <- sequenceInfo$S
	seqMat       <- sequenceInfo$seqInfo

	newickFile   <- file.path(fileDir, NewickTreeFile)
	mergeMatInfo <- newickFormatToMergeMatrixL(newickFile)
	taxaNames    <- mergeMatInfo$leaves
	mergeMatrix  <<- mergeMatInfo$dMat
	mergeVector  <<- c(t(mergeMatrix))

	seqMat       <- ChangeTaxaOrder(seqMat, taxaNames)

	patternInfo  <- getUniquePatternsL(seqMat, pxInv)
	patternMat   <<- patternInfo$patMat
	 pConst	     <<- kTaxa+1
	 Countings   <<- kTaxa+2

	InvOrder     <<- patternInfo$InvOrder
	rearrPxInv   <<- (1:cStat)[order(InvOrder)]

	Beta0        <<- beta
	Beta	     <<- beta
	alpha        <<- 1-Beta
	parFile      <- file.path(fileDir, parMatrixFile)
	parMatrix    <- as.matrix(read.table(parFile, header=FALSE))#read.csv(parFile, header=FALSE, sep="\t")
	#parMatrix    <- as.matrix(parMatrix)
	edges        <<- 2*kTaxa-2

	if(isReversible) {
		freePar	<<- cStat*(cStat-1)/2
		numPar	<<- freePar+cStat+1
	}else{
		freePar <<- cStat*(cStat-1)
		numPar	<<- freePar+1
	}

	if(dim(parMatrix)[2]!=numPar) stop(paste("Input matrix should have", numPar, "values per row"))
	# Initial parameters were set until here

	# Number of edges to optimize per group of rate matrices.
	#########################################################
	sMatRows     <<- split(1:edges, grpSmat)
	sMatGroups   <<- 1:length(sMatRows)

	# Number of edges to optimize per group of pi-vectors.
	#########################################################
	piVecRows    <<- split(1:edges, grpPiVector[-1])
	piVecGroups  <<- 1:length(piVecRows)
	rootVecRow   <<- grpPiVector[1]

	# Loglikelihood given the input parameters
	Loglikelihood <- mleLump(parMatrix    = parMatrix
				, isReversible = isReversible
				, InvProbSame  = InvProbSame
				, numItr       = numItr
				, displayInfo  = verbose)

	write.info(fileDir	= fileDir
		, OutFile	= OutFile
		, taxaNames	= taxaNames
		, stateSpace    = stateSpace
		, grpSmat	= grpSmat
		, grpPiVector	= grpPiVector
		, parMatrix	= parMatrix
		, initialPi	= initialPi
		, pxInv		= pxInv
		, numItr	= Loglikelihood$totalItr
		, outLogL	= Loglikelihood$currentLogL
		, outUncLogL    = Loglikelihood$unconsLogL
		, outparMatrix  = round(Loglikelihood$ParMat, 6)
		, outRootVector = Loglikelihood$rootVector
		, outBeta	= Loglikelihood$Beta
		, outpxInv	= Loglikelihood$pxInv)

	Loglikelihood$currentLogL

}


#########################################################################
#									#
# Write information							#
# This command writes the input and output information in a Vivek's 	#
# template								#
#
#
# Input global parameters
#	kTaxa		k
#	rMerg		rM
#	mergeMatrix
#	Beta0		Initial value of Beta
#	edges		ed
#	numPar		l
#########################################################################

write.info <- function(fileDir, OutFile, taxaNames, stateSpace, grpSmat, grpPiVector, parMatrix
		   , initialPi, pxInv, numItr, outLogL, outUncLogL, outparMatrix, outRootVector
		   , outBeta, outpxInv){

	inOutFile <- file.path(fileDir, OutFile)
		write("Input Tree", inOutFile)
		write(rMerg, inOutFile, append = TRUE)
		write(t(mergeMatrix), inOutFile, ncolumns=2, append=TRUE, sep="\t")
		write(taxaNames, inOutFile, ncolumns=kTaxa, append=TRUE, sep="\t")

		write(" ", inOutFile, append = TRUE)
		write("Input Parameters", inOutFile, append = TRUE)

		write("Edges with the same S-matrix", inOutFile, append = TRUE)
		write(grpSmat, inOutFile, ncolumns=edges, append=TRUE, sep="\t")

		write("Root vector and edges with the same F-vector", inOutFile, append = TRUE)
		write(grpPiVector, inOutFile, ncolumns=edges+1, append=TRUE, sep="\t")

		write("Parameters for S-matrix, F-vector and edge length", inOutFile, append = TRUE)
		write(t(parMatrix), inOutFile, ncolumns=numPar, append=TRUE, sep="\t")

		write("Root-vector", inOutFile, append = TRUE)
		write(initialPi, inOutFile, ncolumns=cStat, append=TRUE, sep="\t")

		write("Beta", inOutFile, append=TRUE)
		write(Beta0, inOutFile, append=TRUE)

		write(paste("Prob of the elements of", paste(stateSpace, collapse=", "), "among invariant sites"), inOutFile, append = TRUE)
		write(pxInv, inOutFile, ncolumns=cStat, append=TRUE, sep="\t")

		write("Number of iterations", inOutFile, append=TRUE)
		write(numItr, inOutFile, append = TRUE)
		write(" ", inOutFile, append = TRUE)
		# All input information has been written so far.

		# New Info starts from here on
		write("Output Parameters", inOutFile, append=TRUE)
		write("Log Likelihood and Unconstrained Log Likelihood", inOutFile, append=TRUE)
		write(paste(outLogL, outUncLogL, sep="\t"), inOutFile, append=TRUE)

		write("Parameters for S-matrix, F-vector and edge length", inOutFile, append = TRUE)
		write(t(outparMatrix), inOutFile, ncolumns=numPar, append=TRUE, sep="\t")

		write("Root-vector", inOutFile, append = TRUE)
		write(outRootVector, inOutFile, ncolumns=cStat, append=TRUE, sep="\t")
		
		write("Beta", inOutFile, append=TRUE)
		write(outBeta, inOutFile, append=TRUE)	

		write(paste("Prob of the elements of", paste(stateSpace, collapse=", "), "among invariant sites"), inOutFile, append = TRUE)
		write(outpxInv, inOutFile, ncolumns=cStat, append=TRUE, sep="\t")		

}


#########################################################################################################################################
#																	#
# Obtaining the maximum likelihood estimates of the parameters										#
#																	#
# Global Input Parameters														#
#	Group		Group
# 	nSites		n														#
#	kTaxa		k														#
#	pConst		k+1 (Index in the unique pattern Vector that indicates the likelihood of the site given that is Invariant)	#
#	Countings	k+2 (Index in the unique pattern Vector that indicates where the countings are)					#
#	cStat		c														#
#	rMerg		rM  rows of the merge Matrix											#
#	mergeMatrix															#
#	mergeVector															#
#	patternMat	Matrix with the unique patterns 										#
#			(it also contains the frequency for each site, and the probability of the invariant sites)			#
#	InvOrder 	Real order in which the constant sites will appear in patternMat (unique pattern matrix). This matters if	#
#			cStat > 9 beacuse the invariant sites in patternMat will have an order: 1, 10, 11..., 2, ... cStat		#
#	rearrPxInv	Order in which the invariant site probabilities should be presented in the output file (1, 2, ..., cStat)	#
#	alpha
#	Beta
#	rootVector
#	edges		ed
#	freePar		s or r
#	numPar		l
#
#	Number of edges to optimize for the rate matrices. Now these are global output parameters as well
#	sMatRows     <<- split(1:edges, grpSmat)
#	sMatGroups   <<- 1:length(sMatRows)
#
#	piVecRows    <<- split(1:edges, grpPiVector[-1])
#	piVecGroups  <<- 1:length(piVecRows)
#	rootVecRow   <<- grpPiVector[1]
#########################################################################################################################################

mleLump <- function(parMatrix, isReversible, InvProbSame, numItr, displayInfo){

	## Compute unconstrained log-likelihood
	counts	   = patternMat[Countings,]
	unconsLogL = sum(counts * log(counts/nSites))

	# Verifying the condition of global stationarity to make InvProbSame true
	if(Beta >0 && InvProbSame==TRUE) { if(length(piVecRows)>1 || rootVecRow!= names(piVecRows)[1]) stop("Invariant sites probabilities should be estimated separately") }

	if(isReversible){ P.vectors <- t(apply(parMatrix, 1, RevPmat))
		 } else { P.vectors <- t(apply(parMatrix, 1, NoRevPmat)) }

	# Not rearranging anything in here it is better to have the array like {p11, p21, p31, p41, ... pcc} 
	# as it's going to give me faster calculations.
	####################################################################################################
	P.array		<- array(P.vectors, c(edges, cStat, cStat))
	prevLogL	<- getTotalLogLikelihood(P.array)
	print("Initial Log likelihood"); print(prevLogL)

	ParamMatrix	<<- parMatrix
	Results = optimizeAllParamL(prevLogL	= prevLogL
				 , isReversible	= isReversible
				 , InvProbSame	= InvProbSame
				 , printInfo	= displayInfo
				 , numItr	= numItr
				 , currentItr	= 0)


	newParMatrix	<- Results$ParMat
	newRootVector	<- Results$rootVector
	newAlpha	<- Results$alpha
	newBeta		<- Results$Beta
	newPxInv	<- Results$pxInv
	finalLikelihood	<- Results$likelihood
	totalItr	<- Results$totalItr


	# Obtaining the normalized parameter matrix
	if(isReversible){ newParMatrix <- t(apply(newParMatrix, 1, NormalizeRevPars))
		   }else{ newParMatrix <- t(apply(newParMatrix, 1, NormalizeNoRevPars)) }

	newPxInv	<- newPxInv[order(InvOrder)]

	list(ParMat      = newParMatrix
	   , rootVector  = newRootVector
	   , alpha	 = newAlpha
	   , Beta	 = newBeta
	   , pxInv	 = newPxInv
	   , currentLogL = finalLikelihood
	   , unconsLogL	 = unconsLogL
	   , totalItr    = totalItr)
	
}



#########################################
#					#
# Algorithms related to Optim		#
#					#
#########################################

#########################################################
#							#
#	Optimizing all the parameters			#
#	(Main command)					#
#							#
#							#
# Global input parameters:				#
#	nSites		n				#
#	kTaxa		k				#
#	pConst		k+1				#
#	Countings	k+2				#
#	cStat		c
#	rMerg		rM
#	mergeMatrix
#	patternMat
#	rearrPxInv
#	Beta
#	rootVector
#	edges		ed
#	freePar		s or r
#	numPar		l
#	sMatRows	sMatRows
#	sMatGroups      1:length(sMatRows)
#
#	piVecRows	piVecRows
#	piVecGroups	1:length(piVecRows)
#
#	rootVecRow	rootVecRow
#	ParamMatrix	parameter matrix this is going to be modified through each step
#########################################################
optimizeAllParamL = function(prevLogL, isReversible, InvProbSame, printInfo, numItr, currentItr){

	InitLogL   = prevLogL
	currentItr = currentItr+1

	print(paste("Iteration number", currentItr))
	print("Optimizing all parameters")

	
	#################################
	## Optimize edge lengths	#
	#################################
	initialEdges <- ParamMatrix[, numPar]

	paramNew = constrOptim(theta			= ParamMatrix[,numPar]
				, f			= UpdateEdgeLength
				, grad			= UpdateEdgeGradient
				, ui			= diag(edges)
				, ci			= rep(0, edges)
				, mu			= 1e-04
				, InputParMatrices	= ParamMatrix			## Parameter passed to FUN
				, isReversible		= isReversible			## Parameter passed to FUN
				, method 		= "BFGS"
				, control 		= list(maxit=1, fnscale=-1)	## One iteration + maximization function
				)


	if( paramNew$value == -(1e10) ) {
		ParamMatrix[, numPar] <<- initialEdges
		newLikelihood         <<- InitLogL
	} else {
		ParamMatrix[, numPar] <<- paramNew$par
		newLikelihood         <<- paramNew$value
	}

	if(printInfo) {print("Edges"); print(newLikelihood)}

	#########################################
	## Optimize the Pi-vectors		#
	## (only in reversible R matrices)	#
	#########################################

	if(isReversible){
		sapply(piVecGroups, optimizePIGroup, isReversible, InvProbSame, printInfo)
	}	

	#########################################################
	## Optimize the S-matrices				#
	## (R-matrices for non reversible transition matrices)	#
	#########################################################

	sapply(sMatGroups, optimizeRateGroup, isReversible, InvProbSame, printInfo)

	#############################################################################################
	# You need the NewP.array for the following calculations and for the loop		    #
	#############################################################################################
	if(isReversible==TRUE){ NewPmatrices <- t(apply(ParamMatrix, 1, RevPmat))
		       } else { NewPmatrices <- t(apply(ParamMatrix, 1, NoRevPmat))}

	NewP.array <- array(NewPmatrices, c(edges, cStat, cStat))


	#########################################################################
	## Optimizing the Root-vector 						#
	## This step is only done if the root-Vector was not part of any of the	#
	## Pi-groups (if reversible)						#
	## R-matrices (if not reversible)					#
	#########################################################################

	if(isReversible) rateGroup=piVecGroups else rateGroup=sMatGroups
	if(!(rootVecRow %in% rateGroup)){ 

		updateRoot <- apply(patternMat, 2, UpdateRootVector, P.array = NewP.array)

		updateRoot	 <- apply(updateRoot, 1, sum)
		rootVector	<<- updateRoot/sum(updateRoot)
		newLikelihood	 <- getTotalLogLikelihood(NewP.array)

		if(printInfo) {print("Root-Vector"); print(newLikelihood)}

	}

	###########################################################
	## 							  #
	## Optimizing the invariant sites			  #
	##							  #
	###########################################################

	if(Beta > 0){
		if(InvProbSame) {patternMat[pConst, patternMat[pConst,] > 0] <<- rootVector[order(rearrPxInv)]}
		
		updateInvParameters <- apply(patternMat, 2, UpdateInvariantParam, P.array = NewP.array)

		alpha <<- alpha*sum(updateInvParameters[1, ])/nSites
		Beta  <<- Beta*sum(updateInvParameters[2, ])/nSites

		if(!InvProbSame){
			pxInv 		  <- updateInvParameters[3, updateInvParameters[3, ]>0]
			pxInv 		  <- pxInv/sum(pxInv)
			patternMat[pConst, patternMat[pConst, ] > 0] <<- pxInv
		}

		newLikelihood	<- getTotalLogLikelihood(NewP.array)
		if(printInfo) {print("Invariant parameters"); print(newLikelihood)}
	}

	FinalLogL = newLikelihood

	if(currentItr==numItr || abs(FinalLogL - InitLogL) < 0.001) {
		return(list(ParMat=ParamMatrix, rootVector=rootVector, alpha=alpha, Beta=Beta, pxInv=patternMat[pConst, patternMat[pConst,] >0], likelihood=FinalLogL, totalItr=currentItr))
	}else{
		Recall(newLikelihood, isReversible, InvProbSame, printInfo, numItr, currentItr)

	}
}


#################################################################
#								#
# Codes related to the optimization of Invariant Sites		#
#								#
# Input global parameters					#
#	kTaxa		k					#
#	pConst		k+1					#
#	Countings	k+2					#
#	rMerg		rM					#
#	mergeMatrix						#
#	mergeVector						#
#	alpha							#
#	Beta							#
#	rootVector						#
#################################################################

UpdateInvariantParam <- function(uniquePatVector, P.array){

	patternVector	<- uniquePatVector[ 1:kTaxa ]
	ConstantProb	<- uniquePatVector[ pConst  ]
	patCounts	<- uniquePatVector[Countings]

	alphaDeriv	<- 0
	betaDeriv	<- 0
	pxInvDeriv	<- 0

	if(ConstantProb==0) {
		alphaDeriv	<- patCounts/alpha
	} else {
		siteLikelihood  <- getSiteProbability(uniquePatVector, P.array)
		denominator	<- alpha * siteLikelihood + Beta * ConstantProb
		alphaDeriv	<- patCounts * siteLikelihood / denominator
		betaDeriv	<- patCounts * ConstantProb / denominator
		pxInvDeriv	<- patCounts * Beta * ConstantProb/ denominator
	}
	
	c(alphaDeriv, betaDeriv, pxInvDeriv)
}
	

#################################################################
#								#
# Codes related to the optimization of the root vector		#
#								#
#################################################################

#################################################################################
#										#
# Obtaining the lagrange multiplier to update the pi-vector at the root node  	#
# The code is basically the same as 						#
#					getSiteProbability			#
# The difference is that the site likelihoods are not added at the end		#
#
# Input global parameters
#	rMerg		rM
#	mergeMatrix
#	mergeVector
#	rootVector
#################################################################################

getLagrangeMultiplier <- function(uniquePatVector, P.array){

	index1   	<- mergeMatrix[rMerg, 1]
	index2   	<- mergeMatrix[rMerg, 2]

	ProbMat1 <- P.array[match(index1, mergeVector), , ]
	ProbMat2 <- P.array[match(index2, mergeVector), , ]
	
	diag(rootVector) %*% (getSubTreeLikelihood(ProbMat1, index1, uniquePatVector, P.array)
			    * getSubTreeLikelihood(ProbMat2, index2, uniquePatVector, P.array))
}

#########################################################################################################
#													#
# Optimizing Pi_0											#
# (This function has to be made for a single site, then we can use apply to get it for all sites)	#
#													#
# Input global parameters										#
#													#
#	kTaxa		k										#
#	pConst		k+1										#
#	Countings	k+2										#
#	rMerg		rM										#
#	mergeMatrix
#	mergeVector
#	alpha
#	Beta
#	rootVector
#########################################################################################################

UpdateRootVector <- function(uniquePatVector, P.array){

	patternVector	<- uniquePatVector[1:kTaxa]
	ConstantProb	<- uniquePatVector[pConst]
	patCounts	<- uniquePatVector[Countings]

	temporals	<- getLagrangeMultiplier(patternVector, P.array)
	
	if(ConstantProb==0) { updateRoot = patCounts * temporals/sum(temporals)}
		       else { updateRoot = patCounts * temporals * alpha / (alpha*sum(temporals) + Beta*ConstantProb)}

	updateRoot
}

######################################################### 
# Codes related to the optimization of Pi and S		#
#							#
#########################################################

#########################################################################
# Optimizing the rate matrices						#
#									#
# Input Global parameters						#
#	sMatRows							#
#	ParamMatrix							#
#	rootVecRow							#
#	newLikelihood Important for not reversible cases		#
#	Group								#
#									#
# These go inside the UpdateRateMatrixParam				#
#	kTaxa		k						#
#	pConst		k+1						#
#	Countings	k+2						#
#	cStat		c						#
#	rMerg		rM						#
#	mergeMatrix							#
#	mergeVector							#
#	patternMat							#
#	rearrPxInv							#
#	alpha								#
#	Beta								#
#	rootVector							#
#	edges		ed						#
#	freePar		r or s						#
#	numPar		l						#
#									#
#########################################################################

optimizeRateGroup <- function(RateGroup, isReversible, InvProbSame, printInfo){

	matricesInGroup <- sMatRows[[RateGroup]]
	totalMatInGroup	<- length(matricesInGroup)

	if(rootVecRow==RateGroup) setVector=TRUE else setVector=FALSE

	initialLikelihood <- newLikelihood
	initialParameters <- ParamMatrix[matricesInGroup[1], 1:freePar]

	paramNew = optim(par 			= initialParameters
			, fn 			= UpdateRateMatrixParam
			, gr 			= NULL
			, paramType		= 1				## Parameter passed to FUN
			, rateMatrixRows	= matricesInGroup		## Parameter passed to FUN
			, edSubIndx		= totalMatInGroup		## Parameter passed to FUN
			, isReversible		= isReversible			## Parameter passed to FUN
			, setRootVector		= setVector			## Parameter passed to FUN
			, invSitesEqualVarSites	= InvProbSame			## Parameter passed to FUN
			, method = "L-BFGS-B"			
			, control = list(maxit=1, fnscale=-1)			## One iteration + maximization function
			, lower = rep(10^(-5), freePar)				## Each param >= 10^(-5)
			)
		
	if( paramNew$value == -(1e10) ) {
		newParameters  <<- initialParameters
		newLikelihood  <<- initialLikelihood
	} else {
		newParameters <<- paramNew$par
		newLikelihood <<- paramNew$value
	}

	ParamMatrix[matricesInGroup, 1:freePar] <<- matrix(rep(newParameters, totalMatInGroup), totalMatInGroup, freePar, byrow=TRUE)
	
	if(totalMatInGroup==1){
		if(Group!="N"){
		if(isReversible){
			ParamMatrix[matricesInGroup, 1:freePar] <<- Gy.to.S(ParamMatrix[matricesInGroup, ],Group = Group)
		} else {
			ParamMatrix[matricesInGroup, 1:freePar] <<- Gy.to.R(ParamMatrix[matricesInGroup, ],  Group = Group)
		}
	}else{


		if(Group!="N"){
			if(isReversible){
				ParamMatrix[matricesInGroup, 1:freePar] <<- t(apply(ParamMatrix[matricesInGroup, ], 1, Gy.to.S, Group = Group))
			} else {
				ParamMatrix[matricesInGroup, 1:freePar] <<- t(apply(ParamMatrix[matricesInGroup, ], 1, Gy.to.R, Group = Group))
				}
			}
		}
	}
	#ParamMatrix<<-ParamMatrix


	if(!isReversible) {
		if(setVector) {
			rootVector <<- GetStatDist(newParameters)
			if(InvProbSame) {patternMat[pConst, patternMat[pConst,] >0] <<- rootVector[order(rearrPxInv)]}
		}
	}

	if(printInfo) {
		if(isReversible) print(paste("S-matrix group", RateGroup)) else print(paste("R-matrix group", RateGroup))
		print(newLikelihood)
	}

	return(0)

}

#########################################################################
# Optimizing the Pi- vectors matrices					#
#									#
# Input Global parameters						#
#	piVecRows							#
#	ParamMatrix							#
#	rootVecRow							#
#									#
# These go inside the UpdateRateMatrixParam				#
#	kTaxa		k						#
#	pConst		k+1						#
#	Countings	k+2						#
#	cStat		c						#
#	rMerg		rM						#
#	mergeMatrix							#
#	mergeVector							#
#	patternMat							#
#	rearrPxInv							#
#	alpha								#
#	Beta								#
#	rootVector							#
#	edges		ed						#
#	freePar		r or s						#
#	numPar		l						#
#									#
#########################################################################

optimizePIGroup <- function(PiGroup, isReversible, InvProbSame, printInfo){

	pisInGroup	<- piVecRows[[PiGroup]]
	totalPisInGroup	<- length(pisInGroup)

	if(rootVecRow==PiGroup) setVector=TRUE else setVector=FALSE
						
	paramNew = optim(par 			= ParamMatrix[pisInGroup[1], (numPar-cStat):(numPar-1)]
			, fn 			= UpdateRateMatrixParam
			, gr 			= NULL
			, paramType		= 2				## Parameter passed to FUN
			, rateMatrixRows	= pisInGroup			## Parameter passed to FUN
			, edSubIndx		= totalPisInGroup		## Parameter passed to FUN
			, isReversible		= isReversible			## Parameter passed to FUN
			, setRootVector		= setVector			## Parameter passed to FUN
			, invSitesEqualVarSites	= InvProbSame			## Parameter passed to FUN
			, method 		= "L-BFGS-B"
			, control 		= list(maxit=1, fnscale=-1)	## One iteration + maximization function
			, lower 		= rep(10^(-5), cStat)		## Each param >= 10^(-5)
			)

	paramNew$par <- paramNew$par/sum(paramNew$par)
	ParamMatrix[pisInGroup, (numPar-cStat):(numPar-1)] <<- matrix(rep(paramNew$par, totalPisInGroup), totalPisInGroup, cStat, byrow=TRUE)
	if(setVector) { 
		rootVector      <<- paramNew$par
		if(InvProbSame) {patternMat[pConst, patternMat[pConst,] >0] <<- paramNew$par[order(rearrPxInv)]}
	}

	newLikelihood <<- paramNew$value
	if(printInfo) print(paste("Pi-vector group", PiGroup)); print(newLikelihood)

	return(0)
}


#################################################################
#								#
# Codes related to the optimization of the rate matrices	#
#								#
# Input Global parameters					#
#								#
#	ParamMatrix	InputParMatrices			#
#								#
#	kTaxa		k					#
#	pConst		k+1					#
#	Countings	k+2					#
#	cStat		c					#
#	rMerg		rM					#
#	mergeMatrix						#
#	mergeVector						#
#	patternMat						#
#	rearrPxInv						#
#	alpha							#
#	Beta							#
#	rootVector						#
#	edges		ed					#
#	freePar		r or s					#
#	numPar		l					#
#								#
#################################################################

UpdateRateMatrixParam <- function(parameterSet, paramType, rateMatrixRows, edSubIndx, isReversible, setRootVector, invSitesEqualVarSites){

	if(!isReversible){
		########################
		# Non-reversible cases #
		########################
		ParamMatrix[rateMatrixRows, 1:freePar]		<<- matrix(rep(parameterSet, edSubIndx), edSubIndx, freePar, byrow=TRUE)
		P.matrices					 <- t(apply(ParamMatrix, 1, NoRevPmat))
		PmatArray					 <- array(P.matrices, c(edges, cStat, cStat))
		if(setRootVector) {
			rootVector <<- GetStatDist(ParamMatrix[rateMatrixRows[1],])
			if(invSitesEqualVarSites) {patternMat[pConst, patternMat[pConst,] > 0] <<- rootVector[order(rearrPxInv)]}
		} 
	}else{
		####################
		# Reversible cases #
		####################
		if(paramType==1){
			#######################
			# S matrix parameters #
			#######################
			ParamMatrix[rateMatrixRows, 1:freePar]		<<- matrix(rep(parameterSet, edSubIndx), edSubIndx, freePar, byrow=TRUE)
			P.matrices 					 <- t(apply(ParamMatrix, 1, RevPmat))
			PmatArray					 <- array(P.matrices, c(edges, cStat, cStat))
		}
		if(paramType==2){
			#########################
			# Pi vector parameters  #
			#########################
			parameterSet							 <- parameterSet/sum(parameterSet)
			ParamMatrix[rateMatrixRows, (freePar+1):(freePar+cStat)]	<<- matrix(rep(parameterSet, edSubIndx), edSubIndx, cStat, byrow=TRUE)
			P.matrices 							 <- t(apply(ParamMatrix, 1, RevPmat))
			PmatArray							 <- array(P.matrices, c(edges, cStat, cStat))

			if(setRootVector) {
				rootVector <<- parameterSet
				if(invSitesEqualVarSites) {patternMat[pConst, patternMat[pConst,] > 0] <<- parameterSet[order(rearrPxInv)]}
			} 
		}
	}

	if(any(parameterSet < 1e-04 )) -(1e10) else{

		if(any(P.matrices <= 0)) -(1e10) else  getTotalLogLikelihood(PmatArray)

	}
}


#################################################################
#								#
# Codes related to the optimization of the edge lengths		#
#								#
#################################################################

	#################################################################################
	# This function calculates the edge lengths from a matrix of parameters		#
	#										#
	# Input global parameters							#
	#	kTaxa		k							#
	#	cStat		c
	#	rMerg		rM
	#	mergeMatrix
	#	mergeVector
	#	patternMat
	#	alpha
	#	Beta
	#	rootVector
	#	edges		ed
	#	freePar		r or s
	#	numPar		l
	#################################################################################
	UpdateEdgeLength <- function(EdgeVector, InputParMatrices, isReversible){
	
		InputParMatrices[, numPar] <- EdgeVector

		if(isReversible==TRUE){ P.matrices <- t(apply(InputParMatrices, 1,   RevPmat))
			       } else { P.matrices <- t(apply(InputParMatrices, 1, NoRevPmat))}

		PmatArray	<- array(P.matrices, c(edges, cStat, cStat))

		if(any(EdgeVector < 10^(-5) )) -(1e10) else getTotalLogLikelihood(PmatArray)
	}

	#################################################################################
	# This function is the gradient of UpdateEdgeLength				#
	#										#
	# Input global parameters							#
	#	kTaxa		k							#
	# 	pConst		k+1							#
	#	Countings	k+2							#
	#	cStat		c
	#	rMerg		rM
	#	mergeMatrix
	#	mergeVector
	#	patternMat
	#	alpha
	#	Beta
	#	rootVector
	#	edges		ed
	#	freePar		r or s
	#	numPar		l
	#################################################################################
	UpdateEdgeGradient <- function(EdgeVector, InputParMatrices, isReversible){
	
		InputParMatrices[, numPar] <- EdgeVector

		# Calculating the transition matrices without any partial derivative
		# Calculating the partial derivatives for each transition matrix
		####################################################################
		if(isReversible==TRUE){ P.matrices <- t(apply(InputParMatrices, 1, RevPmat))
					P.partials <- t(apply(InputParMatrices, 1, RevDeriv))
			       } else { P.matrices <- t(apply(InputParMatrices, 1, NoRevPmat))
					P.partials <- t(apply(InputParMatrices, 1, NoRevDeriv))
		}

		PmatArray		<- array(P.matrices, c(edges, cStat, cStat))
		PderArray		<- array(P.partials, c(edges, cStat, cStat))
		gradVec			<- rep(0, edges)


		siteLikelihood		<- apply(patternMat[1:kTaxa, ], 2, getSiteProbability, P.array = PmatArray)

		# P(x_h|inv) = patternMat[pConst,]
		denominator		<- alpha*siteLikelihood + Beta*patternMat[pConst,]
		DepatternMat		<- rbind(patternMat, siteLikelihood, denominator)

		for(i in 1:edges){

			PgradArray	<- PmatArray
			PgradArray[i,,]	<- PderArray[i,,]

			gradVec_i	<- apply(DepatternMat, 2, getPartialDerivLogL, DerivP.array = PgradArray)  # PgradArray is the one that changes for each partial derivative
			gradVec[i]	<- sum(gradVec_i)
		}
		gradVec
	}

	#################################################################################
	#										#
	# Obtain the partial derivatives of the loglikelihood with respect to  		#
	# the edge lengths								#
	#										#
	#	kTaxa		k							#
	#	pConst		k+1
	#	Countings	k+2
	#	rMerg		rM
	#	mergeMatrix
	#	mergeVector
	#	alpha
	#	rootVector
	#################################################################################
	
	getPartialDerivLogL <- function(uniquePatVector, DerivP.array){

		patternVector	<- uniquePatVector[ 1:kTaxa ]
		ConstantProb	<- uniquePatVector[ pConst  ]
		patCounts	<- uniquePatVector[Countings]
		siteLikelihood	<- uniquePatVector[kTaxa + 3]
		denominator     <- uniquePatVector[kTaxa + 4]

		derivLikelihood  <- getSiteProbability(patternVector, DerivP.array)

		if(ConstantProb==0) {partialLogL = patCounts*derivLikelihood/siteLikelihood
			     } else {partialLogL = patCounts*alpha*derivLikelihood/denominator}

		partialLogL

	}


#########################################################################################
#											#
#	Obtain the total log-likelihood for a phylogenetic tree 			#
#											#
# Input global parameters								#
#	kTaxa		k								#
#	pConst		k+1								#
#	Countings	k+2								#
#	rMerg		rM								#
#	mergeMatrix									#
#	mergeVector									#
#	patternMat	<< The function is applied on this matrix given the P.array	#
#	alpha
#	Beta
#	rootVector
#########################################################################
							
getTotalLogLikelihood <- function(P.arrayMatrix) {
	sum(apply(patternMat, 2, getSiteLikelihood, P.arrayMatrix = P.arrayMatrix))
}

#################################################
#						#
# Obtaining the Log likelihood for a given site	#
# 						#
# Input global parameters			#
#	kTaxa		k			#
#	pConst		k+1			#
#	Countings	k+2			#
#	rMerg		rM			#
#	mergeMatrix				#
#	mergeVector				#
#	alpha					#
#	Beta					#
#	rootVector				#
#################################################

getSiteLikelihood <- function(uniquePatVector, P.arrayMatrix){

	patternVector	<- uniquePatVector[1:kTaxa]
	ConstantProb	<- uniquePatVector[pConst]
	patCounts	<- uniquePatVector[Countings]

	siteLikelihood  <- getSiteProbability(uniquePatVector, P.arrayMatrix)
	
	if(ConstantProb==0) {logLikelihood = patCounts*log(alpha * siteLikelihood)
		     } else {logLikelihood = patCounts*log(alpha * siteLikelihood + Beta * ConstantProb)}
	logLikelihood
}


#########################################################################
#									#
# Obtaining the probability P(x_h|theta) for a given variant site	#
# 									#
# Input global parameters
# 	rMerg		rM
#	mergeMatrix
#	mergeVector
#	rootVector
#########################################################################

getSiteProbability <- function(uniquePatVector, P.array){

	index1   	<- mergeMatrix[rMerg, 1]
	index2   	<- mergeMatrix[rMerg, 2]

	ProbMat1 <- P.array[match(index1, mergeVector), , ]
	ProbMat2 <- P.array[match(index2, mergeVector), , ]
	
	rootVector %*% (getSubTreeLikelihood(ProbMat1, index1, uniquePatVector, P.array)
		      * getSubTreeLikelihood(ProbMat2, index2, uniquePatVector, P.array))

}

#################################################################################################################
##	Recursive function for traversing the phylogenetic tree							#
##														#
##	Input													#
##		Local												#
##			condProbMatrix: Conditional probability matrix for a specific edge			#
##			mergeMatIndex:  Node number (internal/ terminal)					#
##			site:           Site Index								#
##		Global												#
##			uniquePatMat:   Matrix of aligned sequences with rows corresponding to leaf nodes	#
##					The unique pattern matrix is used here to fast the calculations		#
##														#
##	Output													#
##		Likelihood for a sub-tree									#
##
##
##	Global input parameters
##
##	mergeMatrix
##	mergeVector
#################################################################################################################
							
getSubTreeLikelihood <- function(condProbMatrix, mergeMatIndex, patternVector, P.arrayMatrix) {

	if(mergeMatIndex < 0) {
		base           <- patternVector[abs(mergeMatIndex)]
		tempLikelihood <- condProbMatrix[, base[1]]
	} else {
		index1         <- mergeMatrix[mergeMatIndex, 1]
		index2         <- mergeMatrix[mergeMatIndex, 2]

		prob1	       <- P.arrayMatrix[match(index1, mergeVector), , ]
		prob2	       <- P.arrayMatrix[match(index2, mergeVector), , ]

		Left	       <- Recall(prob1, index1, patternVector, P.arrayMatrix)
		Right	       <- Recall(prob2, index2, patternVector, P.arrayMatrix)
		tempLikelihood <- condProbMatrix %*% (Left * Right)
	}
	tempLikelihood
}


#########################################################################
#									#
# Codes related to extracting information from the Newick tree file	#
# These codes are ALMOST the same as Vivek's MaxR			#
#									#
#########################################################################

#################################################
# From a Newick tree to a Merge matrix		#
#
# Global Input parameters
#	kTaxa	k
#	rMerg	rM
#################################################
newickFormatToMergeMatrixL = function(phylipFile) {

	fileInfo = getModNewickFileL(phylipFile)
	newickFile = fileInfo$newick
	fileLen = length(newickFile)

	distMatrix = matrix(0, rMerg, 2)

	i = 1
	distMatrixRow = 1				## Current row of distMatrix to populate
	inVector = NULL					## Queue of nodes

	while(i <= fileLen) {

		if((newickFile[i] == "(") || (newickFile[i] == ",")) {

		}## Do nothing if symbol = "(" or ","
		else if(newickFile[i] == ")") {
			maxCount = length(inVector)
			node1 = inVector[maxCount]
			node2 = inVector[maxCount-1]
			inVector = inVector[1:(maxCount-2)]
			inVector = c(inVector, as.numeric(distMatrixRow))
			#distMatrix[distMatrixRow, 1] = min(node1, node2)
			#distMatrix[distMatrixRow, 2] = max(node1, node2)
			distMatrix[distMatrixRow, 1] = node2				
			distMatrix[distMatrixRow, 2] = node1
			distMatrixRow = distMatrixRow + 1

		} ## Link two nodes
		else {
			inVector = c(inVector, as.numeric(newickFile[i]))
		} ## Read the leaf node

		i = i + 1

	} ## The entire Newick tree has been read

	return(list(dMat = distMatrix, leaves = fileInfo$leaf))

}

#########################################################################
#	Convert the Newick format file into a vector such that		#
#	(, ")", "," and leaf nodes can be read separately		#
#
# 
# Global input parameters
#
#	kTaxa	k
#########################################################################

getModNewickFileL = function(inFile) {
	pFile = scan(inFile, what="character")

	modFile = gsub(" ", "", pFile)
	tempFile = gsub("\\(", "", modFile)
	tempFile = gsub("\\)", "", tempFile)
	leafNodes = unlist(strsplit(tempFile, ","))

	tempFile = modFile
	#for(i in 1:kTaxa) {
	for(i in 1:length(leafNodes)) {
		tempFile = sub(leafNodes[i], (-i), tempFile)
	}

	tempFile = gsub("\\(", "\\(:", tempFile)
	tempFile = gsub("\\)", ":\\)", tempFile)
	tempFile = gsub(",", ":,:", tempFile)

	treeComponents = unlist(strsplit(tempFile, split=":"))

	return(list(newick = treeComponents, leaf=leafNodes))
}


#########################################################################
##	Change taxa-order in input matrix				#
##									#
##	Input								#
##		seqMat: Input sequence matrix				#
##		distMatOrder: Taxa names in the order in which they	#
##				  occur in the mergeMatrix		#
##									#
##	Output								#
##		A matrix of sequences with the rows corresponding	#
##		to taxa sequences in the correct order			#
#########################################################################
ChangeTaxaOrder <- function(seqMat, distMatOrder) {

	## Remove blank spaces from taxa names
	taxaNames <- gsub("[ ]*", "", rownames(seqMat))

	indexes             <- match(distMatOrder, taxaNames)
	revSeqMat           <- seqMat[indexes, ]
	rownames(revSeqMat) <- distMatOrder
	return(revSeqMat)
}


#########################################################################
#									#
# Codes related to extracting information from the sequence file	#
#									#
#########################################################################

#################################################################
# Get a sequence matrix with numeric alphabet from an alignment	#
#################################################################
getSequencesL <- function(seqFile){

	n          <- as.numeric(scan(seqFile, what="character")[2])
	seqInfo    <- scan(seqFile, what="character", skip=1, sep="\r")
	k          <- length(seqInfo)

	namesTaxa  <- c(substr(seqInfo, 0, 10))
	namesTaxa  <- gsub(" ", "", namesTaxa)

	seqInfo    <- substr(seqInfo, 11, n+10)

	seqMat     <- NULL
	for(i in 1:k){
		seqMat.i = unlist(strsplit(seqInfo[i], NULL))
		seqMat   = rbind(seqMat, seqMat.i)
	}
	row.names(seqMat) <- NULL

	stateSpace <- names(table(seqMat))
	c          <- length(stateSpace)
	numericS   <- 1:c

	seqMat     <- apply(seqMat, 2, DirectRecode, stateSpace=stateSpace, numericS=numericS, c=c)
	row.names(seqMat) <- namesTaxa

	list(n=n, k=k, c=c, S=stateSpace, seqInfo=seqMat)

}

#########################################################################
# Direct recoding							#
#	From an alphabetical state space to a numerical state space	#
#########################################################################
DirectRecode <- function(siteCol, stateSpace, numericS, c){
	for(i in 1:c){
		siteCol = gsub(stateSpace[i], numericS[i], siteCol)
	}
	siteCol
}

#################################################################################
##	Obtain the unique patterns and their count				#
##										#
## Input Global parameters							#
## 	kTaxa	 Number of taxa
##	cStat		c
##										#
#################################################################################
							
getUniquePatternsL <- function(seqMat, pxInv){

	patterns            <- apply(seqMat, 2, paste, collapse=",")

	patternNameAndCount <- table(patterns)
	patternNames        <- names(patternNameAndCount)

	ConstantPatterns    <- apply(matrix(1:cStat, kTaxa, cStat, byrow=TRUE), 2, paste, collapse=",")
	ConstantPattIndex   <- match(ConstantPatterns, patternNames)

	nUniq		    <- length(patternNames)
	isConstant	    <- match(1:nUniq, ConstantPattIndex)
	isConstant[is.na(isConstant)] <- 0
	isConstProb	    <- sapply(isConstant, PxInvinIsConst, pxInv=pxInv)

	isConstOrder	    <- isConstant[isConstant>0] 	# This one will be needed when printing the optimised Invariant sites in a file

	patternMatrix       <- sapply(patternNames, divide, split=",")
	patternCounts       <- as.integer(patternNameAndCount)
	patternMatrix	    <- rbind(patternMatrix, isConstProb, patternCounts)

	rownames(patternMatrix) <- c(rownames(seqMat), "P(const)", "Counts")
	colnames(patternMatrix) <- NULL
	list(patMat = patternMatrix, InvOrder=isConstOrder)
}

#########################################################################
# Spliting the name of a unique pattern into a numeric vector		#
#########################################################################

divide <- function(patternName, split){
	patternName <- as.integer(unlist(strsplit(patternName, split)))
	
}

#################################################################################
# Getting the probabilities of P(x_h|inv) in the isConstant-indicator		#
#################################################################################

PxInvinIsConst <- function(xIndex, pxInv){
	if(xIndex>0) {xIndex = pxInv[xIndex]} else {xIndex = 0}
}


####################################################
#						   #
# Internal codes of the BT generating sequence     #
#						   #
####################################################

#########################################
# Scanning the initial OutputTemplate 	#
#
# Global parameters are defined here
#	kTaxa		k
#	mergeMatrix
#	mergeVector
#	rMerg		rM
#	Beta
#	Pi0
#	TaxaNames
#	edges		ed
#########################################
scanTemplate <- function(Dir, File){
	TempFile  	<-  file.path(Dir, File)
	FileInfo  	<-  scan(TempFile, what="character", sep="\r")
	kTaxa     	<<- as.numeric(FileInfo[2])+1
	edges		<<- 2*kTaxa-2
	rMerg		<<- kTaxa-1
	TaxaNames    	<<- unlist(strsplit(FileInfo[kTaxa+2], split="\t"))
	mergeMatrix	<<- matrix(as.numeric(unlist(strsplit(as.vector(FileInfo[3:(kTaxa+1)]), split="\t"))), rMerg, 2, byrow=TRUE)
	mergeVector	<<- c(t(mergeMatrix))
	ncols     	<-  length(as.numeric(unlist(strsplit(as.vector(FileInfo[ 3*kTaxa+21]), split="\t"))))
	parMatrix    	<<- matrix(as.numeric(unlist(strsplit(as.vector(FileInfo[(3*kTaxa+21):(5*kTaxa+18)]), split="\t"))), edges, ncols, byrow=TRUE) 
	Beta      	<<- as.numeric(FileInfo[5*kTaxa+22])
	ProbInv   	<<-  as.numeric(unlist(strsplit(as.vector(FileInfo[5*kTaxa+24]), split="\t")))
	Pi0       	<<-  as.numeric(unlist(strsplit(as.vector(FileInfo[5*kTaxa+20]), split="\t")))  
	cStat     	<<- length(Pi0)
}


#################################################
# Defining the root-state of a site 		#
# 
# Input global parameters
#	cStat	c
#	Pi0
#	ProbInv
#################################################
getRootElement <- function(isInv){
	x0 = NULL
	if(isInv==0) x0=sample(1:cStat, 1, prob=Pi0) else x0=sample(1:cStat, 1, prob=ProbInv)
	x0
}

#################################################
# Defining the evolutionary process of a site	#
#
# Input global parameters
#	cStat	c
#	kTaxa	k
#	rMerg	rM
#	mergeMatrix
#	mergeVector
#################################################
siteEvolution <- function(indRootSeq, NodeMatrix, P.array){
	isInv       <- indRootSeq[1]
	x0          <- indRootSeq[2]
	siteOutput  <- NULL
	if(isInv==1) siteOutput = rep(x0, kTaxa) else {
		EvolMatrix = matrix(0, rMerg, 2)
		for(i in 1:cStat){
			if(x0==i) {EvolMatrix[rMerg, 1] = sample(cStat, 1, P.array[NodeMatrix[rMerg, 1],,i], replace=TRUE)
				   EvolMatrix[rMerg, 2] = sample(cStat, 1, P.array[NodeMatrix[rMerg, 2],,i], replace=TRUE)
			}
		}
		for(i in rMerg:1){
			if(mergeMatrix[i, 1]>0) {
				EvolMatrix[mergeMatrix[i, 1], 1] = sample(cStat, 1, P.array[NodeMatrix[(mergeMatrix[i,1]),1],,EvolMatrix[i,1]], replace=TRUE)
				EvolMatrix[mergeMatrix[i, 1], 2] = sample(cStat, 1, P.array[NodeMatrix[(mergeMatrix[i,1]),2],,EvolMatrix[i,1]], replace=TRUE)
			}
			if(mergeMatrix[i, 2]>0) {
				EvolMatrix[mergeMatrix[i, 2], 1] = sample(cStat, 1, P.array[NodeMatrix[(mergeMatrix[i,2]),1],,EvolMatrix[i,2]], replace=TRUE)
				EvolMatrix[mergeMatrix[i, 2], 2] = sample(cStat, 1, P.array[NodeMatrix[(mergeMatrix[i,2]),2],,EvolMatrix[i,2]], replace=TRUE)
			}
		}
		siteOutput = c(t(EvolMatrix))
		siteOutput = siteOutput[mergeVector<0]
	}
	siteOutput
}

#################################################
# Generation of a Bootstraped sequence		#
#################################################
generateBootstrap <- function(fileDir, templateFile, outFile, replicates, numSites, stateSpace, isReversible){

	FileInfo  <- scanTemplate(fileDir, templateFile)
	#
	# From previous command I got this info
	#print(	kTaxa)
	#print(	TaxaNames)
	#print(	mergeMatrix)
	#print(	mergeVector)
	#print(	rMerg)
	#print(	parMatrix)
	#print(	Beta)
	#print(	ProbInv)
	#print(	Pi0)
	#print(	cStat)

	if(isReversible==TRUE){ freePar <<- cStat*(cStat-1)/2
				numPar  <<- freePar+cStat+1 
			P.vectors <- t(apply(parMatrix, 1, RevPmat))
		       } else { freePar <<- cStat*(cStat-1)
				numPar  <<- freePar+1
			P.vectors <- t(apply(parMatrix, 1, NoRevPmat))
		}

	P.vectors   <- t(apply(P.vectors, 1, rearrangePvector))
	P.array     <- array(P.vectors, c(edges, cStat, cStat))
	NodeMatrix  <- matrix(1:edges, rMerg, 2, byrow=TRUE)

	TaxaNames   <- paste(TaxaNames, "         ", sep="")
	TaxaNames   <- substr(TaxaNames, 1, 10)
	Initials    <- paste(kTaxa, numSites, sep="\t")

	numS        <- 1:cStat
	if(cStat > 9) numS[1:9]=paste(0, numS[1:9], sep="")
	outDir      <- file.path(fileDir, outFile)
	unlink(outDir)

	for(i in 1:replicates){

		PropInv     <- rbinom(numSites, 1, Beta)
		RootSeq     <- sapply(PropInv, getRootElement)
		RootSeq     <- rbind(PropInv, RootSeq)

		BTSequences <- apply(RootSeq, 2, siteEvolution
				     , NodeMatrix  = NodeMatrix
				     , P.array     = P.array)
	
		BTSequences <- apply(BTSequences, 2, InverseRecode
				     , stateSpace = stateSpace
				     , numericS   = numS)	
	
		BTSequences <- apply(BTSequences, 1, paste, collapse="")
		BTSequences <- paste(TaxaNames, BTSequences, sep="")
		if(replicates==1){
			BTSequences <- c(Initials, BTSequences)
		} else {
			BTSequences <- c(Initials, BTSequences, "")
		}

		sapply(BTSequences, write, file=outDir, append=TRUE)

	}
	list(ProbabilityArray = P.array, mergeMatrix=mergeMatrix)
}


#################################################################
# Inverse recoding						#
#	From a numerical state space to the initial state Space	#
#
# Input global parameters
#	cStat	c
#################################################################

InverseRecode <- function(siteCol, stateSpace, numericS){
	if(cStat <= 9){
		for(i in 1:cStat){
			siteCol = gsub(numericS[i], stateSpace[i], siteCol)
		}
	}else{
		siteCol[siteCol <= 9] = paste(0, siteCol[siteCol <= 9], sep="")
		for(i in 1:cStat){
			siteCol = gsub(numericS[i], stateSpace[i], siteCol)
		}
	}
	siteCol
}