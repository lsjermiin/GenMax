######################################################################################							
##	R functions for maximum-likelihood estimation of rooted phyolgenetic trees
##	for ACGT and recoded versions such as RY, etc.						
##							
##	Authors: Victor Vera-Ruiz, John Robinson (based on MaxR.R of Vivek Jayaswal
##	and John Robinson. Version from 2020)						
##	Version: 1.0.0						
##							
## 	LICENSE:						
##      Copyright (C) <2011>  <Vivek Jayaswal>							
##							
##	This library is free software; you can redistribute it and/or modify it 						
##	under the terms of the GNU Lesser General Public License as published by 						
##	the Free Software Foundation; either version 2.1 of the License, or (at 						
##	your option) any later version.						
##							
##	This library is distributed in the hope that it will be useful, but 						
##	WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 						
##	or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public 						
##	License for more details.						
##							
##	You should have received a copy of the GNU Lesser General Public License 						
##	along with this library; if not, write to the Free Software Foundation Inc., 						
##	59 Temple Place, Suite 330, Boston, MA 02111-1307 USA 						
######################################################################################							

######################################################################################
##v	Recoding a DNA sequence according to the following ambitides
##v			R = {A, G}				
##v			Y = {C, T}				
##v			S = {C, G}				
##v			W = {A, T}				
##v			K = {A, C}				
##v			M = {G, T}				
##v			B = {C, G, T}
##v			D = {A, G, T}
##v			H = {A, C, T}
##v			V = {A, C, G}
##v	Input
##v		Seq: A DNA sequence with a state space of four nucleotides
##v			  {A, C, G, T}
##v		Group: 	  A character command containing any of the following letters
##v			         "R", "Y", "RY",
##V				 "S", "W", "SW",
##V				 "K", "M", "KM",
##V				 "B", "D", "H", "V"
##V	output
##v		A recoded DNA sequence according to the given notation	 
##v  
######################################################################################							

recode = function(Seq, Group){
	if(Group=="R") {Seq=gsub("A", "R", Seq); Seq=gsub("G", "R", Seq)}					
	if(Group=="Y") {Seq=gsub("C", "Y", Seq); Seq=gsub("T", "Y", Seq)}					
	if(Group=="RY"||Group=="YR") {Seq=gsub("A", "R", Seq); Seq=gsub("G", "R", Seq); Seq=gsub("C", "Y", Seq); Seq=gsub("T", "Y", Seq)}
	if(Group=="W") {Seq=gsub("A", "W", Seq); Seq=gsub("T", "W", Seq)}
	if(Group=="S") {Seq=gsub("C", "S", Seq); Seq=gsub("G", "S", Seq)}
	if(Group=="WS"||Group=="SW") {Seq=gsub("A", "W", Seq); Seq=gsub("T", "W", Seq); Seq=gsub("C", "S", Seq); Seq=gsub("G", "S", Seq)}
	if(Group=="K") {Seq=gsub("G", "K", Seq); Seq=gsub("T", "K", Seq)}
	if(Group=="M") {Seq=gsub("A", "M", Seq); Seq=gsub("C", "M", Seq)}
	if(Group=="KM"||Group=="MK") {Seq=gsub("A", "M", Seq); Seq=gsub("C", "M", Seq); Seq=gsub("G", "K", Seq); Seq=gsub("T", "K", Seq)}
	if(Group=="B") {Seq=gsub("C", "B", Seq); Seq=gsub("G", "B", Seq); Seq=gsub("T", "B", Seq)}
	if(Group=="D") {Seq=gsub("A", "D", Seq); Seq=gsub("G", "D", Seq); Seq=gsub("T", "D", Seq)}
	if(Group=="H") {Seq=gsub("A", "H", Seq); Seq=gsub("C", "H", Seq); Seq=gsub("T", "H", Seq)}
	if(Group=="V") {Seq=gsub("A", "V", Seq); Seq=gsub("C", "V", Seq); Seq=gsub("G", "V", Seq)} else {Seq=Seq}
	Seq
}


######################################################################################
##v	Reducing the number of columns of the parameter matrix according to the following 
##v	notation
##v			R = {A, G}				
##v			Y = {C, T}				
##v			S = {C, G}				
##v			W = {A, T}				
##v			K = {A, C}				
##v			M = {G, T}				
##v			B = {C, G, T}
##v			D = {A, G, T}
##v			H = {A, C, T}
##v			V = {A, C, G}
##v	Input
##v		isReversible: T/F Command that specifies wether a matrix is reversible
##v
##v		parMatrix: Matrix of edge parameters. 
##v			   If it is reversible, it has 11 columns that are organized in
##v			   the following way					
##v				cols 1-6  correspond to S-matrix elements		
##v				cols 7-10 correspond to pi-vector		
##v				col  11   corresponds to the edge lngth		
##v			   If it is not reversible, it has 13 columns where each element 
##v			   of each row corresponds to an off-diagonal parameter of the  
##v			   instant mutation matrix of that row (node)
##v			   It has the following pattern: 				   	
##v							_ d g j
##v							a _ h k
##v							b e _ l
##v							c f i _
##v		
##v		Group: 	  A character command containing any of the following letters
##v			         "R", "Y", "RY",
##V				 "S", "W", "SW",
##V				 "K", "M", "KM",
##V				 "B", "D", "H", "V"
##V	output
##v		A printed reduced parameter matrix	 
##v  
######################################################################################							

reducing=function(parMatrix, isReversible, Group){
	if(isReversible==TRUE){
		Ss  = parMatrix[, 1:6]
		Pis = parMatrix[, 7:10]
		ls  = parMatrix[, 11]
		if(Group=="R") {parMatrix = cbind(1, 1, 1, Pis[,2], (Pis[,1]+Pis[,3]), Pis[,4], ls)}
		if(Group=="Y") {parMatrix = cbind(1, 1, 1, Pis[,1], Pis[,3], (Pis[,2]+Pis[,4]), ls)}
		if(Group=="RY"||Group=="YR") {parMatrix = cbind(1, (Pis[,1]+Pis[,3]), (Pis[,2]+Pis[,4]), ls)}
		if(Group=="S") {parMatrix = cbind(1, 1, 1, Pis[,1], (Pis[,2]+Pis[,3]), Pis[,4], ls)}
		if(Group=="W") {parMatrix = cbind(1, 1, 1, Pis[,2], Pis[,3], (Pis[,1]+Pis[,4]), ls)}
		if(Group=="WS"||Group=="SW") {parMatrix = cbind(1, (Pis[,2]+Pis[,3]), (Pis[,1]+Pis[,4]), ls)}
		if(Group=="K") {parMatrix = cbind(1, 1, 1, Pis[,1], Pis[,2], (Pis[,3]+Pis[,4]), ls)}
		if(Group=="M") {parMatrix = cbind(1, 1, 1, Pis[,3], (Pis[,1]+Pis[,2]), Pis[,4], ls)}
		if(Group=="KM"||Group=="MK") {parMatrix = cbind(1, (Pis[,3]+Pis[,4]), (Pis[,1]+Pis[,2]), ls)}
		if(Group=="B") {parMatrix = cbind(1, Pi[,1], 1-Pi[,1], ls)}
		if(Group=="D") {parMatrix = cbind(1, Pi[,2], 1-Pi[,2], ls)}
		if(Group=="H") {parMatrix = cbind(1, Pi[,3], 1-Pi[,3], ls)}
		if(Group=="V") {parMatrix = cbind(1, Pi[,4], 1-Pi[,4], ls)} 
		else {parMatrix=parMatrix}
	}else{ 
		ls  = parMatrix[, 13]
		if(Group=="R") {parMatrix = cbind(1, 1, 1, 1, 1, 1, ls)}
		if(Group=="Y") {parMatrix = cbind(1, 1, 1, 1, 1, 1, ls)}
		if(Group=="RY"||Group=="YR") {parMatrix = cbind(1, 1, ls)}
		if(Group=="S") {parMatrix = cbind(1, 1, 1, 1, 1, 1, ls)}
		if(Group=="W") {parMatrix = cbind(1, 1, 1, 1, 1, 1, ls)}
		if(Group=="WS"||Group=="SW") {parMatrix = cbind(1, 1, ls)}
		if(Group=="K") {parMatrix = cbind(1, 1, 1, 1, 1, 1, ls)}
		if(Group=="M") {parMatrix = cbind(1, 1, 1, 1, 1, 1, ls)}
		if(Group=="KM"||Group=="MK") {parMatrix = cbind(1, 1, ls)}
		if(Group=="B") {parMatrix = cbind(1, 1, ls)}
		if(Group=="D") {parMatrix = cbind(1, 1, ls)}
		if(Group=="H") {parMatrix = cbind(1, 1, ls)}
		if(Group=="V") {parMatrix = cbind(1, 1, ls)}
	}
	parMatrix
}
							
######################################################################################							
##v	Tagging the nucleotides according to a defined partition of {A, C, G, T} 						
##v	Part of this code is based on Vivek's code "getSequences"						
##v							
##v	Input						
##v		fileName: Absolute path to the sequence file					
##v		sequenceLength: Total number of nucleotides in the alignment					
##v		Group: The way the nucleotides are going to be grouped.					
##v		       Possible scenaries:					
##v			R	{{A, G}, C, T}		-> {"C", "R", "T"}	
##v			Y	{A, G, {C, T}}		-> {"A", "G", "Y" }	
##v			RY	{{A, G}, {C, T}}	-> {"R", "Y"}		
##v			S	{A, {C, G}, T}		-> {"A", "S", "T"}	
##v			W	{{A, T}, C, G}		-> {"C", "G", "W"}	
##v			SW	{{A, T}, {C, G}}	-> {"S", "W"}		
##v			K	{A, C, {G, T}}		-> {"A", "C", "K"}	
##v			M	{{A, C}, G, T}		-> {"G", "M", "T"}	
##v			KM	{{A, C}, {G, T}}	-> {"K", "M"}		
##v			B	{A, {C, G, T}}		-> {"A", "B"}	
##v			D	{C, {A, G, T}}		-> {"C", "D"}	
##v			H	{G, {A, C, T}}		-> {"G", "H"}	
##v			V	{T, {A, C, G}}		-> {"T", "V"}	
##v			N	{A, C, G, T}		-> {"A", "C", "G", "T"}	
##v							
##v	Output						
##v		A matrix of aligned sequences with the following symbols according to					
##v		the option that was given in Group					
##v			R = A G				
##v			Y = C T				
##v			S = C G				
##v			W = A T				
##v			K = A C				
##v			M = G T				
##v			B = C G T				
##v			D = A G T				
##v			H = A C T				
##v			V = A C G				
##v							
######################################################################################							

Recoding=function(inOutDir, sequenceFileName, paramFileName, isReversible, Group){
	# Reading the files
	seqFile	  = paste(inOutDir, sequenceFileName, sep="/")
	parFile   = paste(inOutDir,    paramFileName, sep="/")

	# Saving info from the first  and second rows
	FirstRow  = scan(seqFile, what="character", nline=1, sep="\t")
	InfoSeqs  = scan(seqFile, what="character", skip=1, sep="\r")
	# Obtaining the taxon names and DNA sequences
	namesTaxa  <- c(substr(InfoSeqs, 0, 10))
	namesTaxa  <- gsub(" ", "", namesTaxa)
	InfoSeqs=substr(InfoSeqs,11,as.numeric(FirstRow[2])+10)
	InfoSeqs = recode(InfoSeqs, Group=Group)

	#Victors input (requires names to have <10 and spaces)
	#InfoSeqs  = scan(seqFile, what="character", skip=1, sep="")
	#InfoSeqs  = matrix(InfoSeqs, ncol=2, byrow=TRUE)	
	# Obtaining the taxon name
	#namesTaxa = InfoSeqs[,1]
	# Recoding the DNA sequences
	#InfoSeqs = InfoSeqs[, 2]
	#NewSeqFile = gsub(".txt", "Recoded.txt", sequenceFileName)
	#NewSeqFile = paste(inOutDir, NewSeqFile, sep="/")

	# Rewriting a new sequence file with the reduced state-space
	NewSeqFile = gsub(".txt", "Recoded.txt", seqFile)
	
	#JR code to give output with names in first 10 places
	write(paste(FirstRow[1],FirstRow[2],sep="\t"),NewSeqFile,append=TRUE)
	for(i in 1:length(namesTaxa)){
		nodeName = paste(namesTaxa[i], "          ", sep="")
		nodeName = substr(nodeName,1,10)
		write(paste(nodeName,InfoSeqs[i],sep=""),NewSeqFile,append=TRUE)
		}

#Victor's code for output
#	namesTaxa  = c(FirstRow[1], namesTaxa)
#	InfoSeqs   = c(FirstRow[2], InfoSeqs)
#	NewSeqFile = gsub(".txt", "Recoded.txt", sequenceFileName)
#	NewSeqFile = paste(inOutDir, NewSeqFile, sep="/")
#	write(paste(namesTaxa, InfoSeqs), NewSeqFile, sep="\t")

	# Reducing the parameterMatrix
	InfoPara   = read.table(parFile, header=FALSE)
	InfoPara   = as.matrix(InfoPara)  
	InfoPara   = reducing(InfoPara, isReversible, Group)
	columns    = dim(InfoPara)[2]

	# Rewriting a New parameter matrix with the reduced state-space	
	NewParFile = gsub(".txt", "Reduced.txt", paramFileName)
	NewParFile = paste(inOutDir, NewParFile, sep="/")
	write(t(InfoPara), ncolumns=columns, NewParFile, sep="\t")
}


######################################################################################
##v	Read the aligned nucleotide/amino acid sequences from a Phylip file in sequential format. 
##v	In fact, all the base/amino acid sequence corresponding to a particular species have to be on 
##v	the same line as the species name.
##	
##	Input
##v		fileName: Absolute path to the sequence file
##v		It is not necessary to have an input for the length of the sequences since this
##v		information already appears on the Phylip file in sequential format.
##v		BE CAREFUL IN THIS CASE THE NAMES OF THE TAXA CAN ONLY HAVE 9 CHARACTERS LEGTH. CHARACTER NUMBER 10 MUST BE A SPACE
##
##	Output
##v		A matrix of aligned sequences 
##v		let be C={c1, c2, ... , cn} be a state space with n elements.
##v		Then, {c1, c2, ..., cn} = {1, 2, ..., n}
##v		The numbers are assigned according to the alphabetical order of the elements in
##v		the state space
######################################################################################

getSequences = function(fileName) {							#v
	
	# Obtain the length of the sequences						#v
	seqlength=as.numeric(scan(fileName, what="character")[2])			#v Temporal change sep="\t"

	info = scan(file=fileName, what="character", skip=1, sep="\r")			#v note that now sep="". This is difference is fundamental

	# Obtain the taxon name
	namesTaxa  <- c(substr(info, 0, 10))
	namesTaxa  <- gsub(" ", "", namesTaxa)

	info    <- substr(info, 11, seqlength+10)
	numTaxa = length(info)
	seqMat  = NULL

	for(i in 1:numTaxa){								#v
		seqMat.i = unlist(strsplit(info[i], NULL))
		seqMat   = rbind(seqMat, seqMat.i)

	}										#v 
	row.names(seqMat) <- NULL

	StateSpace <<- names(table(seqMat))						#v I have converted the commands "StateSpace" and "NumericStateSpace" into global inputs. I will need to use them both in the command generateBootstraps 
	NumericStateSpace <<- 1:length(StateSpace)					#v I could put "r" NumericStateSpace<<-1:r 
	
	for(h in 1:length(StateSpace)){							#v

		## Convert the State Space of {"c1", "c2", ..., "cn"}			#v
		## into a vector of {1, 2, ..., n}					#v
		seqMat=gsub(StateSpace[h], NumericStateSpace[h], seqMat)		#v
	}										#v
	seqMat = matrix(as.numeric(seqMat), numTaxa, seqlength, byrow=FALSE)		#v
	rownames(seqMat) = namesTaxa							#v

	return(seqMat)
}						
							
							
##############################################################################							
##	Obtain the conditional probability matrix for a non-reversible rate-matrix						
##							
##	Input						
##		Local					
##			paramIndex: Node number
##			
##		Global					
##			paramMatrix: 
##v				BE CAREFUL: In this particular case, the number of colummns
##v				of paramMatrix are equal to:
##v					colsirrev + 1  =  r(r-1)+1					
##v
##v				   Example: for the i-th row of ParamMatrix file.
##v
##v					if r = 4, the i-th row has the following array:
##v						
##v					a b c d e f g h i j k l Edgelength 
##v						The a-l values are going to be organized in 
##v						the R matrix as:
##v							_ d g j
##v							a _ h k
##v							b e _ l
##v							c f i -
##
##							
##	Output						
##		Matrix of conditional probabilities with rows corresponding 					
##		to the parent node and columns corresponding to the child					
##		node					
##############################################################################							
							
getUnrevCondProb = function(paramIndex) {					#v

	paramirrev = paramMatrix[paramIndex, 1:colsirrev]			#v

	#paramirrev = Gy.to.R(paramirrev, Group)					# vforlump

	t = paramMatrix[paramIndex, colsirrev+1]				#v

	R = matrix(paramirrev, r, r-1)						#v
	R = rbind(R, 0)								#v
	R = as.vector(R)							#v
	R = matrix(c(0, R), r, r)						#v
	R = R - diag(apply(R, 1, sum))						#v

	eigR = eigen(R)
	Pt = abs(eigR$vectors %*% diag(exp(eigR$values*t)) %*% solve(eigR$vectors))

	return(Pt)
}
							
							
##############################################################################							
##	Obtain the conditional probability matrix						
##							
##	Input						
##		Local					
##			paramIndex: Node number
##v			Group: The way the nucleotides are going to be grouped				
##		Global					
##			paramMatrix: Matrix of edge parameters such that				
##v					cols 1-s correspond to S-matrix elements		
##v					cols (s+1)-(s+r) correspond to pi-vector		
##v					col  (s+r+1) corresponds to edge length		
##							
##	Output						
##		Matrix of conditional probabilities with rows corresponding 					
##		to the parent node and columns corresponding to the child					
##		node					
##############################################################################							
							
getCondProb = function(paramIndex, isReversible) {				#v						
							
	index = paramIndex						
							
	## If the matrix is not reversible, return the relevant P-matrix						
	if(!isReversible) return(getUnrevCondProb(index))	                #v
							
	##v First s parameters correspond to the S-matrix						
	##v Next r parameters correspond to the Pi-vector						
	##  Last parameter corresponds to the time of evolution						
	sMatrix = getSmatrix(paramMatrix[index, 1:s], paramMatrix[index, (s+1):(s+r)])	#v
	###################################################################
	####The next line is used to replace the previous one in Later GenMaxR.R
	####It seems to be needed when using the LumpTest
	##sMatrix = getSmatrix(Gy.to.S(paramMatrix[index, 1:s], paramMatrix[index, (s+1):(s+r)], Group), paramMatrix[index, (s+1):(s+r)]) #vforlump
							
	## Diagonal matrix corresponding to the pi-vector						
	pi = diag(paramMatrix[index, (s+1):(s+r)])					#v					
							
	## Time is the last parameter						
	t = paramMatrix[index, (s+r+1)]							#v
							
	## Symmetrize the rate matrix before obtaining the eigen values 						
	## and eigen vectors						
	ax = matrix(0, nrow = r, ncol = r)						#v
	fx = matrix(0, nrow = r, ncol = r)						#v
							
	Rx = sMatrix %*% pi						
							
	eigVal = eigen((sqrt(pi)) %*% Rx %*% (solve(sqrt(pi))), symmetric = T)$values						
	eigVector = eigen((sqrt(pi)) %*% Rx %*% (solve(sqrt(pi))), symmetric = T)$vectors						
							
	for(j in 1:r) {									#v
		ax = exp(eigVal[j] * (t)) * eigVector[, j] %*% t(eigVector[, j])					
		fx = fx + ax					
	}						
							
	probMatrix = (solve(sqrt(pi))) %*% fx %*% (sqrt(pi))						
							
	return(probMatrix)						
							
}							
							
#######################
# All codes to optim
########################


							
####################################################							
##	Obtain the S-matrix for a given edge						
##	This code has completely been modified
##							
##	Input						
##v		S: Vector of s elements
##v		s=r*(r-1)/2
## 		piVector: Stationary probability.
##v		piVector has r elements
##							
##	Output						
##		S-matrix					
####################################################							
							
getSmatrix = function(S, piVector) {						#v		
    r=length(piVector)								#v
    pi = diag(piVector)								#v
    sMat=matrix(0, r, r)							#v
    sMat=lower.tri(sMat)							#v
    sMat[lower.tri(sMat)]<-S							#v
    sMat=sMat+t(sMat)								#v
    for(i in 1:r) {								#v
        sMat[i, i] = (- sum((sMat[i,  - i]) %*% (diag(pi)[-i])))/diag(pi)[i]
    }
    return(sMat)
}

							
#############################################################################							
##	Obtain the likelihood for a given site						
##							
##	Input						
##		Local					
##			site: Site number to obtain the relevant pattern 				
##				from the sequenceMatrix			
##		Global					
##			mergeMatrix: Rooted tree topology in matrix format				
##			rootVector: Pi-vector for the root				
##							
##	Output						
##		Likelihood for the site					
#############################################################################							
							
getSiteLikelihood = function(site) {							
							
	## The R function hclust returns a merge matrix						
	## that has the top-most bifurcation of phylogenetic						
	## tree as the last row						
	numRows = dim(mergeMatrix)[1]						
							
	index1 = mergeMatrix[numRows, 1]						
	index2 = mergeMatrix[numRows, 2]						
							
	prob1 = globalProbMatrix[match(index1, mergeVector), , ]						
	prob2 = globalProbMatrix[match(index2, mergeVector), , ]						
							
	siteLikelihood = diag(rootVector) %*% (getSubTreeLikelihood(prob1, index1, site) 						
						* getSubTreeLikelihood(prob2, index2, site))	
							
	return(apply(siteLikelihood, 2, sum))						
}							
							
###############################################################################							
##	Recursive function for traversing the phylogenetic tree						
##							
##	Input						
##		Local					
##			condProbMatrix: Conditional probability matrix for 				
##					a specific edge		
##			mergeMatIndex: Node number (internal/ terminal)				
##			site: Site Index				
##		Global					
##			sequenceMatrix: Matrix of aligned sequences with 				
##					rows corresponding to leaf nodes		
##							
##	Output						
##		Likelihood for a sub-tree					
###############################################################################							
							
getSubTreeLikelihood = function(condProbMatrix, mergeMatIndex, site) {							
							
	if(mergeMatIndex < 0) {						
		base = sequenceMatrix[abs(mergeMatIndex), site]					
		## tempLikelihood = condProbMatrix[, base]					
		tempLikelihood = matrix(condProbMatrix[, base[1]], nrow=r, ncol=1)					
		if(length(site) > 1) {					
			for(i in 2:length(site)) tempLikelihood = cbind(tempLikelihood, condProbMatrix[, base[i]])				
		}					
							
	} ##Leaf node						
	else {						
		index1 = mergeMatrix[mergeMatIndex, 1]					
		index2 = mergeMatrix[mergeMatIndex, 2]					
							
		prob1 = globalProbMatrix[match(index1, mergeVector), , ]					
		prob2 = globalProbMatrix[match(index2, mergeVector), , ]					
		left = Recall(prob1, index1, site)					
		right = Recall(prob2, index2, site)					
							
		#tempLikelihood = condProbMatrix %*% matrix((left * right), ncol=1)					
		tempLikelihood = condProbMatrix %*% (left * right)					
							
	} ## Internal node						
							
	return(tempLikelihood)						
}							

							
#########################################################################							
##	Obtain the total log-likelihood for a phylogenetic tree 						
##							
##	Input						
##		Global					
##			sequenceMatrix: Matrix of aligned sequences				
##					with only the unique patterns		
##			paramMatrix: Matrix of parameters per edge				
##			siteCount: Number of occurences of patterns 				
##				   in sequenceMatrix			
##			cstSitesIndex: Index of constant sites in a determined order				
##							
##	Output						
##		Log-likelihood based on all the sites in the 					
##		alignment					
#########################################################################							
							
getLogLikelihood = function(isReversible) {							
							
	numSites = dim(sequenceMatrix)[2]						
	logLikelihood = 0						
							
	steps = 400						
	numSteps = ceiling(numSites/steps)						
	siteLikelihood = NULL						
							
	for(siteIndex in 1:numSteps) {						
							
		siteStart = (siteIndex - 1) * steps + 1					
		siteEnd = siteStart + steps - 1					
		if(siteEnd > numSites) siteEnd = numSites					
							
		siteLimit = siteStart:siteEnd					
		siteLikelihood = c(siteLikelihood, getSiteLikelihood(siteLimit))					
	}						
							
	for(siteIndex in 1:numSites) {						
		index = match(siteIndex, cstSitesIndex)					
							
		if(is.na(index)) logLikelihood = logLikelihood + siteCount[siteIndex] * log(alpha * siteLikelihood[siteIndex])					
		else logLikelihood = logLikelihood + siteCount[siteIndex] * log(alpha * siteLikelihood[siteIndex] + betaVar * probXGivenInv[index])					
	}						
							
	return(logLikelihood)						
							
}							
							
###################################################################							
##	Convert each column of the input data matrix into  						
##	comma-separated values 						
##							
##	Input						
##		x: Column of the sequence data matrix					
##		numTaxa: Number of taxa					
##							
##	Output						
##		Concatenated string of elements separate by commas					
###################################################################							
							
getSeqVector = function(x, numTaxa) {							
	seq = x[1]						
	for(i in 2:numTaxa) seq = paste(seq, x[i], sep=",")						
	return(seq)						
}							
							
#####################################################################							
##	Obtain the unique patterns and their count						
##							
##	Input						
##		seqData: Sequence data matrix					
##							
##	Output						
##		A list containing two elements:					
##		1. mat: Matrix of unique patterns					
##		2. patCount: Number of occurences of patterns in mat					
##v		3. cstSitesIndex: Index of constant sites in a determined order				
##			
#####################################################################							
							
getUniquePatterns=function(seqData){												#v
	numTaxa=dim(seqData)[1]													#v
	patterns=apply(seqData, 2, getSeqVector, numTaxa)									#v
	patternNameAndCount=table(patterns)											#v
	patternNames=names(patternNameAndCount)											#v
	uniquePatternNumber=length(patternNames)										#v
	patternMat = matrix(0, numTaxa, uniquePatternNumber)									#v
																#v
	sitesIndexTemp=rep(0, r)												#v
	names(sitesIndexTemp)=1:r												#v
																#v
	for(j in 1:uniquePatternNumber){											#v
		patternMat[,j]=as.integer(unlist(strsplit(patternNames[j], ",")))						#v
		for(h in 1:r){													#v
			if(length(grep(as.integer(names(sitesIndexTemp))[h], patternMat[,(j)])) == numTaxa) sitesIndexTemp[h]=j	#v
		}														#v
	}															#v
	rownames(patternMat)=rownames(seqData)											#v
	return(list(mat = patternMat, patCount= as.integer(patternNameAndCount), cstIndex = sitesIndexTemp))			#v
}																#v
							
########################################################################							
##	Obtain the maximum-likelihood estimates of the	parameters					
##							
##	Input						
##		Local					
##			displayInfo: TRUE -> Display log likelihoods				
##				     at each optimization step			
##			computeInvSites: TRUE -> Compute invariant 				
##					sites parameters		
##			numItr: Number of iterations				
##			inOutFile: File name for saving input and output parameters
##							
##		Global					
##			paramMatrix: User-defined matrix of 				
##  				     input parameters per edge			
##			sMat: Group to which a S-matrix belongs				
##			piVect: Group to which a pi-vector belongs 				
##				including root vector			
##	Output						
##		MLE parameters saved in paramMatrix					
########################################################################							
							
mleParam = function(displayInfo = TRUE, computeInvSites=FALSE, isInvProbSame=FALSE
				, numItr=20, inOutFile, isReversible) {					# v

	## Determine the conditional probability matrices						
	numEdges = dim(paramMatrix)[1]
							
	for(i in 1:numEdges) globalProbMatrix[i, , ] <<- getCondProb(i, isReversible)			# v				
							
	prevLogL = getLogLikelihood(isReversible)						
							
	## Obtain the number of edges to optimize for rate matrices						
	sMatRows = split(seq(1:numEdges), sMat)
	piVectRows = split(seq(1:numEdges), piVect[-1])						
	rootVectRow = piVect[1]						
							
	## P(X|inv) = P(X|var) is possible if and only if GTR model over the entire tree						
	## with pi_root = pi_GTR_Model						
	if(isInvProbSame && computeInvSites) {						
		if((length(piVectRows) > 1)					
			|| (rootVectRow != names(piVectRows)[1])) stop("Inv sites prob should be estimated separately")				
	}						
							
	## Save the input parameters in the file						
	write("Input Tree", inOutFile)						
	write(dim(mergeMatrix)[1], inOutFile, append = TRUE)						
	write(t(mergeMatrix), inOutFile, ncolumns = 2, append = TRUE, sep="\t")						
	write(rownames(sequenceMatrix), inOutFile, ncolumns = dim(sequenceMatrix)[1], append = TRUE, sep="\t")						
							
	write(" ", inOutFile, append = TRUE)						
	write("Input Parameters", inOutFile, append = TRUE)						
							
	write("Edges with the same S-matrix", inOutFile, append = TRUE)						
	write(sMat, inOutFile, ncolumns = length(sMat), append = TRUE, sep="\t")						
							
	write("Root vector and edges with the same F-vector", inOutFile, append = TRUE)						
	write(piVect, inOutFile, ncolumns = length(piVect), append = TRUE, sep="\t")						
							
	write("Parameters for S-matrix, F-vector and edge length", inOutFile, append = TRUE)						
	write(t(paramMatrix), inOutFile, ncolumns = dim(paramMatrix)[2], append = TRUE, sep="\t")						
							
	write("Root-vector", inOutFile, append = TRUE)						
	write(rootVector, inOutFile, ncolumns = r, append = TRUE, sep="\t")						#v
							
	write("Beta", inOutFile, append = TRUE)						
	write(betaVar, inOutFile, append = TRUE)						
							
	write("Prob of the elements of the State-Space among invariant sites", inOutFile, append = TRUE)		#v
	write(probXGivenInv, inOutFile, ncolumns = r, append = TRUE, sep="\t")						
							
	write("Number of iterations", inOutFile, append = TRUE)						
	write(numItr, inOutFile, append = TRUE)						
							
	## Compute unconstrained log-likelihood						
	totalSites = sum(siteCount)						
	unconsLogL = sum(siteCount * log(siteCount/totalSites))						
							
	itr = TRUE						
	for(i in 1:numItr) {						
		print(i)					
		print(prevLogL)					
		currentParamMatrix = paramMatrix					
							
		optimizeAllParam(currentParamMatrix					
					, numEdges		
					, sMatRows		
					, piVectRows		
					, rootVectRow		
					, displayInfo		
					, computeInvSites		
					, isInvProbSame		
					, isReversible)		
							
		currentLogL = getLogLikelihood(isReversible)					
		if(abs(currentLogL - prevLogL) < 0.001) itr = FALSE					
		prevLogL = currentLogL					
	}						
							
	## Save the output parameters in the file						
	write(" " , inOutFile, append = TRUE)						
	write("Output Parameters" , inOutFile, append = TRUE)						
							
	write("Log Likelihood and Unconstrained Log Likelihood", inOutFile, append = TRUE)						
	info = paste(currentLogL, unconsLogL, sep = "\t")						
	write(info, inOutFile, append = TRUE)						
							
	write("Parameters for S-matrix, F-vector and edge length", inOutFile, append = TRUE)						
							
	## Obtain the normalized rate matrix and time						
	rateTime = getNormalizedValues(paramMatrix, isReversible)						#v						
	write(t(rateTime), inOutFile, ncolumns = dim(rateTime)[2], append = TRUE, sep="\t")						
							
	write("Root-vector", inOutFile, append = TRUE)						
	write(rootVector, inOutFile, ncolumns = r, append = TRUE, sep="\t")						
							
	write("Beta", inOutFile, append = TRUE)						
	write(betaVar, inOutFile, append = TRUE)						
							
	write("Prob of the elements of the State-Space among invariant sites", inOutFile, append = TRUE)	#v
	write(probXGivenInv, inOutFile, ncolumns = r, append = TRUE, sep="\t")						
							
	return(currentLogL)						
							
}							
							
#############################################################							
##	Update all the parameters - all edge lengths and 						
##	all the unique rate matrices						
##							
##	Input						
##		curParamMat: Current estimate of rate matrices					
##		numEdges: Number of edges in the tree					
##		sMatrixRows: List of edges that have a common 					
##			     rate matrix				
##		piVectorRows: List of edges that have a common 					
##			     pi-vector				
##		rootVectorRow: Group to which the root vector					
##				belongs			
##		printInfo: TRUE -> Display log likelihoods					
##			   at each optimization step				
##		computeInvSites: TRUE -> Compute invariant sites					
##		isInvProbDiff: TRUE -> P(X|inv) <> P(X|var)					
##				   FALSE -> P(X|inv) = P(X|var)			
##							
##							
##	Output						
##		Updated parameters					
#############################################################							
							
optimizeAllParam = function(curParamMat, numEdges, sMatrixRows							
				, piVectorRows, rootVectorRow, printInfo			
				, computeInvSites, isInvProbSame, isReversible) {			
							
	## Parameter optimization is done stepwise using the current 						
	## values - (1) Optimize all the edge lengths, (2) optimize						
	## the rate matrix parameters per edge, and (3) optimize						
	## the root node vector						
							
	print("Optimize all parameters")						
							
	###########################						
	## Optimize edge lengths						
	###########################						
							
	if(isReversible) edgeLengthCol = s+r+1 				#v						
	else edgeLengthCol = colsirrev+1				#v						
							
	paramNew = optim(par = paramMatrix[ ,edgeLengthCol]		## Parameter to be optimized				
			, getUpdatedParam				## Function to be used for parameter optimization
			, gr = NULL				
			, paramType = 1 				## Parameter passed to FUN
			, rateMatrixRows = NULL			## Parameter passed to FUN	
			, setRootVector = NULL			## Parameter passed to FUN	
			, invSitesEqualVarSites = NULL	## Parameter passed to FUN			
			, isRevMatrix = isReversible		## Parameter passed to FUN		
			, method = "L-BFGS-B"				
			, control = list(maxit=1, fnscale=-1)	## One iteration + maximization function			
			, lower = c(rep(10^(-5), numEdges))	## Edge_length >= 10^(-5)			
			)				
							
	paramMatrix[, edgeLengthCol] <<- paramNew$par						
							
	if(printInfo) {						
		print("Edge")					
		print(paramNew$value)					
	}						
							
	##########################						
	## Optimize the S-matrix						
	##########################						
							
	numRateMatrices = length(sMatrixRows)						
							
	for(i in 1:numRateMatrices) {						
							
		rowIndexes = sMatrixRows[[i]]					
							
		if(isReversible) numParam = s					#v
		else numParam = 2*s						#v
							
		paramNew = optim(par = paramMatrix[rowIndexes[1], 1:numParam]	## Parameter to be optimized				
				, fn = getUpdatedParam			## Function to be used for parameter optimization
				, gr = NULL			
				, paramType = 2 			
				, rateMatrixRows = rowIndexes		## Parameter passed to FUN	
				, setRootVector = NULL			## Parameter passed to FUN
				, invSitesEqualVarSites = NULL	## Parameter passed to FUN		
				, isRevMatrix = isReversible		## Parameter passed to FUN	
				, method = "L-BFGS-B"			
				, control = list(maxit=1, fnscale=-1)	## One iteration + maximization function		
				, lower = c(rep(10^(-5), numParam))		## Each param >= 10^(-5)	
				)			
							
		for(rIndex in 1:length(rowIndexes)) paramMatrix[rowIndexes[rIndex], 1:numParam] <<- paramNew$par[1:numParam]					
							
		if(printInfo) {					
			print("S-matrix")				
			print(paramNew$value)				
		}					
							
	} ## All S-matrices have been considered						
							
	######################						
	## Optimize pi-vector						
	######################						
							
	numRateMatrices = length(piVectorRows)						
	piGrp = names(piVectorRows)						
							
	if(isReversible) {						
							
		for(i in 1:numRateMatrices) {					
							
			rowIndexes = piVectorRows[[i]]				
							
			if(rootVectorRow == piGrp[i]) setVector = TRUE				
			else setVector = FALSE				
							
			paramNew = optim(par = paramMatrix[rowIndexes[1], (s+1):(s+r)]	##v Parameter to be optimized			
					, fn = getUpdatedParam		
					, gr = NULL		
					, paramType = 3 		
					, rateMatrixRows = rowIndexes		
					, setRootVector = setVector		
					, invSitesEqualVarSites = isInvProbSame		
					, isRevMatrix = isReversible		
					, method = "L-BFGS-B"		
					, control = list(maxit=1, fnscale=-1)	##  One iteration + maximization function	
					, lower = c(rep(10^(-5), r))		##v Each param >= 10^(-5)
					)		
							
			piRate = paramNew$par/sum(paramNew$par)				
			for(rIndex in 1:length(rowIndexes)) paramMatrix[rowIndexes[rIndex], (s+1):(s+r)] <<- piRate				
			if(setVector) rootVector <<- piRate				
			if(isInvProbSame) {				
				probXGivenInv <<- piRate			
				names(probXGivenInv) <<- NULL			
			}				
							
			if(printInfo) {				
				print("Pi-vector")			
				print(paramNew$value)			
			}				
							
		} ## All pi-vector have been considered					
							
	} ## Pi-vectors are calculated only if the R-matrix is reversible						
							
	###############################						
	## Optimize root node vector						
	###############################						
							
	if(!(rootVectorRow %in% piGrp) || !isReversible) {						
		if(printInfo) print("Root vector")					
		getUpdatedRootVector(isReversible)					
	}						
							
	######################################						
	## Optimize invariant sites parameters						
	######################################						
							
	if(computeInvSites) {						
		if(printInfo) print("Inv sites")					
							
		if(isInvProbSame) {					
			probXGivenInv <<- rootVector				
							
			## Avoid output of the form [1] 2				
			## V7 				
			## -4332.84 				
			names(probXGivenInv) <<- NULL				
		}					
							
		getInvSitesParam(isInvProbSame, isReversible)					
	}						
							
}							
							
#############################################################							
##	Function called by "optim" for obtaining the new 						
##	estimates of the parameters						
##							
##	Input						
##		paramList: Current estimate of rate matrix					
##		paramType: Parameters to be updated					
##			   1 -> All the edge lengths				
##			   2 -> Pi-vector along an edge				
##			   3 -> S-matrix along an edge				
##		rateMatrixRows: Indexes of edges that have 					
##				the same rate matrix			
##	Output						
##		Revised estimates of the relevant parameters					
#############################################################							
							
getUpdatedParam = function(paramList, paramType, rateMatrixRows, setRootVector							
				, invSitesEqualVarSites, isRevMatrix) {			
							
	numEdges = dim(paramMatrix)[1]						
							
	if(isRevMatrix) edgeLengthCol = s+r+1			#v			
	else edgeLengthCol = r*(r-1)+1				#v		
							
	if(paramType == 1) {						
		paramMatrix[, edgeLengthCol] <<- paramList					
		for(i in 1:numEdges) globalProbMatrix[i, , ] <<- getCondProb(i, isRevMatrix)					
							
		return(getLogLikelihood(isRevMatrix))					
							
	} ## Edge lengths have been updated						
							
	if(paramType == 2) {						
		if(isRevMatrix) numParam = s			#v
		else numParam = colsirrev			#v
							
		for(rIndex in 1:length(rateMatrixRows)) {					
			paramMatrix[rateMatrixRows[rIndex], 1:numParam] <<- paramList				
			globalProbMatrix[rateMatrixRows[rIndex], , ] <<- getCondProb(rateMatrixRows[rIndex], isRevMatrix) 	# v				
		}					
							
		return(getLogLikelihood(isRevMatrix))					
							
	} ## S-matrix along an edge has been updated						
							
	if(paramType == 3) {						
		## Check that sum(pi_i) = 1 and each 0 < pi_i < 1					
		revPiVector = paramList/sum(paramList)					
							
		for(rIndex in 1:length(rateMatrixRows)) {					
			paramMatrix[rateMatrixRows[rIndex], (s+1):(s+r)] <<- revPiVector				
			globalProbMatrix[rateMatrixRows[rIndex], , ] <<- getCondProb(rateMatrixRows[rIndex], isRevMatrix)	# v				
		}					
							
		if(setRootVector) rootVector <<- revPiVector					
		if(invSitesEqualVarSites) probXGivenInv <<- revPiVector					
							
		return(getLogLikelihood(isRevMatrix))					
							
	} ## Pi-vector along an edge has been updated						
							
	return(0)						
}							
							
#############################################################################							
##	Obtain the Lagrange multiplier for updation of pi-vector at the 						
##	root node						
##							
##	Input						
##		Local					
##			site: Site number to obtain the relevant pattern 				
##				from the sequenceMatrix			
##		Global					
##			mergeMatrix: Rooted tree topology in matrix format				
##			rootVector: Pi-vector for the root				
##							
##	Output						
##		pi_i deriv[likelihood(site)]/ likelihood(site) 					
##
##
##		for pi_i = ith element of the State Space					
#############################################################################							
							
getLagrangeMultiplier = function(site) {							
							
	## The R function hclust returns a merge matrix						
	## that has the top-most bifurcation of phylogenetic						
	## tree as the last row						
	numRows = dim(mergeMatrix)[1]						
							
	index1 = mergeMatrix[numRows, 1]						
	index2 = mergeMatrix[numRows, 2]						
							
	prob1 = globalProbMatrix[match(index1, mergeVector), , ]						
	prob2 = globalProbMatrix[match(index2, mergeVector), , ]						
							
	siteLikelihood = 0						
	siteLikelihood = diag(rootVector) %*% matrix((getSubTreeLikelihood(prob1, index1, site) 						
						* getSubTreeLikelihood(prob2, index2, site)), ncol=1)	
							
	return(siteLikelihood)						
}							
							
#########################################################################							
##	Obtain the updated value of pi-vector for the root node						
##							
##	Input						
##		Global					
##			sequenceMatrix: Matrix of aligned sequences				
##					with only the unique patterns		
##			siteCount: Number of occurences of patterns 				
##				   in sequenceMatrix			
##							
##	Output						
##		Updated value of the pi-vector at the root node					
#########################################################################							
							
getUpdatedRootVector = function(isReversible) {							
							
	numSites = dim(sequenceMatrix)[2]						
	updateMatrix = rep(0, r)							#v						
							
	## Obtain conditional probability matrices						
	numEdges = dim(paramMatrix)[1]						
	for(i in 1:numEdges) globalProbMatrix[i, , ] <<- getCondProb(i, isReversible)	#v
							
	for(siteIndex in 1:numSites) {							
											
		temp = getLagrangeMultiplier(siteIndex)					
							
		index = match(siteIndex, cstSitesIndex)					
		if(is.na(index)) {					
			temp = temp/ sum(temp)				
			updateMatrix = updateMatrix + siteCount[siteIndex] * temp				
		}					
		else{					
			numr = alpha * temp				
			denom = alpha * sum(temp) + betaVar * probXGivenInv[index]				
			updateMatrix = updateMatrix + siteCount[siteIndex] * (numr/denom)				
		}					
	}						
	
	newVector=rep(0, r)								#v
	for(h in 1:r){									#v
		newVector[h]=updateMatrix[h]/sum(updateMatrix)				#v
	}										#v
	rootVector <<- newVector							#v
}							
							
#########################################################################							
##	Obtain the estimates for the invatriant sites parameters						
##							
##	Input						
##		Global					
##			sequenceMatrix: Matrix of aligned sequences				
##					with only the unique patterns		
##			siteCount: Number of occurences of patterns 				
##				   in sequenceMatrix			
##			cstSitesIndex: Index of constant sites in the 				
##						order from site 1 until site r. The order is		# v
##						assigned according to the original			# v
##						alphabetical order					# v
##			alpha: P(var_site)				
##			betaVar: P(inv_site)				
##			probXGivenInv: Vector of P(X|inv) where X=Set of elements of the state space	#v
##							
##	Output						
##		Updated values of the invariant sites parameters					
#########################################################################							
							
getInvSitesParam = function(invSitesEqualVarSites, isReversible) {							
							
	## Obtain conditional probability matrices						
	numEdges = dim(paramMatrix)[1]						
	for(i in 1:numEdges) globalProbMatrix[i, , ] <<- getCondProb(i, isReversible)		# v				
							
	numSites = dim(sequenceMatrix)[2]						
							
	## Obtain P(X|variable), where X = Set of elements of the state space 			# v						
	probXGivenVar = rep(0, r)								# v
	for(h in 1:r){ 										# v						
		probXGivenVar[h] = getSiteLikelihood(cstSitesIndex[h])				# v			
	}	
							
	alphaDeriv = 0						
	betaDeriv = 0						
							
	probXGivenInvDeriv = rep(0, r)								# v
							
	for(siteIndex in 1:numSites) {						
							
		index = match(siteIndex, cstSitesIndex)					
							
		if(is.na(index)) alphaDeriv = alphaDeriv + siteCount[siteIndex]/ alpha					
		else {					
			denom = alpha * probXGivenVar[index] + betaVar * probXGivenInv[index]				
			alphaDeriv = alphaDeriv + siteCount[siteIndex] * probXGivenVar[index]/ denom				
			betaDeriv = betaDeriv + siteCount[siteIndex] * probXGivenInv[index]/ denom				
			probXGivenInvDeriv[index] = probXGivenInvDeriv[index] + siteCount[siteIndex] * betaVar * probXGivenInv[index]/ denom				
		}					
							
	} ## All sites have been considered						
							
	## Assign new values to the global parameters						
	totalSites = sum(siteCount)						
	alpha <<- alpha * alphaDeriv / totalSites						
	betaVar <<- betaVar * betaDeriv / totalSites						
							
	if(!invSitesEqualVarSites) probXGivenInv <<- probXGivenInvDeriv/ sum(probXGivenInvDeriv)						
							
}							
							
##################################################################							
##	Convert a Newick format file to a "merge matrix"						
##	format file. This function is ONLY meant for reading 						
##	rooted trees						
##							
##	Input						
##		phylipFile: Newick format file					
##							
##	Output						
##		A list containing the following components					
##		1. Matrix in the same format as "merge matrix"					
##		2. Leaf nodes in the order in which they were					
##		   converted into -ve numbers					
##################################################################							
							
newickFormatToMergeMatrix = function(phylipFile) {							
							
	fileInfo = getModNewickFile(phylipFile)						
	newickFile = fileInfo$newick						
	fileLen = length(newickFile)						
							
	## Number of edges for a rooted tree = 2n-2, where n = number of 						
	## leaf nodes						
	distMatrix = matrix(0, ncol=2, nrow=(fileInfo$numLeafNode - 1))						
							
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
							
	return(list(dMat = distMatrix						
			, leaves = fileInfo$leaf))				
							
}							
							
##################################################################							
##	Convert the Newick format file into a vector such that						
##	(, ")", "," and leaf nodes can be read separately						
##							
##	Input						
##		inFile: Newick format file					
##							
##	Output						
##		A list containing the following components					
##		1. Newick tree with leaf nodes replaced by -ve					
##		   numbers and each symbol being a separate 					
##		   character					
##		2. Number of leaf nodes					
##		3. Names of leaf node in the order in which they					
##		   were converted into -ve numbers					
##################################################################							
							
getModNewickFile = function(inFile) {							
							
	## Read the Newick tre topology from inFile						
	pFile = scan(inFile, what="character")						
							
	## Remove whitespaces						
	modFile = gsub(" ", "", pFile)						
							
	tempFile = gsub("\\(", "", modFile)						
	tempFile = gsub("\\)", "", tempFile)						
	leafNodes = unlist(strsplit(tempFile, ","))						
							
	tempFile = modFile						
	for(i in 1:length(leafNodes)) {						
		node = -1 * i					
		tempFile = sub(leafNodes[i], node, tempFile)					
	}						
							
	tempFile = gsub("\\(", "\\(:", tempFile)						
	tempFile = gsub("\\)", ":\\)", tempFile)						
	tempFile = gsub(",", ":,:", tempFile)						
							
	treeComponents = unlist(strsplit(tempFile, split=":"))						
							
	return(list(newick = treeComponents						
			, numLeafNode = length(leafNodes)				
			, leaf = leafNodes))				
}							
							
##################################################################							
##	Change taxa-order in input matrix						
##							
##	Input						
##		seqMat: Input sequence matrix					
##		distMatOrder: Taxa names in the order in which they					
##				  occur in the mergeMatrix			
##							
##	Output						
##		A matrix of sequences with the rows corresponding					
##		to taxa sequences in the correct order					
##################################################################							
							
changeTaxaOrder = function(seqMat, distMatOrder) {							
							
	## Remove blank spaces from taxa names						
	taxaNames = gsub("[ ]*", "", rownames(seqMat))						
							
	indexes = match(distMatOrder, taxaNames)						
							
	revSeqMat = seqMat[indexes, ]						
	rownames(revSeqMat) = distMatOrder						
							
	return(revSeqMat)						
}							
							
##################################################################							
##	Convert the merge matrix to a Newick tree						
##							
##	Input						
##		hMat: Merge Matrix					
##		leaves: Vector of species names					
##							
##	Output						
##		A rooted Newick tree					
##################################################################							
							
mergeMatrixToNewickFormat = function(hMat, leaves) {							
							
	numRows = dim(hMat)[1]						
							
	newickTree = paste("("						
				, paddedVal(hMat[numRows, 1])			
				, ","			
				, paddedVal(hMat[numRows, 2])			
				, ")"			
				, sep="")			
							
	for(i in (numRows-1):1) {						
							
		addRow = paste("("					
				, paddedVal(hMat[i, 1])			
				, ","			
				, paddedVal(hMat[i, 2])			
				, ")"			
				, sep="")			
							
		rowPattern = paste("@", i, "@", sep="")					
		newickTree = sub(rowPattern, addRow, newickTree)					
							
	} ## The mergeMatrix has been considered						
							
	## Replace leaf indexes with leaf names						
	for(i in 1:length(leaves)) {						
							
		index = paste("#", -i, "#", sep="")					
		newickTree = sub(index, leaves[i], newickTree)					
							
	}						
							
	return(newickTree)						
}							
							
##################################################################							
##	Distinguish between leaf nodes and internal nodes						
##							
##	Input						
##		item: Node index					
##							
##	Output						
##		Node index appropriately padded to distinguish between					
##		leaf nodes and internal nodes 					
##################################################################							
							
paddedVal = function(item) {							
							
	if(item < 0) newItem = paste("#", item, "#", sep="")						
	else newItem = paste("@", item, "@", sep="")						
							
	return(newItem)						
}							
							
							
###########################################################################							
##	Initialize the global model parameters						
##							
##	Input						
##		inOutDir: Directory of source files					
##		sequenceFileName: Seqeunce file in PHYLIP format					
##		paramFileName: Input parameters file for S-matrix, pi-vector					
##				   and edge length			
##		newickTreeTopFile: Newick tree topology for a rooted tree					
##		initRootVect: Root vector			
##							
##	Output						
##		The following global variables are initialized					
##		1. mergeMatrix: Matrix of tree topology					
##		2. mergeVector: mergeMatrix in vector form					
##		3. sequenceMatrix: Matrix of unique patterns					
##		4. siteCount: Number of sites of each pattern					
##		5. cstSitesIndex: Index of constant sites					
##		6. alpha, betaVar: P(site is variable), P(site is invariable)					
##		7. probXGivenInv: P(X|inv), where X={A,C,G,T}					
##		8. rootVector: Root vector					
##		9. paramMatrix: K x 11 matrix with each row corresponding					
##		   to 6 S-matrix elements + 4 Pi-vector elements +edge length					
###########################################################################							
							
initializeGlobalParam = function(inOutDir, sequenceFileName, paramFileName			#v I removed the option seqLength for now. It lied between sequenceFileName and paramFileName							
				, newickTreeTopFile, initRootVect, isReversible
				, Group) {			
							
	paramFile = paste(inOutDir, paramFileName, sep="/")						
	sequenceFile = paste(inOutDir, sequenceFileName, sep="/")						

	## Define the r variable		v
	     r 	  <<- length(initRootVect) 	# v
	     s 	  <<- r*(r-1)/2			# v	
	colsirrev <<- r*(r-1)			# v This variable will only be required for the non-reversible-rate-matrix cases. 	

	## Defining the group to reduce the parameters

	 
						
	## Read the aligned sequences						
	inputSeqMat = getSequences(sequenceFile)			 #v						

							
	## Obtain the merge matrix. 						
	treeTopFile = paste(inOutDir, newickTreeTopFile, sep="/")						
	hc = newickFormatToMergeMatrix(treeTopFile)						
	mergeMatrix <<- hc$dMat						
							
	## Convert the merge matrix into a vector						
	mergeVector <<- NULL						
	numRows = dim(mergeMatrix)[1]						
	for(i in 1:numRows) mergeVector <<- c(mergeVector, as.vector(mergeMatrix[i, ]))						
							
	## Change order of taxa in inputSeqMat						
	inputSeqMat = changeTaxaOrder(inputSeqMat, hc$leaves)						
							
	## Obtain list of <pattern, pattern_count>						
	x = getUniquePatterns(inputSeqMat)						
	sequenceMatrix <<- x$mat						
	siteCount <<- x$patCount						
	cstSitesIndex <<- x$cstIndex						
							
	alpha <<- 1						
	betaVar <<- 1 - alpha						
	probXGivenInv <<- rep(1/r, r)	#v					
							
	## Set the root vector						
	rootVector <<- initRootVect						
							
	paramMatrix <<- read.table(paramFile, header=FALSE)		#v I replaced the original read.csv(paramFile, header=FALSE, sep="\t") to avoid problems between tabultaions and blank spaces 						
	paramMatrix <<- as.matrix(paramMatrix)						
							
	if(isReversible && dim(paramMatrix)[2] != s+r+1) stop("Input parameter matrix should have s+r+1 values per row")   	     # v						
	if(!isReversible && dim(paramMatrix)[2] != colsirrev+1) stop("Input parameter matrix should have r(r-1)+1 values per row")   # v						
							
}							
							
###########################################################################							
##	Initialize the global model parameters and call the mleParam 						
##	function.						
##							
##	Input						
##		inOutDir: Directory of output file					
##		outFile: Output file name					
##		grpSmat: Grouping of S-matrices					
##		grpPiVector: Grouping of pi-vectors with the first element 					
##				being	the root vector		
##		numIterations: Number of iterations					
##		invProbSame: Is P(X|inv) = P(X|var)					
##		computeInv: Estimate P(site is invariable)					
##		beta: Probability of a site being invariant					
##		pxInv: P(X|inv)					
##		verbose: Display stepwise iteration information
##							
##	Output						
##		The following global variables are initialized					
##		1. sMat: Grouping of S-matrices					
##		2. piVect: Grouping of pi-vectors					
##		3. globalProbMatrix: The condition prob matrices per edge					
###########################################################################							
							
estimateMleParam = function(inOutDir, outFile, grpSmat, grpPiVector							
					, numIterations, invProbSame, computeInv		
					, beta = 0, pxInv, verbose, isReversible) {	#v
							
	startTime = Sys.time()						
							
	## Set the value of beta to 0.75 * proportion_of_inv_sites						
	if(computeInv) {						
		if(beta == 0) {					
			betaVar <<- 0.75 * sum(siteCount[cstSitesIndex])/ sum(siteCount)				
			alpha <<- 1 - betaVar						
			probXGivenInv <<- rep(1/r, r)						#v
		}					
		else {					
			betaVar <<- beta				
			alpha <<- 1 - betaVar				
			probXGivenInv <<- pxInv				
		}					
	}						
							
	sMat <<- grpSmat						
	piVect <<- grpPiVector						
							
	globalProbMatrix <<- array(0, c(dim(paramMatrix)[1], r, r))				#v
							
	logL = mleParam(displayInfo = verbose, computeInvSites = computeInv						
			, isInvProbSame = invProbSame, numItr = numIterations				
			, inOutFile = paste(inOutDir, outFile, sep="/")				
			, isReversible)							#v
							
	endTime = Sys.time()						
	print("Execution time")						
	print(endTime-startTime)						
							
	return(logL)						
}							
							
##############################################################################							
##	Generate parametric bootstrap replicates						
##							
##	Input						
##		numReplicates: Number of replicates					
##		numSites: Total number of sites					
##		paramFileName: Input file with parameters for generarting					
##				   bootstrap replicates			
##		outFileName: Output file for storing bootstrap replicates					
##							
##	Output						
##		Bootstrap replicates in a user-defined file in PHYLIP format					
##############################################################################							
##JR has made two changes to get correct output:
##on l1734 changed sep=" " to sep="" to keep 10 spaces for taxa names
##on 1741 added if(numReplicates>1) so that single bootstrap replicates could be obtained							
#################################################################################

generateBootstraps = function(numReplicates, numSites, paramFileName				
					, outFileName, isReversible, StateSpaceVector=NULL) {	
	
	if(exists("r")==FALSE){
		num.taxa    <- scan(paramFileName, skip=1, nlines=1)+1		# v
		skip.lines  <- 19 + 5*num.taxa					# v [23 constant rows] + [ 2(2Taxa-2) rows (Two matrices of parameters)] + [(taxa-1) rows (These are the rows of the tree topology)]+ [1 row (that indicates the number of rows in the tree topology)]	
		initialRoot <- scan(paramFileName, skip=skip.lines, nlines=1)	# v
		r	    <- length(initialRoot)				# v
		s	    <- r*(r-1)/2					# v
		colsirrev   <- r*(r-1)						# v
		NumericStateSpace <- 1:r
	}									# V
										# v
	r	  <<- r								# v
	s	  <<- s								# v
	colsirrev <<- colsirrev							# v
	NumericStateSpace <<- NumericStateSpace					# v

	if(exists("StateSpace")==TRUE) {StateSpace=StateSpace} 									#v
	else {															#v
		if(is.null(StateSpaceVector)==TRUE) {										#v
			if(r==4) {StateSpace = c("A","C","G","T")} else {stop("you need to provide a StateSpaceVector")}	#v
		} else {StateSpace=sort(StateSpaceVector)}									#v
	}

	StateSpace <<- StateSpace
	lineCount = 1 	## Ignore the first line of paramFileName
							
	## Set the global parameters						
	mergeMatrixRows = scan(paramFileName, skip=lineCount, nlines=1)		
	mergeMatrix <<- matrix(0, ncol=2, nrow=mergeMatrixRows)			
							
	lineCount = lineCount + 1						
	for(i in 1:mergeMatrixRows) {						
		temp = scan(paramFileName, skip=lineCount, nlines=1, what="character")		
		temp = as.numeric(gsub("\\+0i", "", temp))					
		mergeMatrix[i, ] <<- temp					
		lineCount = lineCount + 1					
	}						
							
	## Convert the merge matrix into a vector						
	mergeVector <<- NULL						
	numRows = dim(mergeMatrix)[1]						
	for(i in 1:numRows) mergeVector <<- c(mergeVector, as.vector(mergeMatrix[i, ]))		
							
	leafNodeLabels = scan(paramFileName, what="character", skip=lineCount, nlines=1)	
							
	## Ignore the input parameters						
	lineCount = lineCount + 8 + numRows * 2 + 13 						
							
	if(isReversible) {numParam = s+r+1} else {numParam = colsirrev+1}		#v Gotta fix this
							
	paramMatrix <<- matrix(0, ncol=numParam, nrow=2*mergeMatrixRows)		
	for(i in 1:(2*mergeMatrixRows)) {						
		paramMatrix[i, ] <<- scan(paramFileName, skip=lineCount, nlines=1)	
		lineCount = lineCount + 1					
	}						
							
	lineCount = lineCount + 1						
	rootVector <<- scan(paramFileName, skip=lineCount, nlines=1)		
							
	lineCount = lineCount + 2						
	betaVar <<- scan(paramFileName, skip=lineCount, nlines=1)		
	alpha <<- 1 - betaVar						
							
	lineCount = lineCount + 2						
	probXGivenInv <<- scan(paramFileName, skip=lineCount, nlines=1)		
							
	## Set the conditional probability matrices				
	globalProbMatrix <<- array(0, c(dim(paramMatrix)[1], r, r))					#v
							
	for(i in 1:dim(paramMatrix)[1]) globalProbMatrix[i, , ] <<- getCondProb(i, isReversible)	#v
							
	## Delete outFile if it already exists						
	unlink(outFileName)			

	if(r<=9){
		for(replicates in 1:numReplicates) {								#v If the state space <=9, I keep the same code a Vivek's was. However, 					
														#v I've had to do some modifications for states spaces >= 10. (See below)
			bsMat = getReplicate(numSites, leafNodeLabels)						#v LeafNodeLabels are the corresponding taxa names
			write(paste(dim(bsMat)[1], dim(bsMat)[2], sep="\t"), outFileName, append=TRUE)		#v This is only saying the dimension of each replicate matrix
			for(j in 1:length(leafNodeLabels)) {					
								
				## Obtain 10 character leaf nodes 						#v He makes this to keep the taxa name
				nodeName = paste(leafNodeLabels[j], "          ", sep="")			#v with ten characters at most!
				nodeName = substr(nodeName, 1, 10)				
									
				## Obtain the concatenated sequence				
				temp = bsMat[j, ]								#v
				for(h in 1:r){									#v
					temp=gsub(NumericStateSpace[h], StateSpace[h], temp)			#v
				}										#v
				
				sequence = temp[1]				
				for(k in 2:numSites) sequence = paste(sequence, temp[k], sep="")		
				
				write(paste(nodeName, sequence, sep=""), outFileName, append=TRUE)		#v I have put at least one blanket space between the taxa name and the sequence
				
			} ## All sequences for a replicate have been saved					
						
			## Blank space between two replicates					
			if(numReplicates>1)write(" ", outFileName, append=TRUE)					
								
		} ## All replicates have been generated						
	} else {
		NumStatSpacVect 	<- paste(NumericStateSpace)							#v These steps are made in order to keep the sequence in the state space as: 
		NumStatSpacVect[1:9]	<- paste(0, NumStatSpacVect[1:9], sep="")					#v ("01", "02", ..., "09", "10", "11", ... )
															
		for(replicates in 1:numReplicates) {						
			bsMat = getReplicate(numSites, leafNodeLabels)							#v This line generates a bootstrap-aligned-sequence matrix. Each number corresponds to an element in the state space. 
			write(paste(dim(bsMat)[1], dim(bsMat)[2], sep="\t"), outFileName, append=TRUE)			#v This is only saying the dimension of each replicate matrix.
															#v Note that length(leafNodeLabels) = dim(bsMat)[1] and numSites = dim(bsMat)[2] 
			bsMat = paste(bsMat)										#v
			for(j in 1:length(bsMat)){									#v
				if(as.numeric(bsMat[j])<=9) bsMat[j]=paste(0, bsMat[j], sep ="") else bsMat[j]=bsMat[j]	#v
			}												#v

			for(j in 1:r) {bsMat = gsub(NumStatSpacVect[j], StateSpace[j], bsMat)}				#v Recoding from numbers to letters.
			bsMat = matrix(bsMat, length(leafNodeLabels), numSites)						#v bsMat is returned to its matrix form. Yes, everything has been checked so 
			
			for(j in 1:length(leafNodeLabels)) {					
								
				## Obtain 10 character leaf nodes 							#v I am keeping the same first ten characters to the name of the taxa	
				nodeName = paste(leafNodeLabels[j], "          ", sep="")				
				nodeName = substr(nodeName, 1, 10)				
									
				## Obtaining the concatenated sequence												 	
				temp = bsMat[j, ]									#v					
				sequence = temp[1]									#v I removed the gsub part from here as it's been done above
				for(k in 2:numSites) sequence = paste(sequence, temp[k], sep="")			#v	
				
				write(paste(nodeName, sequence, sep=" "), outFileName, append=TRUE)
				
			} ## All sequences for a replicate have been saved					
						
			## Blank space between two replicates					
			write(" ", outFileName, append=TRUE)					
								
		} ## All replicates have been generated						

	}							
}							
							
##############################################################################							
##	Obtain a single bootstrap replicate						
##							
##	Input						
##		Local					
##			numSites: Total number of sites				
##			leafNames: Vector of leaf names such that -1 corresponds 				
##					to the first leaf and so on		
##		Global					
##			globalProbMatrix: Conditional probability matrix per edge				
##						in the same order as the mergeMatrix	
##			mergeMatrtix: Merge matrix (rooted tree topology)				
##			alpha: Proportion of variable sites				
##			rootVector: Nucleotide frequencies at the root				
##							
##	Output						
##		Matrix of sequences such that each column corresponds 					
##		to a particular site.					
##############################################################################							
##v  The convertion to characters does NOT happen here
###################################################################################

							
getReplicate = function(numSites, leafNames){							
							
	NV = rbinom(1, numSites, alpha)						
	numRows = dim(mergeMatrix)[1]						
							
	seqInt = matrix(0,numRows,NV)						
	seqInt[numRows,] = sample(1:r, NV, prob=rootVector, replace=TRUE)								#v
							
	seqLeaf = matrix(0, numRows+1, NV)						
							
	for (j in 1:numRows){						
							
		index1 = mergeMatrix[numRows+1-j, 1]					
		index2 = mergeMatrix[numRows+1-j, 2]					
							
		prob1 = globalProbMatrix[match(index1, mergeVector), , ]					
		prob2 = globalProbMatrix[match(index2, mergeVector), , ]					
							
		if (index1<0){					
        		for (i in 1:NV) seqLeaf[abs(index1),i] = sample(1:r, 1, replace=FALSE, prob1[seqInt[numRows+1-j,i],])		#v			
		} 					
		else {     					
		       for (i in 1:NV) seqInt[index1,i] = sample(1:r, 1, replace=FALSE, prob1[seqInt[numRows+1-j,i],])			#v		
		}					
        							
		if(index2<0){					
			for (i in 1:NV) seqLeaf[abs(index2),i] = sample(1:r, 1, replace=FALSE, prob2[seqInt[numRows+1-j,i],])		#v		
          	}						
	       else {						
		       for (i in 1:NV) seqInt[index2,i] = sample(1:r, 1, replace=FALSE, prob2[seqInt[numRows+1-j,i],])			#v		
		}					
							
	}						
							
	invSites = rmultinom(1, numSites-NV, probXGivenInv)						

	for(h in 1:r){															#v
		state=matrix(NumericStateSpace[h], numRows+1, invSites[h])								#v
		seqLeaf=cbind(seqLeaf, state)												#v
	}																#v
							
	rownames(seqLeaf) = leafNames						
	return(seqLeaf)						
							
}							
							
##############################################################################							
##	Obtain the PHYLIP file from a file with multiple bootstrap replicates						
##							
##	Input						
##		replicateFile: File with multiple bootstrap replicates					
##		phylipFile: File for a single bootstrap replicate					
##		replicate: Replicate number					
##							
##	Output						
##		A PHYLIP format file					
##############################################################################							
							
getPhylipFile = function(replicateFile, phylipFile, replicate) {							
							
	numSpecies = scan(replicateFile, nlines=1)[1]						
	linesToSkip = (replicate - 1) * (numSpecies + 2) - 1		## First line is <numTaxa, numSequence>, last line is blank				
							
	for(i in 1:(numSpecies+1)) {						
							
		info = scan(replicateFile, what="character", skip=linesToSkip+i, nlines=1, sep=":")					
		write(info, phylipFile, append=TRUE)					
							
	} ## PHYLIP file generated for a given replicate						
							
}							
							
##############################################################################							
##	Perform maximum-likelihood analysis of multiple datasets and save 						
##	the output per dataset in a separate out file						
##							
##	Input						
##		inOutDir: Directory of source files					
##		seqLength: Number of sites per sequence					
##		paramFileName: Input parameters file for S-matrix, pi-vector					
##				   and edge length			
##		newickTreeTopFile: Newick tree topology for a rooted tree					
##		numIterations: Number of iterations					
##		invProbSame: Is P(X|inv) = P(X|var)					
##		computeInv: Estimate P(site is invariable)					
##		verbose: Display stepwise iteration information					
##		numRep: Number of datasets to analyze					
##		replicateFile: File that contains the datasets					
##							
##							
##	Output						
##		Output files containing maximum-likelihood parameters per dataset					
##############################################################################							
							
analyzeData = function(inOutDir, paramFileName, newickTreeTopFile							# I removed the command seqLength. I appeared between inOutDir and paramFileName
				, numIterations, invProbSame, computeInv			
				, verbose, numRep, replicateFile			
				, isReversible = TRUE) {			
							
	phylipFile = paste(inOutDir, "bsRep.txt", sep="/")						
	repFile = paste(inOutDir, replicateFile, sep="/")						
							
	absParamFile = paste(inOutDir, paramFileName, sep="/")						
							
	## Read paramFile to obtain the following vectors						
	## 1. grpSmat 2. grpPiVector 3. inBeta 4. inPGivenInv						
	## 5. S-matrix + Pi-vector + time per edge						
	lineCount = 1						
	numRows = scan(absParamFile, skip=lineCount, nlines=1)						
							
	lineCount = lineCount + 1 + numRows + 4						
	grpSmat = scan(absParamFile, skip=lineCount, nlines=1)						
							
	lineCount = lineCount + 2						
	grpPiVector = scan(absParamFile, skip=lineCount, nlines=1)						
							
	lineCount = lineCount + 2 + numRows*2 + 13						
	tempParamFile = paste(inOutDir, "tempParamFile.txt", sep="/")						
							
	## Delete an existing tempParamFile						
	unlink(tempParamFile)						
							
	for(edge in 1:(numRows*2)) {						
		tempVector = scan(absParamFile, skip=lineCount, nlines=1, what="character")					
		tempVector = as.numeric(gsub("\\+0i", "", tempVector))					
							
		write(tempVector, tempParamFile, ncolumns=length(tempVector), append=TRUE, sep="\t")					
		lineCount = lineCount + 1					
	} ## All edges have been considered						
							
	lineCount = lineCount + 1						
	initRootVect = scan(absParamFile, skip=lineCount, nlines=1)						
							
	lineCount = lineCount + 2						
	inBeta = scan(absParamFile, skip=lineCount, nlines=1)						
							
	lineCount = lineCount + 2						
	inPGivenInv = scan(absParamFile, skip=lineCount, nlines=1)						
							
	for(replicates in 1:numRep) {						
							
		## Obtain the phylip file for a particular dataset					
		unlink(phylipFile)					
		getPhylipFile(repFile, phylipFile, replicates)					
							
		## Initialize the global parameters					
		initializeGlobalParam(inOutDir, "bsRep.txt", "tempParamFile.txt"			#v Since I removed the option seqLength for the original "initializeGlobalParam" function. I have to delete it here as well. Originally, It lied between sequenceFileName and paramFileName		
						, newickTreeTopFile, initRootVect, isReversible)	
							
		## Name of the unique output file per dataset					
		outFile = paste("Output", replicates, ".out", sep="")					
							
		## Obtain the maximum-likelihood estimates per dataset					
		estimateMleParam(inOutDir, outFile, grpSmat, grpPiVector					
					, numIterations, invProbSame, computeInv		
					, beta = inBeta		
					, pxInv = inPGivenInv		
					, verbose, isReversible)					#v
							
	} ## All replicates have been considered						
							
}							
							
##############################################################################							
##	Obtain the mean and standard deviation of the bootstrap replicates 						
##	for each parameter						
##							
##	Input						
##		fileDir: Directory where .out files are saved, one per replicate					
##		numTaxa: Number of taxa					
##		meanStdFile: Output file for mean and std deviation per parameter					
##		deltaLogFile: Output file containing deltaLogL					
##							
##	Output						
##		Two output files, one containing the mean and std deviation					
##		and the other containing deltaLogL value per replicate					
##############################################################################							
							
computeStatistics = function(fileDir, numTaxa, meanStdFile, deltaLogFile, isReversible) {							
							
	absMeanFile = paste(fileDir, meanStdFile, sep="/")						
	absDeltaFile = paste(fileDir, deltaLogFile, sep="/")						
							
	outNames = dir(fileDir, pattern = ".out")						
	nFiles = length(outNames)						
							
	numRows = numTaxa - 1						
							
	## Initialize parameters						
	deltaLogL = rep(0, nFiles)						
	allBeta = rep(0, nFiles)						
							
	allRootVector = matrix(0, nrow=nFiles, ncol=r)			#v						
	allProbXGivenInv = matrix(0, nrow=nFiles, ncol=r)		#v				
							
	## For each edge, create a matrix of s+r+1 columns		#v
	if(isReversible) edgeLenCount = s+r+1				#v
	else edgeLenCount = colsirrev+1					#v
							
	arrayValues = array(0, c(2*numRows, nFiles, edgeLenCount))						
							
	for(indexes in 1:nFiles) {						
							
		fileName = paste(fileDir, outNames[indexes], sep="/")					
		skipLines = 1 + numRows + 20 + 2*numRows					
							
		## Read log-likelihood					
		logLikelihoods = scan(file = fileName, what = double(0), skip = skipLines, nlines = 1, sep="\t")					
		deltaLogL[indexes] = logLikelihoods[2] - logLikelihoods[1]					
							
		## Read the edge parameters					
		skipLines = skipLines + 2					
							
		for(temp in 1:(2*numRows)) {					
							
			tempData = scan(file = fileName, skip = skipLines, nlines = 1				
					, sep = "\t", what="character")		
							
			tempData = as.numeric(gsub("\\+0i", "", tempData))				
			arrayValues[temp, indexes, ] = tempData				
							
			skipLines = skipLines + 1				
							
		} ## All edges have been read					
							
		## Read root vector					
		skipLines = skipLines + 1					
		allRootVector[indexes, ] = scan(file = fileName, what = double(0), skip = skipLines, nlines = 1)					
							
		## Read beta					
		skipLines = skipLines + 2					
		allBeta[indexes] = scan(file = fileName, what = double(0), skip = skipLines, nlines = 1)					
							
		## Read P(X|inv)					
		skipLines = skipLines + 2					
		allProbXGivenInv[indexes, ] = scan(file = fileName, what = double(0), skip = skipLines, nlines = 1)					
							
	} ## All output files have been read						
							
	## Save all deltaLogLikelihoods in a file						
	write(deltaLogL, absDeltaFile, ncolumns=1)						
							
	## Save avg(beta) and std_dev(beta)						
	saveScalar(allBeta, "Beta", absMeanFile)						
							
	## Save avg(rootVector) and std_dev(rootVector)						
	saveVector(allRootVector, "Root Vector", absMeanFile)						
							
	## Save avg[P(X|inv)] and std_dev[P(X|inv)]						
	saveVector(allProbXGivenInv, "P(X|inv)", absMeanFile)						
							
	## Save edge parameters						
	for(i in 1:(2*numRows)) {						
		saveVector(as.matrix(arrayValues[i, ,])					
				, paste("Edge", i, sep = ":")			
				, absMeanFile)			
	} ## All edges have been considered						
							
}							
							
##############################################################################							
##	Obtain the mean and standard deviation of the bootstrap replicates 						
##	for edge parameters, root vector, and P(X|inv)						
##							
##	Input						
##		inputMat: Matrix of parameters					
##		info: Matrix description					
##		fileName: Output file					
##							
##	Output						
##		Output file containing the mean and std deviation for the 					
##		appropriate parameters					
##############################################################################							
							
saveVector = function(inputMat, info, fileName) {							
							
	meanVal = apply(inputMat, 2, mean)						
	stdVal = sqrt(apply(inputMat, 2, var))						
							
	temp = NULL						
							
	for(i in 1:length(meanVal)) temp = c(temp, paste(meanVal[i], stdVal[i], sep="+/-"))						
	write(info, fileName, append=TRUE)						
	write(temp, fileName, append=TRUE, sep="\t", ncolumns=length(temp))						
							
}							
							
##############################################################################							
##	Obtain the mean and standard deviation of the bootstrap replicates 						
##	for scalar parameters - log-likelihood and beta						
##							
##	Input						
##		inputMat: Matrix of parameters					
##		info: Matrix description					
##		fileName: Output file					
##							
##	Output						
##		Output file containing the mean and std deviation for the 					
##		appropriate parameters					
##############################################################################							
							
saveScalar = function(inputVector, info, fileName) {							
							
	meanVal = mean(inputVector)						
	stdVal = sqrt(var(inputVector))						
	write(info, fileName, append=TRUE)						
	write(paste(meanVal, stdVal, sep="+/-"), fileName, append=TRUE)						
}							
							
##############################################################################							
##	Calculation of standardized rate matrix and the edge lengths for a 						
##	time-reversible process						
##							
##	Input						
##		inFile: File containing the input and output parameters					
##							
##	Output						
##v		A rate matrix with s+r+1 elements per row. The first s elements 					
##v		correspond to the elements s1, s2, ... with f=...?. The next r elements					
##v		correpsond to the pi-vector and the last element is the edge					
##v		length					
##############################################################################							
							
getStandardRateMatrix = function(inFile) {							
							
	## Determine the number of lines to exclude						
	numTaxa = scan(inFile, skip=1, nlines=1)						
	linesToExclude = 2 + numTaxa + 8 + 2*numTaxa + 13								#v Don't touch this line. It is right 					
							
	rateTimeMatrix = matrix(0, nrow=2*numTaxa, ncol=r+s+1)								#v
							
	for(lines in 1:(2*numTaxa)) {						
		data = scan(inFile, skip=linesToExclude, nlines=1)					
		linesToExclude = linesToExclude + 1					
							
		tempSmat = getSmatrix(data[1:s], data[(s+1):(s+r)])							#v					
		# tempSmat = getSmatrix(Gy.to.S(data[1:s], data[(s+1):(s+r)], Group), data[(s+1):(s+r)])			#vforlump
		tempRmat = tempSmat %*% diag(data[(s+1):(s+r)])								#v
		tempStandardize = -sum(data[(s+1):(s+r)] * diag(tempRmat))						#v
							
		## 1. Divide the S-matrix elements by the standardization term					
		## 2. Multiply time by the standardization term. However, this term has to be 					
		##    divided by two as the current estimate of time is from node to node while					
		##    the "accepted" definition of time splits the edge into two edges with a "root"					
		##    in the middle					
		rateTimeMatrix[lines, ] = c(data[1:s]/tempStandardize, data[(s+1):(s+r)], 0.5 * tempStandardize)	#v					
							
	} ## All edges have been considered						
							
	rateTimeMatrix[, 1:s] = rateTimeMatrix[, 1:s]/rateTimeMatrix[,s]						#v
							
	return(rateTimeMatrix)						
}							
							
######################################################################################							
##	Calculation of standardized rate matrix and the edge lengths. The edge						
##	lengths are twice of that returned by std programs like PhyML						
##							
##	Input						
##v		pMatrix: Matrix of (r+s+1) or (r(r-1)+1) elements per edge.
##v			 Notice that if r=3, then r+s+1 = r(r-1)+1.
##v			 Therefore, the input isReversible is necessary. 		
##v		isReversible: TRUE/FALSE
##							
##	Output						
##v		Rate and time matrix with (r+s+1)/(2s) elements per row. For a time-reversible					
##v		process, the first s elements correspond to the elements s1, s2, ... , the next 					
##v		r elements correpsond to the pi-vector and the last element is the edge					
##v		length. For an unreversible matrix, the first 2s=r(r-1)+1 elements represent					
##v		the rate matrix elements and the last element is the edge length					
######################################################################################							
							
getNormalizedValues = function(pMatrix, isReversible) {												#v							
							
	numElements = dim(pMatrix)						
	stdMatrix = matrix(0, nrow=numElements[1], ncol=numElements[2])						
	tempRmat = NULL						
	piVect = NULL						
	if(isReversible){															#v						
		for(i in 1:numElements[1]) {													#v
			tempSmat = getSmatrix(pMatrix[i, 1:s], pMatrix[i, (s+1):(s+r)])							#v		
			# tempSmat = getSmatrix(Gy.to.S(pMatrix[i, 1:s], pMatrix[i, (s+1):(s+r)], Group), pMatrix[i, (s+1):(s+r)]) 		#vforlump
			tempRmat = tempSmat %*% diag(pMatrix[i, (s+1):(s+r)])									#v	
			piVect = pMatrix[i, (s+1):(s+r)]											#v				
			tempStandardize = -sum(piVect * diag(tempRmat))										#v
			stdMatrix[i, ] = c(pMatrix[i, 1:s]/tempStandardize, pMatrix[i, (s+1):(s+r)], tempStandardize*pMatrix[i, s+r+1])		#v
		}																#v
	} else {																#v
		for(i in 1:numElements[1]){													#v
			tempRup=tempRlo=matrix(0, r, r)												#v
			tempRup=lower.tri(tempRup); tempRup[lower.tri(tempRup)]=pMatrix[i, 1:(colsirrev/2)]; tempRup=t(tempRup)			#v
			tempRlo=lower.tri(tempRlo); tempRlo[lower.tri(tempRlo)]=pMatrix[i, (colsirrev/2+1):colsirrev]				#v
			tempRmat=tempRup+tempRlo												#v
			tempRmat = tempRmat - diag(apply(tempRmat, 1, sum))									#v			
			## Obtain the equilibrium/stationary probability									#v
			piVect = eigen(t(tempRmat))$vector[,r]											#v 
			piVect = piVect/sum(piVect)												#v
			tempStandardize = -sum(piVect * diag(tempRmat))										#v
			stdMatrix[i, ] = c(pMatrix[i, 1:colsirrev]/tempStandardize, tempStandardize*pMatrix[i, colsirrev+1])			#v
		}																#v
	} ##v All rows have been considered in both cases.					
	return(stdMatrix)															#v		
}							
							
######################################################################################							
##	Reduce the number of rate matrices based on equilibrium frequencies or zero						
##	edge lengths						
##							
##	Input						
##		inFile: Input file containing the input and output parameters					
##			  for a given number of distinct rate matrices				
##		threshHeight: Cut-off value below which the top two (and not the topmost)					
##				  rows of the the height matrix are used to determine the 			
##				  the number of distinct rate matrices			
##		isRev: TRUE/FALSE					
##		distMethod: Method used for computing pairwise distances					
##							
##	Output						
##		A list containig two elements					
##		1. index1: Vector of rate matrices obtained using equilibrium frequencies					
##		2. index2: Vector of rarte matrices obtained using zero edge					
######################################################################################							
							
reduceModelComplexity = function(inFile, threshHeight = 0.025, isRev, distMethod, distCalc) {							
							
	numTaxa = scan(inFile, skip=1, nlines=1)						
							
	## Read the number of unique pi-vectors						
	linesToExclude = 2 + numTaxa + 6						
	allPis = scan(inFile, skip=linesToExclude, nlines=1)						
							
	## Exclude the root						
	root = allPis[1]						
	allPis = allPis[-1]						
							
	## Read the marginal probabilities						
	linesToExclude = 2 + numTaxa + 8 + 2*numTaxa + 13						
							
	if(isRev) edgeMatrix = matrix(0, nrow=2*numTaxa, ncol=s+r+1)		##v Assume a reversible model per edge					
	else edgeMatrix = matrix(0, nrow=2*numTaxa, ncol=colsirrev+1)		##v Assume a non-reversible model per edge				
							
	for(lines in 1:(2*numTaxa)) {						
							
		if(isRev) data = scan(inFile, skip=linesToExclude, nlines=1)					
		else {					
			data = scan(inFile, skip=linesToExclude, nlines=1, what="character")				
			data = as.numeric(gsub("\\+0i", "", data))				
		}					
							
		linesToExclude = linesToExclude + 1					
		edgeMatrix[lines, ] = data					
	} ## All edges have been considered						
							
	## Lower model complexity by considering the 						
	## the two closest edges in terms of pi-vectors						
							
	nonZeroEdge = combinePiVectors(edgeMatrix, allPis, 1, isRev, distMethod, distCalc)						
	grpVectorNonZeroEdge = nonZeroEdge$grps						
							
	if((length(nonZeroEdge$ht) > 1) && ((nonZeroEdge$ht[2] - nonZeroEdge$ht[1]) < threshHeight)) {						
		info = combinePiVectors(edgeMatrix, allPis, topRows=2, isRev, distMethod, distCalc)$grps					
		compareAndWrite(info)					
	}						
							
	## Lower model complexity by combining the parent and child edges						
	## If the edge-reduction necessarily requires the combining of 						
	## of an edge linked to the root, then ignore the grouping						
							
	grpVectorZeroEdge = getZeroEdge(edgeMatrix, allPis, isRev)						
	if(length(unique(grpVectorZeroEdge)) > length(unique(grpVectorNonZeroEdge))) grpVectorZeroEdge = grpVectorNonZeroEdge						
							
	return(list(index1 = grpVectorNonZeroEdge, index2 = grpVectorZeroEdge))						
							
}							
							
######################################################################################							
##	Reduce the number of rate matrices based on equilibrium frequencies 						
##							
##	Input						
##		eMatrix: Matrix containing 2n-2 rows and s+r+1 columns with the columns 					
##			   (s+1):(s+r) corresponding to pi-vectors				
##		allPisVect: Current assignment of rate matrices to different edges					
##				excluding the root			
##		topRows: Number of rows to consider for combining the rate matrices					
##		isRev: TRUE/FALSE					
##		distMethod: Method used for computing pairwise distances					
##							
##	Output						
##		A list containig two elements					
##		1. grps: Vector of rate matrices obtained using equilibrium frequencies					
##		2. ht: Height vector used for combining the rate matrices					
######################################################################################							
##v Gotta modify this one later on
							
combinePiVectors = function(eMatrix, allPisVect, topRows=1, isRev, distMethod							
				, distCalc) {			
							
	## Determine the two closest edges based on pi-vectors						
							
	uniquePis = unique(allPisVect)						
	uniquePiIndexes = match(uniquePis, allPisVect)						
							
	## Use only pi-vectors						
	if(distCalc == 1) {						
							
		if(!isRev) {					
			estPiVect = getPiVector(eMatrix)				
			piMatrix = estPiVect[uniquePiIndexes, ]				
							
		} 					
		else {					
			piMatrix = eMatrix[uniquePiIndexes, (s+1):(s+r)]	#v				
		}					
	}						
							
	##v Use a vector of all (s+r) / (colsirrev) elements if the matrices are symmetric						
	if(distCalc == 2) {						
							
		tempRateMatrix = getNormalizedValues(eMatrix, isReversible)	#v					
							
		edgesOfInterest = colsirrev					#v					
		if(isRev) edgesOfInterest = r+s					#v
							
		piMatrix = tempRateMatrix[uniquePiIndexes, 1:edgesOfInterest]					
	}						
							
	## Use matrix norms to compute distances between matrices						
	## TO BE WRITTEN						
							
	## Determine the number of different heights to try for vectors						
	if(distMethod == "Euclidean") {						
		heightDiff = hclust(dist(piMatrix))$height					
		heightMatrix = hclust(dist(piMatrix))$merge					
	}						
							
	if(distMethod == "Aitchison") {						
							
		stopifnot(distCalc == 1)					
							
		modDist = as.dist(aitchDist(piMatrix))					
		heightDiff = hclust(modDist)$height					
		heightMatrix = hclust(modDist)$merge					
	}						
							
	## Ensure that group indexes are non-numeric						
	grpIndexes = paste("X", allPisVect, sep="")						
							
	for(i in 1:topRows) {						
		elementA = heightMatrix[i, 1]					
		elementB = heightMatrix[i, 2]					
							
		## Substitute for element A					
		if(elementA < 0) {					
			temp = grpIndexes[uniquePiIndexes[abs(elementA)]]				
			grpIndexes = gsub(paste("^",temp,"$", sep=""), i, grpIndexes)				
		}					
		else {					
			grpIndexes = gsub(paste("^",elementA,"$", sep=""), i, grpIndexes)				
		}					
							
		## Substitute for element B					
		if(elementB < 0) {					
			temp = grpIndexes[uniquePiIndexes[abs(elementB)]]				
			grpIndexes = gsub(paste("^",temp,"$", sep=""), i, grpIndexes)				
		}					
		else {					
			grpIndexes = gsub(paste("^",elementB,"$", sep=""), i, grpIndexes)				
		}					
							
	} ## The top rows of the height matrix have been considered						
							
	## Re-set the values for grouping pi-vectors						
	xIndexes = grep("X", grpIndexes)						
							
	if(length(xIndexes) > 0) {						
							
		## Consider the elements that do not have "X" prepended to them					
		labels = sort(unique(as.numeric(grpIndexes[-xIndexes])))					
		for(i in 1:length(labels)) grpIndexes = gsub(paste("^",labels[i],"$", sep=""), i, grpIndexes)					
							
		## Consider the elements that have an "X" prepended to them					
		currentI = length(labels)					
		labels = sort(unique(grpIndexes[xIndexes]))					
		for(i in 1:length(labels)) grpIndexes = gsub(paste("^",labels[i],"$", sep=""), currentI+i, grpIndexes)					
							
	}						
	else {						
		labels = sort(unique(as.numeric(grpIndexes)))					
		for(i in 1:length(labels)) grpIndexes = gsub(paste("^",labels[i],"$", sep=""), i, grpIndexes)					
	}						
							
	grpIndexes = c(length(grpIndexes)+1, grpIndexes)						
							
	return(list(grps=grpIndexes, ht=heightDiff))						
}							
							
######################################################################################							
##	Reduce the number of rate matrices based on shortness of an edge such that 						
##	the shortest edge and its parent edge are assigned the same rate matrix 						
##							
##	Input						
##v		eMatrix: Matrix containing 2n-2 rows and s+r+1 columns with the columns 					
##v			   (s+1):(s+r) corresponding to pi-vectors				
##		piGrps: Current assignment of rate matrices to different edges					
##				excluding the root			
##		isRev: TRUE/FALSE					
##							
##	Output						
##		Vector of rate matrices obtained using "zero" edge					
######################################################################################							
							
getZeroEdge = function(eMatrix, piGrps, isRev) {							
							
	numEdges = dim(eMatrix)[1]						
	qMatrix = NULL						
							
	paramMatrix <<- eMatrix	## Ensure that the correct paramMatrix is used					
							
	## Obtain the edge lengths 						
	if(isRev) edgeLengthCol = s+r+1							#v					
	else edgeLengthCol = colsirrev+1						#v						
							
	edgeVector = getNormalizedValues(paramMatrix, isReversible=isRev)[, edgeLengthCol]	#v I shortened to isRev. Let us see if it works						
							
	## Consider only the real part of normalized edge lengths						
	if(!isRev) {						
		edgeVector = as.character(edgeVector)					
		edgeVector = as.numeric(gsub("\\+0i", "", edgeVector))					
	}						
							
	## Consider all edges except the edges linked to the root						
	for(i in 1:(numEdges-2)) {						
		parent = ceiling(i/2)					
		parentEdge = match(parent, mergeVector)					
							
		if(piGrps[i] != piGrps[parentEdge]) qMatrix = c(qMatrix, edgeVector[i])					
		else qMatrix = c(qMatrix, 1)					
							
	} ## All edges have been considered						
							
	grpIndexes = NULL						
							
	edge = match(min(qMatrix) ,qMatrix)						
	parentEdge = match(ceiling(edge/2), mergeVector)						
							
	grpIndexes = gsub(paste("^", piGrps[edge], "$", sep="")						
				, piGrps[parentEdge]			
				, piGrps)			
							
	## Re-set the labels						
	labels = sort(unique(as.numeric(grpIndexes)))						
	for(i in 1:length(labels)) grpIndexes = gsub(paste("^",labels[i],"$", sep=""), i, grpIndexes)						
							
	grpIndexes = c(length(grpIndexes)+1, grpIndexes)						
							
	return(grpIndexes)						
}							
							
######################################################################################							
##	Main function for estimating the minimum number of rate matrices required						
##	to model the evolutionary process						
##							
##	Input						
##		ioDir: Directory for input and output files					
##		seqFile: Sequence file					
##		numSites: Number of sites per species					
##		pFile: Input parameters file					
##		treeTop: Tree topology file containing a rooted Newick format tree					
##		initF: Initial value of F-vector at the root					
##		isRev: Is the matrix reversible? Currently supports only reversible matrices					
##		numItr: Number of iterations					
##		numTaxa: Number of taxa					
##		threshHeightDiff: Cut-off value for height matrix. If height difference					
##					between the top two rows is less than the cut-off		
##					value, then both rows are considered for reducing 		
##					the number of distinct rate matrices		
##		threshLogDiff: Cut-off value for difference in log-likelihood for 					
##				   rate matrices obtained using pi-vectors and zero edge.			
##				   If the difference is less than the cut-off value, then			
##				   both sets of rate matrices are analyzed to determine 			
##				   the best set			
##		distMethod: Method used for computing pairwise distances					
##				Acceptable values include "Euclidean" and "Aitchison"			
##							
##	Output						
##		A set of files for further analysis and obtaining the best set of 					
##		rate matrices as per AIC or another criterion					
##		1. Log.txt: Contains log-likelihoods for all the evaluated cases					
##		2. FileX_Y.txt: Contains the input and output parameters for a particular					
##		   			case. Here, X>=0 and 1<=Y<=max(rate_matrices)		
######################################################################################							
							
iterativeModelSearch = function(ioDir, seqFile, numSites, pFile, treeTop							
				, initF, isRev, numItr, numTaxa			
				, threshHeightDiff = 0.025, threshLogDiff = 1			
				, distMethod = "Euclidean"			
				, distCalc = 1) {

							
	## Estimate parameters with the maximum number of parameters						
	grpSmat = 1:(2*numTaxa - 2)						
	grpPiVector = c(2*numTaxa - 1, grpSmat)						
							
	fileName = "MaxRates.txt"						
							
	initializeGlobalParam(inOutDir = ioDir						
				, seqFile					#v Since I removed the option seqLength for the original "initializeGlobalParam" function. I have to delete it here as well. Originally, It lied between seqFile and paramFileName with the name numSites
				, paramFileName = pFile			
				, newickTreeTopFile = treeTop			
				, initRootVect = initF			
				, isReversible = isRev)			
							
	estimateMleParam(inOutDir = ioDir, outFile = fileName, grpSmat = grpSmat, grpPiVector = grpPiVector						
			   , numIterations = numItr, invProbSame = FALSE				
			   , computeInv = TRUE, beta = 0, verbose = TRUE				
			   , isReversible = isRev)							#v			
							
	startRate = 2*numTaxa - 3	 ## (2n-2) - 1 					
	cont = TRUE						
	index = 0						
							
	while(cont) {						
							
		if(!file.exists("StartPoint.txt")) index = 0					
		else {					
			piVect = scan("StartPoint.txt", sep="\t", skip=index, nlines=1)				
			startRate = length(unique(piVect)) - 2				
							
			## Obtain the new file				
			index = index + 1				
			fileName = paste("File", index, "_", startRate+1, ".txt", sep="")				
							
			initializeGlobalParam(inOutDir = ioDir				
						, seqFile				#v Since I removed the option seqLength for the original "initializeGlobalParam" function. I have to delete it here as well. Originally, It lied between seqFile and paramFileName with the name numSites
						, paramFileName = pFile	
						, newickTreeTopFile = treeTop	
						, initRootVect = initF	
						, isReversible = isRev)	
							
			grpPiVector = piVect				
			grpSmat = grpPiVector[-1]				
							
			estimateMleParam(inOutDir = ioDir, outFile = fileName, grpSmat = grpSmat, grpPiVector = grpPiVector				
						, numIterations = numItr, invProbSame = FALSE	
						, computeInv = TRUE, beta = 0, verbose = TRUE	
						, isReversible = isRev)					#v
							
		}					
							
		for(itr in startRate:1) {					
							
			x = reduceModelComplexity(fileName, threshHeightDiff, isRev, distMethod, distCalc)				
							
			fileNameNoExt = paste("File", index, "_", itr, sep="")				
			fileName = paste(fileNameNoExt, ".txt", sep="")				
							
			## Rate matrix reduction based on pi-vectors				
							
			grpPiVector = as.numeric(x[[1]])				
			grpSmat = grpPiVector[-1]				
			fileName1 = paste(fileNameNoExt, "NonZero.txt", sep="")				
							
			initializeGlobalParam(inOutDir = ioDir				
						, seqFile				#v Since I removed the option seqLength for the original "initializeGlobalParam" function. I have to delete it here as well. Originally, It lied between seqFile and paramFileName with the name numSites
						, paramFileName = pFile	
						, newickTreeTopFile = treeTop	
						, initRootVect = initF	
						, isReversible = isRev)	
							
			log1 = estimateMleParam(inOutDir = ioDir, outFile = fileName1, grpSmat = grpSmat, grpPiVector = grpPiVector				
						, numIterations = numItr, invProbSame = FALSE	
						, computeInv = TRUE, beta = 0, verbose = TRUE	
						, isReversible = isRev)					#v
							
			## Rate matrix reduction based on edge length				
							
			grpPiVector = as.numeric(x[[2]])				
			grpSmat = grpPiVector[-1]				
			fileName2 = paste(fileNameNoExt, "Zero.txt", sep="")				
							
			initializeGlobalParam(inOutDir = ioDir				
						, seqFile				#v Since I removed the option seqLength for the original "initializeGlobalParam" function. I have to delete it here as well. Originally, It lied between seqFile and paramFileName with the name numSites	
						, paramFileName = pFile	
						, newickTreeTopFile = treeTop	
						, initRootVect = initF	
						, isReversible = isRev)	
							
			log2 = estimateMleParam(inOutDir = ioDir, outFile = fileName2, grpSmat = grpSmat, grpPiVector = grpPiVector				
						, numIterations = numItr, invProbSame = FALSE	
						, computeInv = TRUE, beta = 0, verbose = TRUE	
						, isReversible = isRev)					#v
							
			write(itr, "Log.txt", append=TRUE)				
			write(paste("Log1 = ", log1, "Log2 = ", log2, sep="\t"), "Log.txt", append=TRUE)				
							
			if(log2 > log1) {				
				unlink(fileName1)			
				file.rename(fileName2, fileName)			
			}				
			else {				
				unlink(fileName2)			
				file.rename(fileName1, fileName)			
			}				
							
			if((abs(log1 - log2) < threshLogDiff) && (abs(log1 - log2) > 0.09)) {				
							
				if(log2 < log1) info = as.numeric(x[[2]])			
				else info = as.numeric(x[[1]])			
							
				compareAndWrite(info)			
			}				
							
		} ## All rate matrices considered starting from a given value					
							
		cont = FALSE					
							
		if(file.exists("StartPoint.txt")) {					
			temp = scan("StartPoint.txt", sep="?")				
			if(length(temp) > index) cont = TRUE				
		}					
							
	} ## All indexes have been considered						
							
}							
							
######################################################################################							
##	Add a rate matrix set to the file provided the the rate matrix set is not already 						
##	present in the file. If the input rate matrix set is exactly the same as one 						
##	of the RMS already in the file, the new RMS is not added. This function could 						
##	be modified such that <1,1,2,2,2,2> is treated as equal to <2,2,1,1,1,1>						
##							
##	Input						
##		dataToWrite: New RMS					
##		outputFile: File name where the RMS are stored					
##							
##	Output						
##		The revised outputFile					
######################################################################################							
							
compareAndWrite = function(dataToWrite, outputFile="StartPoint.txt") {							
							
	if(file.exists(outputFile)) {						
							
		tempRows =  read.table(outputFile, sep="\t", header=FALSE)					
							
		if(sum(apply(tempRows, 1, isEqual, dataToWrite)) == 0) {					
			write(dataToWrite, outputFile, sep="\t", ncolumns=length(dataToWrite), append=TRUE)				
		}					
							
	}						
	else {						
		write(dataToWrite, outputFile, sep="\t", ncolumns=length(dataToWrite), append=TRUE)					
	}						
							
}							
							
##########################################							
##	Check whether two vector are equal						
##							
##	Input						
##		x: First RMS					
##		vect: Second RMS					
##							
##	Output						
##		True/False value					
##########################################							
							
isEqual = function(x, vect) {							
							
	masterValue = as.numeric(vect)[-1]						
	presentValue = as.numeric(x)[-1]						
							
	currentCfg = getNewConfig(masterValue, presentValue)						
	diff = masterValue - currentCfg						
							
	if(sum(diff) == 0) return(TRUE)						
	return(FALSE)						
							
	#return(all(x == vect))						
}							
							
############################################################							
##	Obtaint the pi-vector for a non-symmetric rate matrix 						
##							
##	Input						
##		rateMatrix: Non-symmetric rate matrix					
##							
##	Output						
##		A vector of equilibrium frequencies					
############################################################							
							
getPiVector = function(rateMatrix) {
							
	piVectMatrix = matrix(0, nrow=1, ncol=r)		#v

	for(i in 1:dim(rateMatrix)[1]) {
		R.pars = rateMatrix[i, 1:colsirrev]		#v
		R = matrix(R.pars, r, r-1)			#v
		R = rbind(R, 0)					#v
		R = as.vector(R)				#v
		R = matrix(c(0, R), r, r)			#v
		R = R - diag(apply(R, 1, sum))
		eigR = eigen(t(R))				## Obtain the left eigen vector of R

		## Save the first eigen vector after scaling
		piVectMatrix = rbind(piVectMatrix, eigR$vector[, r]/sum(eigR$vector[, r]))	#v why does he takes the last eigen vector?

	} ## All rows have been considered

	return(piVectMatrix[-1, ])

}

							
############################################################							
##	Compute all pairwise Aitchison's distance						
##							
##	Input						
##		inMatrix: Input matrix					
##							
##	Output						
##		A matrix of pairwise distances					
############################################################							
							
aitchDist = function(inMatrix) {							
							
	numVect = dim(inMatrix)[1]						
	numCol = dim(inMatrix)[2]						
							
	distMatrix = matrix(0, nrow=numVect, ncol=numVect)						
							
	for(i in 1:numVect) {						
							
		for(j in 1:numVect) {					
							
			if(i == j) distMatrix[i, j] = 0				
			if(i > j) distMatrix[i, j] = distMatrix[j, i]				
			if(i < j) distMatrix[i, j] = getDist(inMatrix[i, ], inMatrix[j, ])				
							
		} ## All cols considered					
							
	} ## All rows considered						
							
	return(distMatrix)						
}							
							
############################################################							
##	Calculation of Aitchison's distance between two 						
##	vectors						
##							
##	Input						
##		vectOne: First vector					
##		vectTwo: Second vector					
##							
##	Output						
##		A scalar distance between the two vectors					
############################################################							
							
getDist = function(vectOne, vectTwo) {							
							
	numElements = length(vectOne)						
	distVal = 0						
							
	for(i in 1:numElements) {						
							
		for(j in 1:numElements) {					
							
			distVal = distVal + (log(vectOne[i]/vectOne[j]) - log(vectTwo[i]/vectTwo[j]))^2				
							
		} ## Inner loop ends					
							
	} ## Outer loop ends						
							
	distVal = sqrt(distVal/(2*numElements))						
							
	return(distVal)						
							
}							
							
#####################################################################							
##	Return the parameter estimates as a list of values						
##							
##	Input						
##		fileName: Output file name					
##		numTaxa: Number of taxa					
##		isRev: TRUE -> reversible rate matrix					
##		       FALSE -> non-reversible rate matrix					
##							
##	Output						
##		A list containing the following elements 					
##		logL: log likelihood 					
##		unconsLogL: Unconstrained log likelihood					
##		rateMatrix: Edge-specific parameters including 					
##			    edge lengths				
##		rootVect: Marginal frequencies at the root node					
###		probInv: Probability of a site being invariable					
##		condProbGivenInv: P(X|inv), where X = {A, C, G, T} 					
#####################################################################							
							
obtainParameterEstimates = function(fileName, numTaxa, isRev) {							
							
	if(isRev) totalCol = s+r+1						#v
	else totalCol = colsirrev+1						#v
							
	outMatrix = matrix(0, nrow=1, ncol=totalCol)						
							
	lineToSkip = (numTaxa - 1) + 10 + 2*(numTaxa - 1) + 11						
							
	tempLine = scan(fileName, skip = lineToSkip, nlines = 1, sep="\t")						
	logLikelihood = tempLine[1]						
	unconsLogLikelihood = tempLine[2]						
							
	lineToSkip = lineToSkip + 2						
							
	for(i in 1:(2*numTaxa - 2)) {						
							
		tempLine = scan(fileName, what="character", skip = lineToSkip, nlines = 1, sep="\t")					
		tempLine = removeIota(tempLine)					
							
		outMatrix = rbind(outMatrix, tempLine)					
		lineToSkip = lineToSkip + 1					
							
	} ## All edges have been considered						
							
	lineToSkip = lineToSkip + 1						
	rootVector = scan(fileName, skip = lineToSkip, nlines = 1, sep="\t")						
							
	lineToSkip = lineToSkip + 2						
	beta = scan(fileName, skip = lineToSkip, nlines = 1, sep="\t")						
							
	lineToSkip = lineToSkip + 2						
	probGivenInv = scan(fileName, skip = lineToSkip, nlines = 1, sep="\t")						
							
	outMatrix = outMatrix[-1, ]						
							
	return(list(logL = logLikelihood						
			, unconsLogL = unconsLogLikelihood				
			, rateMatrix = outMatrix				
			, rootVect = rootVector				
			, probInv = beta				
			, condProbGivenInv = probGivenInv))				
}							
							
#####################################################################							
##	Function for removing the complex part when a 12-parameter						
##	rate matrix is used						
##							
##	Input						
##		inputVector: Edge specific parameters					
##							
##	Output						
##		A vector of real numbers					
#####################################################################							
							
removeIota = function(inputVector) {							
							
	newInputVector = NULL						
							
	for(i in 1:length(inputVector)) newInputVector = c(newInputVector, gsub("\\+0i", "", inputVector[i]))						
							
	return(as.double(newInputVector))						
}							
							
							
###################################################################################################							
##	Function for optimal model search for one or more datasets						
##							
##	Input						
##		ioDir: Directory for input and output files					
##		seqFile: Sequence file					
##		numSites: Number of sites per species					
##		pFile: Input parameters file					
##		treeTop: Tree topology file containing a rooted Newick format tree					
##		initF: Initial value of F-vector at the root					
##		isRev: Is the matrix reversible? Currently supports only reversible matrices					
##		numItr: Number of iterations					
##		numTaxa: Number of taxa					
##		threshHeightDiff: Cut-off value for height matrix. If height difference					
##					between the top two rows is less than the cut-off		
##					value, then both rows are considered for reducing 		
##					the number of distinct rate matrices		
##		threshLogDiff: Cut-off value for difference in log-likelihood for 					
##				   rate matrices obtained using pi-vectors and zero edge.			
##				   If the difference is less than the cut-off value, then			
##				   both sets of rate matrices are analyzed to determine 			
##				   the best set			
##		distMethod: Method used for computing pairwise distances					
##				Acceptable values include "Euclidean" and "Aitchison"			
##		numReplicates: Number of replicates					
##							
##							
##	Output						
##		For each data set, a directory is created with the the relevant files					
###################################################################################################							
							
modelSearch = function(ioDir, seqFile, numSites, pFile, treeTop							
				, initF, isRev, numItr, numTaxa			
				, threshHeightDiff = 0.025, threshLogDiff = 1			
				, distMethod = "Euclidean"			
				, distCalc = 1			
				, numReplicates) {			

	setwd(ioDir) 	#v										
	srcFileName = paste(ioDir, seqFile, sep="/")						
							
	for(repIndex in 1:numReplicates) {						
							
		skipLines = (repIndex-1) * (numTaxa + 2)					
							
		temp = scan(srcFileName, what='character', sep="?", skip=skipLines, nlines=numTaxa + 1)					
		write(temp, "Replicate.txt")					
							
		## Heuristic to search for the best set of rate matrices					
		iterativeModelSearch(ioDir, seqFile, numSites, pFile, treeTop, initF					
				, isRev, numItr, numTaxa, threshHeightDiff, threshLogDiff			
				, distMethod, distCalc)	
							
		## Save replicate files in a separate directory					
		newDir = paste(ioDir, "/", repIndex, sep="")					
		dir.create(newDir)					
							
		filesToMove = dir(pattern="^File") 			#v "^File"					
		filesToMove = c(filesToMove, "MaxRates.txt", "Log.txt", "StartPoint.txt")					
		file.copy(from = filesToMove, to = newDir)					
		file.remove(filesToMove)					
							
	} ## All replicates have been considered						
							
}


##########################################################
##	Obtain the bias-corrected AIC value
## 
##	Input
##		logVal: Log-likelihood value
##		numRates: Number of distint rate matrices
##		nSites: Number of sites
##		nTaxa: Number of taxa
##		isRev: TRUE/FALSE
##
##	Output
##		AIC value
##########################################################

getAIC = function(logVal, numRates, nsites, nTaxa, isRev) {
	
	if(isRev) {
		maxDeg = (2*nTaxa - 2)*(r+s-1) + (r-1) + r		#v Formula is: number of eges x number of free parameters into the rate matrix + number of parameters in the root vector + alpha=1 (this goes for the ration of variant sites) + r-1
		df = maxDeg - (2*nTaxa - 2 - numRates)*(r+s-2)		#v 
	}
	else {
		maxDeg = (2*nTaxa - 2)*colsirrev + (r-1) + r		#v Formula is: number of eges x number of free parameters into the rate matrix + number of parameters in the root vector + alpha=1 (this goes for the ration of variant sites) + r-1 
		df = maxDeg - (2*nTaxa - 2 - numRates)*(colsirrev-1)	#v 
	}
	
	aic = -2*logVal + 2*df

	if((nsites/df) <= 40) {
		aic = aic + 2*df*(df+1)/(nsites - df - 1)	
	}

	return(aic)

}

##############################################################
##	Obtain the best AIC value per bootstrap replicate
##
##	Input
##		masterDir: Input directory
##		numReplicates: Number of bootstrap replicates
##		numSites: Number of sites
##		numSpecies: Number of taxa
##		isRev; TRUE/FALSE
##
##	Output
##		- A file containing the best AIC value per 
##		  rate per bootstrap replicate
##		- A file containing the best AIC value over 
##		  all possible rates per bootstrap replicate
##############################################################

getSummary = function(masterDir, numReplicates, numSites, numSpecies, isRev=TRUE) {
	
	setwd(masterDir) #v
	numRates = 2*numSpecies - 2
	bestAicPerReplicate = matrix("X", nrow=numReplicates, ncol=5)

	for(replicates in 1:numReplicates) {

		repDir = paste(masterDir, "/", replicates, sep="")

		if(length(dir(repDir)) > 0) {

			replicateMatrix = getBestAicPerRate(mastDir = masterDir
									, nTaxa = numSpecies
									, nSites = numSites
									, bsDataSet = replicates
									, isRev)

			## Save best AIC per rate per replicate in a file
			writeToFile(replicateMatrix, replicates, masterDir)

			## Obtain the overall best AIC per replicate
			allAics = as.numeric(replicateMatrix[, 1])
			minAic = min(allAics)
			index = match(minAic, allAics)

			bestAicPerReplicate[replicates, ] = c(replicateMatrix[index, ], index)
		}
		
	} ## All replicates have been considered

	writeToFile(bestAicPerReplicate, 0, masterDir)

}

##########################################################
##	Write the AIC values to a file
##
##	Input
##		dataMatrix: Data matrix 
##		fileKey: Bootstrap replicate number
##		outDir: Output directory
##
##	Output
##		Text file
##########################################################

writeToFile = function(dataMatrix, fileKey, outDir) {

	numCols = dim(dataMatrix)[2]	

	if(fileKey != 0) fileName = paste("Best_AIC_Per_Rate_BS_", fileKey, ".txt", sep="")
	else fileName = "Best_AIC_Per_BS.txt"

	for(i in 1:dim(dataMatrix)[1]) {
		write(dataMatrix[i, ], fileName, ncolumns=numCols, sep="\t", append=TRUE)
	} ## All rows have been considered

}

##############################################################
##	For a given data set, obtain the best AIC per U, 
##	where U denotes the number of distinct rate matrices
##
##	Input
##		mastDir: Input directory
##		numTaxa: Number of taxa
##		numSites: Number of sites
##		bsDataSet: Bootstrap replicate number
##		isRev: TRUE/FALSE
##
##	Output
##		A Zx4 matrix such that Z denotes the value of 
##		U in increasing order and columns denote the AIC
##		value, log-likelihood, pattern and filename 
##		in that order
##############################################################

getBestAicPerRate = function(mastDir, nTaxa, nSites, bsDataSet, isRev) {

	replicates = bsDataSet
	setwd(mastDir)	#v
	currentDir = paste(mastDir, "/", replicates, sep="")
	filesToCheck = paste(currentDir, dir(path=currentDir, pattern="File"), sep="/")
	filesToCheck = c(paste(currentDir, "MaxRates.txt", sep="/"), filesToCheck)

	maxRates = 2*nTaxa - 2

	## AIC, log-likelihood, pattern, file name
	aicMatrix = matrix("X", nrow=maxRates, ncol=4)

	for(fileIndex in 1:length(filesToCheck)) {

		## Read the S-matrices
		linesToSkip = 2 + maxRates/2 + 4 
		numSmatrices = scan(filesToCheck[fileIndex], what='character', skip=linesToSkip, nlines=1, sep="\t")
		numUniqueMatrices = length(unique(numSmatrices))
		
		## Obtain the log-likelihood and uncosntrained log-likelihood
		linesToSkip = linesToSkip + 4 + maxRates + 11
		logL = scan(filesToCheck[fileIndex], , what='character', skip=linesToSkip, nlines=1, sep="\t")
		logL = logL[1]
		
		currentAIC = getAIC(as.numeric(logL), numUniqueMatrices, nSites, nTaxa, isRev)

		storedAIC = aicMatrix[numUniqueMatrices, 1]

		if((storedAIC == "X") || (as.numeric(storedAIC) > as.numeric(currentAIC))) {
			aicMatrix[numUniqueMatrices, ] = c(currentAIC, logL, toString(numSmatrices), filesToCheck[fileIndex])
		}

	} ## All files have been considered

	return(aicMatrix)

}

################################################################
##	Obtain the number of distinct rate matrix sets 
##	that are analyzed
##
##	Input
##		masterDir: Input directory
##		numReplicates: Number of replicates
##		maxParam: Maximum number of distinct rate matrices 
##
##	Output
##		Number of distinct RMS analyzed
################################################################

numFilesAnalyzed = function(masterDir, numReplicates, maxParam) {

	numFiles = NULL

	for(replicates in 1:numReplicates) {

		uniqueRMS = list()
		for(i in 1:(maxParam-2)) uniqueRMS[[i]] = 0

		currentDir = paste(masterDir, "/", replicates, sep="")
		filesToCheck = paste(currentDir, dir(path=currentDir, pattern="Log"), sep="/")
	
		inputData = scan(filesToCheck, sep="?", what="character")
		lineNum = 1		

		while(lineNum < length(inputData)) {

			#print(lineNum)
			numUniqueSmat = as.integer(inputData[lineNum]) - 1	## Subtract 1 as U=1 is the GTR model
				
			if(numUniqueSmat > 0) {

				temp = unlist(strsplit(inputData[lineNum + 1], "\t"))
				logL1 = as.numeric(temp[2])
				logL2 = as.numeric(temp[4])

				currentElements = uniqueRMS[[numUniqueSmat]]
				if(min(abs(logL1 - currentElements)) >= 0.1) uniqueRMS[[numUniqueSmat]] = c(currentElements, logL1)

				currentElements = uniqueRMS[[numUniqueSmat]]
				if(min(abs(logL2- currentElements)) >= 0.1) uniqueRMS[[numUniqueSmat]] = c(currentElements, logL2)

			} ## RMS has two or more unique rate matrices

			lineNum = lineNum + 2

		} ## All lines have been read

		temp = length(unlist(uniqueRMS)) - (maxParam - 2)
		numFiles = c(numFiles, temp)

	} ## All replicates have been considered

	## For each replicate, two additional files are analyzed - 
	## (1) A unique rate matrix per edge and (2) the same rate matrix per edge
	return(numFiles + 2)
}

################################################################
## Obtaint the AIC values for a set of files
##
## Input
##	inputFileNames: Vector of file names
##	nTaxa: Number of taxa
##	nSites: Number of sites
##
## Output
##	A vector of AIC values, one per file	
################################################################

getAicPerGroup = function(inputFileNames, nTaxa, nSites) {

	currentAic = NULL
	maxRates = 2*nTaxa - 2

	for(fileIndex in 1:length(inputFileNames)) {

		## Read the S-matrices
		linesToSkip = 2 + maxRates/2 + 4 
		numSmatrices = scan(inputFileNames[fileIndex], what='character', skip=linesToSkip, nlines=1, sep="\t")
		numUniqueMatrices = length(unique(numSmatrices))
		
		## Obtain the log-likelihood and uncosntrained log-likelihood
		linesToSkip = linesToSkip + 4 + maxRates + 11
		logL = scan(inputFileNames[fileIndex], , what='character', skip=linesToSkip, nlines=1, sep="\t")
		logL = logL[1]
		
		currentAic = c(currentAic
					, getAIC(as.numeric(logL), numUniqueMatrices, nSites, nTaxa))

	} ## All files have been considered

	return(currentAic)

}

######################################################
## Eliminate whitespaces as prefix/suffix. This 
## function is a copy of a function of the same name 
## provided as part of LIMMA package
##
## Input
##	x: String
##
## Output
##	String without whitespaces as prefix/suffix
######################################################

trimWhiteSpace = function(x) 
{
    sub("[ \t\n\r]*$", "", sub("^[ \t\n\r]*", "", x))
}

######################################################
## Obtain the difference between two RMS. Both the RMS 
## are input as comma-separated values
##
## Input
##	currentConfig: Observed RMS
##	masterConfig: True RMS
##
## Output
##	The number of changes required to convert 
##	the observed RMS to true RMS
######################################################

getPatternIndex = function(currentConfig, masterConfig) {
	
	mstCfg = as.numeric(trimWhiteSpace(unlist(strsplit(masterConfig, split=","))))
	currentCfg = as.numeric(trimWhiteSpace(unlist(strsplit(currentConfig, split=","))))
	
	currentCfg = getNewConfig(mstCfg, currentCfg)
	diff = mstCfg - currentCfg

	indexes = length(which(diff != 0))

	return(indexes)
}

######################################################
## Obtain an RMS configuration that is as close as 
## possible to the true RMS configuration
##
## Input
##	master: True RMS
##	pat: Observed RMS
##
## Output
##	A vector representing an RMS
######################################################

getNewConfig = function(master, pat) {

	cont = TRUE
	currentPat = pat
	checkedValues = NULL

	while(cont) {

		diff = master - currentPat

		## Initialize parameters
		pos = 1
		temp = 1

		maxPos = length(diff)
		startPosVect = 1	## A vector of starting positions for runs

		i = 1			## Determine the number of elements in the list
		runList = list()	
		runList[[i]] = 1	## Number of elements per run

		while(temp < maxPos) {

			temp = pos + 1

			while((diff[temp] == diff[pos]) 
				&& (pat[temp] == pat[pos])
				&& (temp <= maxPos)) {

				runList[[i]] = runList[[i]] + 1
				temp = temp + 1

			}

			i = i + 1

			if(temp <= maxPos) {
				pos = temp
				startPosVect = c(startPosVect, pos)
				runList[[i]] = 1
			}

		} ## End of WHILE loop --> All elements of the difference vector have been considered

		runListVector = unlist(runList)
		names(runListVector) = startPosVect

		bestRun = sort(runListVector, decreasing = TRUE)
		bestStartPos = as.numeric(names(bestRun))

		## Substitue pattern from master to currenPat
		## if the value to be subsituted (in currentPat)
		## and the new value (in master) have not already 
		## been used
		removeIndexes = which(master[bestStartPos] %in% checkedValues)
		removeIndexes = c(removeIndexes, which(currentPat[bestStartPos] %in% checkedValues))
		if(length(removeIndexes) > 0) bestRun = bestRun[-removeIndexes]

		if(length(bestRun) > 0) {
		
			x = modifyPattern(currentPat, master
					, checkedValues
					, as.numeric(names(bestRun)[1]))

			currentPat = x[[1]]
			checkedValues = x[[2]]
		}
		else { cont = FALSE }

	} ## Current RMS as close as possible to true RMS
	
	return(currentPat)
	
}

##########################################################
## Convert source RMS into a form that is close to the 
## true RMS
##
## Input
##	pat: Source RMS
##	mstr: True RMS
##	valVector: Elements mapped between true RMS and 
##		   source RMS
##	posCheck: Start position of a run
##
## Output
##	A list containing two components
##	1. cfg: A vector representing an RMS
##	2. valueVector: Elements that have been mapped 
##			between true RMS and source RMS
##########################################################

modifyPattern = function(pat, mstr, valVector, posCheck) {

	## Change the labels
	temp = paste("X", pat, sep="")
	
	if(length(valVector) > 0) {
	
		for(i in 1:length(valVector)) {
			
			patToSearch = paste("X", valVector[i], sep="")
			
			## Code added on 11-Jan-2011
			## Ensure that the exact pattern is searched
			patToSearch = paste("^", patToSearch, "$", sep="")

			temp = gsub(patToSearch, valVector[i], temp)
		}
	
	} ## Pre-selected labels have been set

	newValue = mstr[posCheck]
	oldValue = paste("X", pat[posCheck], sep="")
	
	temp = gsub(oldValue, newValue, temp)
	tempValVector = c(valVector, newValue)
	
	## Convert the remaining X's to numbers
	
	valuesToAdd = 1:length(unique(pat))
	valuesToAdd = setdiff(valuesToAdd, tempValVector)
		
	xToChange = unique(temp)
	xToChange = setdiff(xToChange, tempValVector)
	
	if(length(xToChange) > 0) {
		
		for(i in 1:length(xToChange)) temp = gsub(xToChange[i], valuesToAdd[i], temp)
		
	} ## All X-values have been converted to numeric values

	return(list(cfg = as.numeric(temp), valueVector = tempValVector))	
}

##############################################################################
##Convert Newick tree file and Newick with rate matrix arrangement to
##merge matrix and rate matrix 
##Input
##treeFile contains one line of form ((a,b),(c,d))
##rateTreeFile contains one line of form ((a_1,b_1)2,(c_1,d_2)2)
##Output mergeMat, rateTree, leaves.
##corresponding to merge  -1 -2 -3 -4 1 2 and grouping 1 1 1 2 2 2
############################################################################
newickFormsToMergeRateMats=function(treeFile,rateTreeFile){
	bt=newickFormatToMergeMatrix(treeFile)
        btl=gsub("_","-",bt$leaves)
	btV=c(t(bt$dMat))
	btrates=scan(rateTreeFile,what="character")
	for(i in 1:length(btl))btrates=sub(bt$leaves[i],btl[i],btrates)
	aa=gsub("\\(","",btrates)
	aa=unlist(strsplit(aa,"_"))
	aa=unlist(strsplit(aa,","))
	aa=gsub(")",":):",aa)
	aa=unlist(strsplit(aa,":"))
	aa=aa[-length(aa)]
	In=which(aa==")")
	aa[In]=1:length(In)
	for(i in 1:length(aa)){
	for(j in 1:length(btl)){
		if(aa[i]==btl[j])aa[i]=-j
		}
	}
	aa1=as.numeric(aa[1:(length(aa)/2)*2-1])
	aa2=as.numeric(aa[1:(length(aa)/2)*2])
	pMat=matrix(as.numeric(outer(aa1,btV,"==")),length(aa)/2,length(aa)/2)
	mTree=matrix(t(pMat)%*%aa1,length(aa)/4,2,byrow=TRUE)
	rTree=matrix(t(pMat)%*%aa2,length(aa)/4,2,byrow=TRUE)
	mV=c(t(mTree))
	rV=c(t(rTree))
	list(mergeMat=mTree,rateTree=rTree,leaves=bt$leaves)
}

########################################################################
##Convert merge matrix, rate matrix, leaves to newickTree withrates
##Output of form ((a_1,b_1)2,(c_1,d_2)2)
##corresponding to merge vector -1 -2 -3 -4 1 2 and grouping 1 1 1 2 2 2
#########################################################################
mergeMatrixToNewickRateFormat = function(hMat, rMat, leaves) {							

paddedRateVal = function(item,ritem) {							
							
	if(item < 0) newItem = paste("#", item, "_",ritem,"#", sep="")						
	else newItem = paste("@", item, "@", sep="")						
							
	return(newItem)						
}							
							
	numRows = dim(hMat)[1]						
							
	newickTree = paste("("						
				, paddedRateVal(hMat[numRows, 1],rMat[numRows,1])			
				, ","			
				, paddedRateVal(hMat[numRows, 2],rMat[numRows,2])			
				, ")"			
				, sep="")			
							
	for(i in (numRows-1):1) {						
							
		addRow = paste("("					
				, paddedRateVal(hMat[i, 1],rMat[i,1])			
				, ","			
				, paddedRateVal(hMat[i, 2],rMat[i,2])			
				, ")"
				, rMat[which(hMat==i,arr.ind=TRUE)]			
				, sep="")			
							
		rowPattern = paste("@", i, "@", sep="")					
		newickTree = sub(rowPattern, addRow, newickTree)					
							
	} ## The mergeMatrix has been considered						
							
	## Replace leaf indexes with leaf names						
	for(i in 1:length(leaves)) {								
		newickTree = sub(-i, leaves[i], newickTree)												
	}						
	newickTree=gsub("#","",newickTree)						
	return(newickTree)						
}
