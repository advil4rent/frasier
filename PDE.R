#  Forever and a day we've been evolving F as:
#	d F(s) /dt = -s F 
#  The goal here is to work out toy software to evolve the PDE:
#	\frac{\partial^2 F}{\partial t \partial s} =  -F(s) -s \partial F/\partial s


#One set of ss for everything
NUMSVALS=1000
C <- 0.01
Smin <- 0.001
deltas <- 0.01
Svals <- Smin*(1+C)^(0:(NUMSVALS-1))
# does this work way better with linear s?
#Svals <- seq(from=Smin, to=Smin+deltas*NUMSVALS,by=deltas)
# this is delta n.
dlogs <- 1
log1pc <- log(1+C)
ns <- 1:NUMSVALS
print(paste("min time constant is", min(1/Svals)))
print(paste("max time constant is", max(1/Svals)))

# let's choose a delta t
DELTFRAC <- .1
delt <- min(DELTFRAC/Svals)
NUMSTEPS <- 100
maxt <- NUMSTEPS*delt
times <- seq(0,maxt,by=delt)

# initializing F to be all 1.  We will decay.
Fvals <- array(1, length(Svals))

#	\frac{\partial^2 F}{\partial t \partial s} =  -F(s) -s \partial F/\partial s
getnewFPDE <- function(OldF, delt){
	# Let's figure out the partial
	partial <- array(NA,length(OldF)) 
	delS <- array(NA,length(OldF)) 
	rhs <- array(NA,length(OldF))
	# now the first index, partial[2] holds OldF[3] - OldF[1]
	for(i in 2:(length(OldF))){
		delS[i] <- (Svals[i] - Svals[i-1])
		# leaving out the delS because i'm going to multiply the rhs
		# by delS in a minute
		#partial[i-1] <- (OldF[i+1] - OldF[i-1])/delS[i-1]
		partial[i] <- (OldF[i] - OldF[i-1])
		# pre-multiplying the delS onto OldF for Euler method.
		rhs[i] <- OldF[i]* delS[i] + partial[i]
	}
	#test <<- partial
	# the Rhs is multiplied by \partial t and \partial s.
	# Using Euler method
	rhs <- rhs*(-delt)
	# now we want to assign newF.  
	newF <- array(NA, length(OldF))
	# This is the slowest time constant
	newF[1] <- 1
	# going forwards
	#... should really go backwards setting F(s=\infty) = 0 as boundary
	#condition.
	for(i in 2:(length(OldF))){
#		newF[i+1] - newF[i-1] - (OldF[i+1] - OldF[i-1])  = rhs[i]
		newF[i] = rhs[i] +  newF[i-1] +  (OldF[i] - OldF[i-1])
	}	
	return(newF)
}

#\frac{\partial^2 F}{\partial t \partial s} =  -F(s) -s \partial F/\partial s
# version of the PDE that evolves as a function of cell number, rather than s.
# Cell number n = \log_{1+c} s + \log_{1+c}[s_1]
# dn/ds = 1/(s * log 1+c) ; 
# s = (1+c)^n = e^(n log 1+c); ds/dn = log (1+c) s (as it should!)
# So, \partial/\partial n = log (1+c) s \partial / \partial s =  and 
# \partial / \partial s = 1/(s log 1+c)  \partial/\partial n
# eq is:
# 1/(s * log(1+c) \frac{\partial^2 F}{\partial t \partial n} 
	# =  -F(s) - 1/log(1+c) \partial F/\partial n
# or
# \frac{\partial^2 F}{\partial t \partial n} =  -s log(1+c)F(s) -s\partial F/\partial n
getnewFPDEn <- function(OldF, delt){
	# Let's figure out the partial
	partial <- array(NA,length(OldF)) 
	
	rhs <- array(NA,length(OldF))
	# now the first index, partial[2] holds OldF[3] - OldF[1]
	for(i in 2:(length(OldF))){
		# leaving out the delS because i'm going to multiply the rhs
		# by delS in a minute
		#partial[i-1] <- (OldF[i+1] - OldF[i-1])/delS[i-1]
		partial[i] <- (OldF[i] - OldF[i-1])
		# pre-multiplying the delS onto OldF for Euler method.
		# dlogs is a constant for n
		rhs[i] <- Svals[i] * (OldF[i]* log1pc * dlogs + partial[i])
	}
	#test <<- partial
	# the Rhs is multiplied by \partial t and \partial s.
	# Using Euler method
	rhs <- rhs*(-delt)
	# now we want to assign newF.  
	newF <- array(NA, length(OldF))
	# This is the slowest time constant
	newF[1] <- 1
	# going forwards
	for(i in 2:(length(OldF))){
#		newF[i+1] - newF[i-1] - (OldF[i+1] - OldF[i-1])  = rhs[i]
		newF[i] = rhs[i] +  newF[i-1] +  (OldF[i] - OldF[i-1])
	}	
	return(newF)
}

OldF <- Fvals
for(time in times){
	NewF <- getnewFPDEn(OldF,delt)
	OldF <- NewF
}
# see if this worked.  If things went according to plan, 
#	 F(s) = e^{-maxt s} 
# So, this plot should give a line with slope 1 and intercept zero.
plot(NewF,exp(-maxt*Svals), type="p", pch=19, col=rgb(0.5,0.5,1.0,0.5))
# to compare the two
abline(0,1)
write.csv(NewF, file = "R_generated_F.csv")
