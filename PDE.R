#  Forever and a day we've been evolving F as:
#	d F(s) /dt = -s F 
#  The goal here is to work out toy software to evolve the PDE:
#	\frac{\partial^2 F}{\partial t \partial s} =  -F(s) -s \partial F/\partial s

#One set of ss for everything
NUMSVALS=1000
C <- 0.01
Smin <- 0.01
Svals <- Smin*(1+C)^(0:(NUMSVALS-1))
# this is delta n.
dlogs <- 1
log1pc <- log(1+C)
ns <- 1:NUMSVALS
# let's choose a delta t
DELTFRAC <- .1
# how much bigger than the smallest time constant should you run?
TIMEFAC <- 20
delt <- DELTFRAC * min(1/Svals)
maxt <- TIMEFAC* min(1/Svals)
times <- seq(0,maxt,by=delt)
print(paste("min time constant is", min(1/Svals)))
print(paste("max time constant is", max(1/Svals)))
print(paste("max time is", maxt))
print(paste("delt is", delt))
print(paste("running for", length(times), "steps"))

# i want to check in on P at stages of evolution to see what's going wrotng
numPchecks=21
pchecktimes = sort(maxt*(0:(numPchecks-1))/(numPchecks-1))

Pchecks <- array(NA,c(numPchecks,length(Svals)))
Fchecks <- array(NA,c(numPchecks,length(Svals)))


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
getnewFPDEn <- function(OldF, delt, revtime=F){
	if(revtime){
		delt= -delt
	}
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

runevolve <- function(){
	# initializing F to be all 1.  We will decay.
	Fvals <<- array(1, length(Svals))
	# for going backwards let's start at maxt.  Should end up at 1
	Pvals <<- exp(-maxt*Svals)

	OldF <<- Fvals
	OldP <<- Pvals
	checktimeind <- 1
	for(time in times){

		if(time >= pchecktimes[checktimeind]){
			print(paste("found index", checktimeind, "at time", time))
			Pchecks[checktimeind,] <<- OldP
			Fchecks[checktimeind,] <<- OldF
			checktimeind <- checktimeind + 1
		}
		NewF <<- getnewFPDEn(OldF,delt)
		NewP <<- getnewFPDEn(OldP,delt,revtime=T)
		
		OldP <<- NewP
		OldF <<- NewF
	}
}

wasFokfig <- function(){
	layout(t(1:2))
	# see if this worked.  If things went according to plan, 
	#	 F(s) = e^{-maxt s} 
	# So, this plot should give a line with slope 1 and intercept zero.
	plot(NewF,exp(-maxt*Svals), type="l", lwd=5, col=rgb(0.5,0.5,1.0,0.5))
	# to compare the two
	#plot(NewP,exp(-futstopt*Svals), type="p", pch=19,
	#     col=rgb(0.5,0.5,1.0,0.5),xlim=c(0,1), ylim=c(0,1)) 
	abline(0,1)

#	Does newF end up in the same configuration that P started with?	
	plot(NewF,Pchecks[1,], type="l", lwd=5, col=rgb(0.5,0.5,1.0,0.5))
	abline(0,1)
}

evolvefig <- function(){
	layout(t(1:2))
	plot(-log(Svals), NewF,type="l", ylim=c(0,1),col="blue")
	for(i in 1:length(Pchecks[,1])){
		lines(-log(Svals),Fchecks[i,])
	}

	plot(-log(Svals), NewP,type="l", ylim=c(0,1),col="blue")
	for(i in 1:length(Pchecks[,1])){
		lines(-log(Svals), Pchecks[i,])
	}
}	

stepbackandforth <- function(starttime =TIMEFAC*min(1/Svals),
			     delt=DELTFRAC * min(1/Svals), nsteps=50,
			     forward = T 
			     ){
	StartF <- exp(-starttime * Svals)
	LeaveF <- array(NA,c(nsteps,length(Svals)))
	ReturnF <- array(NA,c(nsteps,length(Svals)))

	# let's evolve out nsteps
	LeaveF[1,] <- StartF
	out <- !forward
	for(i in 2:nsteps){
		LeaveF[i,] <- getnewFPDEn(LeaveF[i-1,], delt, rev=out)
	}
	# now lets' see if we get back to the start	

	ReturnF[nsteps,] <- LeaveF[nsteps,]
	for(i in (nsteps-1):1){
		ReturnF[i,] <- getnewFPDEn(ReturnF[i+1,], delt, rev=!out)
 	} 
	retval <- list(LeaveF,ReturnF)	
	return(retval)
}

makeoutbackfigure <- function(...){
	test <- stepbackandforth(...)
	# figure out what nsteps was
	nsteps <- length(test[[1]][,1])
	layout(t(1:2))
	plot(log(Svals), test[[1]][1,],
	     type="l", xlab="log s", ylab="F(s)",ylim=c(0,1)
	)
	lines(log(Svals), test[[1]][nsteps,], col="blue")
	plot(log(Svals), test[[2]][nsteps,],
	     type="l", xlab="log s", ylab="F(s)", ,ylim=c(0,1)
	)
	lines(log(Svals), test[[2]][1,], col="blue") 
}

