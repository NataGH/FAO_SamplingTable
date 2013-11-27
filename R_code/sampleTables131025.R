setwd("/Users/marcogarieri/Desktop/FAO/Natalia_Project/R_code/")

rm(list=ls())
gc()



muTab <- read.csv("ExpectedValues.csv",row.names=1)

#n0<-read.csv("ValoriItalia.csv",header=FALSE,row.names=1)
n0 <- muTab[,7]
muTab <- muTab[,-(7:9)]


bounds <- array(NA,c(dim(muTab),2))
dimnames(bounds) <- c(dimnames(muTab),list(c("Lower","Upper")))
bounds[,,1] <- as.matrix(read.csv("lowerBounds.csv",row.names=1))[,-(7:9)]
bounds[,,2] <- as.matrix(read.csv("upperBounds.csv",row.names=1))[,-(7:9)]


sampleTables <- function(n0,muTab, bounds,nIter=100,N=10000,controlCol=NULL,controlRow=NULL){
	
	sdTab <- abs((bounds[,,2]-bounds[,,1])/2)
	
	nr<-nrow(muTab)
	nc<-ncol(muTab)
	
	okTab <- array(NA,c(nr,nc,nIter))
	#browser()
	if(is.null(controlRow)) controlRow <- do.call(rbind,lapply(1:nr,function(i)range(round(rnorm(N,muTab[i,nc],sdTab[i,nc])))))
	if(is.null(controlCol)) controlCol <- do.call(cbind,lapply(1:nc,function(j)range(round(rnorm(N,colSums(muTab),sqrt(colSums(sdTab^2)))))))
	
	iter <- t <- 1
	while(t <= nIter){
		
		sim <- do.call(rbind,lapply(1:nr,function(i){
			rrow <- rep(NA,nc)
			repeat{
				rrow[1:(nc-1)] <- round(rnorm(nc-1,unlist(muTab[i,-nc]),unlist(sdTab[i,-nc])),0)
				rrow[-nc] <- apply(bounds[i,-nc,],1,function(x)sample(x[1]:x[2],1))
				rrow[nc] <- (n0[i]-sum(rrow[-nc]))
				break
#				if(rrow[nc]>=controlRow[i,1] & rrow[nc]<=controlRow[i,2]) break
				#browser()
				}
			return(rrow)
			
			}))
		
		totCol <- colSums(sim)
#		cond <- sapply(1:nc,function(j)(totCol[j]>=controlCol[1,j] & totCol[j]<=controlCol[2,j]))
#		if(all(cond)){
			okTab[,,t]<-sim
			t <- t+1
#			print(t)
#			}
		iter <- iter + 1
        }
        
      return(list(tabs=okTab,iter=iter))
	}
	
tab <- sampleTables(n0,muTab,bounds)


summary(do.call(rbind,lapply(1:dim(tab[[1]])[3],function(t)tab[[1]][,6,t])))




den <- c(0.5,1,2,5,10)
time <- rep(NA,4)
rmse <- cbind(rep(NA,4),rep(NA,4))
reject <- rep(NA,4)
Scenarios <- list(read.csv("Scenarios/Scenario2.csv",row.names=1),read.csv("Scenarios/Scenario2.csv",row.names=1))

for(i in 1:length(den)){
	print(i)
	controlCol <- do.call(cbind,lapply(1:ncol(muTab),function(j)range(round(rnorm(10000,colSums(muTab),sqrt(colSums(sdTab^2))/den[i])))))
	time[i]<-system.time(res <- sampleTables(n0,muTab,sdTab,controlCol=controlCol))[3]
	rmse[i,1] <- sqrt(sum(unlist(lapply(1:dim(res[[1]])[3],function(t)sum((res[[1]][,,t]-Scenarios[[1]])^2)))))
	rmse[i,2] <- sqrt(sum(unlist(lapply(1:dim(res[[1]])[3],function(t)sum((res[[1]][,,t]-Scenarios[[2]])^2)))))
	reject[i] <- res[[2]]
	gc()
	}
	
	
	
plot(1/den,time,type="b",col="blue",main="Time vs. Strenght of the Constraint",xlab="narrowness")
points(1/den,time,col="red",pch=16)