rm(list=ls())
gc()


#setwd("C:/Users/natalia/Documents/FAO/script")

require("msm")


muTab <- read.csv("ExpectedValues.csv",row.names=1,sep=";")
muTab

n0<-read.csv("ExpectedValues.csv",header=TRUE,row.names=1,sep=";")[,8]
n0
#n0 <- muTab[,7]
#n0 <- muTab[,8]
muTab <- muTab[,-(7:8)]
muTab

bounds <- array(NA,c(dim(muTab),2))
dimnames(bounds) <- c(dimnames(muTab),list(c("Lower","Upper")))
bounds[,,1] <- as.matrix(read.csv("lowerBounds.csv",row.names=1,sep=';'))
bounds[,,2] <- as.matrix(read.csv("upperBounds.csv",row.names=1,sep=';'))

sdTab <- abs((bounds[,,2]-bounds[,,1])/2)

sampleTables <- function(n0,muTab, bounds,nIter=100,N=10000,controlCol=NULL,controlRow=NULL,sdev=5,verbose=TRUE){
	
	sdTab <- abs((bounds[,,2]-bounds[,,1])/2)
	
	nr<-nrow(muTab)
	nc<-ncol(muTab)
	
	okTab <- list()
	
	if(is.null(controlRow)) controlRow <- do.call(rbind,lapply(1:nr,function(i){
		if(muTab[i,nc]==0) nc <- which.max(apply(bounds[i,,],1,function(x)diff(range(x))))
		
		range(bounds[i,nc,])
		
		}))
	if(is.null(controlCol)) controlCol <- do.call(cbind,lapply(1:nc,function(j)range(round(rnorm(N,colSums(muTab),sqrt(colSums(sdTab^2)))))))
	
	iter <- t <- 1
	uniqueT <- 1
	
	while(t <= nIter){
		
		sim <- do.call(rbind,lapply(1:nr,function(i){
			
			rrow <- rep(NA,nc)
			repeat{
				
				nc<-ncol(muTab)
				
				### VARSTOCK structural 0
				if(muTab[i,nc]==0){
					nc <- max(which(muTab[i,-nc]!=0))
					rrow[(1:length(rrow))>nc]<-0
					
					maxTol <- which.max(apply(bounds[i,,],1,function(x)diff(range(x))))
					
					rrow[-c(maxTol,(nc+1):length(rrow))] <- round(rtnorm(nc-1, mean=unlist(muTab[i,-c(maxTol,(nc+1):length(rrow))]), sd=sdev, lower=bounds[i,-c(maxTol,(nc+1):length(rrow)),1],upper=bounds[i,-c(maxTol,(nc+1):length(rrow)),2]),0)

					rrow[maxTol] <- (n0[i]-sum(rrow[-c(maxTol,nc+1:length(rrow))]))
					nc<-maxTol
					}else{ 
						###VARSTOCK not structural 0
						
						if(length(rrow[-nc])!=nc-1)browser()
						rrow[-nc] <- round(rtnorm(nc-1, mean=unlist(muTab[i,-nc]), sd=sdev, lower=bounds[i,-nc,1],upper=bounds[i,-nc,2]),0)
						rrow[nc] <- -(n0[i]-sum(rrow[-nc]))
						
						}
				if(rrow[nc]>=min(controlRow[i,]) & rrow[nc]<=max(controlRow[i,]))break
				
				}
				if(verbose)cat("*")
			return(rrow)
			
			}))
			
		if(verbose)cat("\n")
		
		totCol <- colSums(sim)
		cond <- sapply(1:nc,function(j)(totCol[j]>=controlCol[1,j] & totCol[j]<=controlCol[2,j]))
		
		if(all(cond)){
			if(t==1){
				okTab[[uniqueT]] <- sim
				attr(okTab[[uniqueT]],"mult") <- 0
				}
				
			dejavu <- FALSE
			
			for(k in 1:uniqueT){
				if(all(sim==okTab[[k]])){
					attr(okTab[[k]],"mult") <- attr(okTab[[k]],"mult") + 1
					break
					}
				}
				
			if(!dejavu & t<nIter){
				uniqueT <- uniqueT+1
				okTab[[uniqueT]] <- sim
				attr(okTab[[uniqueT]],"mult") <- 1
				}
			
			t <- t+1
			if(verbose)print(t)
			}
		iter <- iter + 1
        }
        
      return(list(tabs=okTab[1:uniqueT],iter=iter))
	}
	
tab <- sampleTables(n0,muTab,bounds,sdev=5)

# tab <- sampleTables(n0,muTab,bounds,sdev=2.5)


## con 4 o 5 non converge...
den <- 1:3
reject <- rmse <- time <- rep(NA,length(den))

Scenario <- read.csv("ValoriItalia.csv",sep=",",row.names=1)[,-(7:9)]

for(i in 1:length(den)){
	print(i)
	controlCol <- do.call(cbind,lapply(1:ncol(muTab),function(j)range(round(rnorm(10000,colSums(muTab),sqrt(colSums(sdTab^2))/den[i])))))
	time[i]<-system.time(res <- sampleTables(n0,muTab,bounds,controlCol=controlCol))[3]
	rmse[i] <- sqrt(sum(unlist(lapply(res[[1]],function(t)sum((t-Scenario)^2)))))
	reject[i] <- res[[2]]
	gc()
	}
plot(1/den,time,type="b",col="blue",main="Time vs. Strenght of the Constraint",xlab="narrowness")
points(1/den,time,col="red",pch=16)





