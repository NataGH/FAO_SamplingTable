rm(list=ls())
gc()


setwd("C:/Users/natalia/Documents/FAO/script")

require(msm)


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




sampleTables <- function(n0,muTab, bounds,nIter=100,N=10000,controlCol=NULL,controlRow=NULL,sdev=10){
	
	sdTab <- abs((bounds[,,2]-bounds[,,1])/2)
	
	nr<-nrow(muTab)
	nc<-ncol(muTab)
	
	okTab <- array(NA,c(nr,nc,nIter))
	#browser()
	if(is.null(controlRow)) controlRow <- do.call(rbind,lapply(1:nr,function(i){
		if(muTab[i,nc]==0) nc <- max(which(muTab[i,-nc]!=0))
#		range(round(rnorm(N,muTab[i,nc],fact*sdTab[i,nc])))
		c(bounds[i,nc,1], bounds[i,nc,2])
		}))
	if(is.null(controlCol)) controlCol <- do.call(cbind,lapply(1:nc,function(j)range(round(rnorm(N,colSums(muTab),sqrt(colSums(sdTab^2)))))))
	iter <- t <- 1
	while(t <= nIter){
		
		sim <- do.call(rbind,lapply(1:nr,function(i){
			
			rrow <- rep(NA,nc)
			repeat{
				
				#rrow[-nc] <- apply(bounds[i,-nc,],1,function(x)sample(x[1]:x[2],1))
				### VARSTOCK structural 0
				if(muTab[i,nc]==0){
					nc <- max(which(muTab[i,-nc]!=0))
					rrow[(1:length(rrow))>nc]<-0
					
					#rrow[-(nc:length(rrow))] <- round(unlist(muTab[i,-(nc:length(rrow))])+unlist(sdTab[i,-(nc:length(rrow))])*abs(rnorm(nc-1)),0)
					rrow[-(nc:length(rrow))] <- round(rtnorm(nc-1, mean=unlist(muTab[i,-(nc:length(rrow))]), sd=sdev, lower=bounds[i,-(nc:length(rrow)),1],upper=bounds[i,-(nc:length(rrow)),2]),0)

					rrow[nc] <- (n0[i]-sum(rrow[-(nc:length(rrow))]))
					}else{ ###VARSTOCK not structural 0
						
						rrow[-(nc:length(rrow))] <- round(rtnorm(nc-1, mean=unlist(muTab[i,-nc]), sd=sdev, lower=bounds[i,-nc,1],upper=bounds[i,-nc,2]),0)
						rrow[nc] <- -(n0[i]-sum(rrow[-nc]))
						#if(i==3)break
						
						}
				if(rrow[nc]>=min(controlRow[i,]) & rrow[nc]<=max(controlRow[i,]))break
				#break
				#if(i>=5)browser()
				#break
				#browser()
				}
				cat("*")
			return(rrow)
			
			}))
		cat("\n")
		totCol <- colSums(sim)
		cond <- sapply(1:nc,function(j)(totCol[j]>=controlCol[1,j] & totCol[j]<=controlCol[2,j]))
		if(all(cond)){
			okTab[,,t]<-sim
			t <- t+1
			print(t)
			}
		iter <- iter + 1
        }
        
      return(list(tabs=okTab,iter=iter))
	}
	
tab <- sampleTables(n0,muTab,bounds,sdev=100)


tab



colMeans(do.call(rbind,lapply(1:dim(tab[[1]])[3],function(t)tab[[1]][,,t])))

summary(do.call(rbind,lapply(1:dim(tab[[1]])[3],function(t)tab[[1]][,6,t])))



###
# Commenti di Luca:
# 1. righe 29-36: come vedi uno dei problemi coi vincoli di riga é 
#    che quando abbiamo zeri strutturali dobbiamo aggiustare anche i vincoli... 
#    ma le varianze sono sempre un po' off... riusciamo a tirare fuori un buon set di priors? 
#    come vedi queste sono fuori almeno di un fattore 10 (i.e. se metto 1 invece di 10 l'algoritmo non funge....)
# 2. riga 52: la TN é anche ottenibile così:
# Z~N(0,1) and X~TN(mu,sigma) iff X = mu+sigma*|Z|
# 3. é ancora un po' awkward come definisco le prior qua perché ho voluto tenere la struttura data dalla uniforme... 
# ma con nuove priors posso implementare le normali...

# Davvero, qua i vincoli hanno ancora qualcosa che non va... riusciamo magari a stringere un pochetto le varianze 
# di cella e allargare quelle di stockvar?
###








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




aa<-round(rtnorm(nc-1, mean=unlist(muTab[i,-(nc:length(rrow))]), sd=sdev, lower=bounds[i,-(nc:length(rrow)),1],upper=bounds[i,-(nc:length(rrow)),2]),0)
aa




