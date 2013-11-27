rm(list=ls())
gc()


require("msm")

#sdTab <- abs((bounds[,,2]-bounds[,,1])/2)
	
	
#################################
### OOP

setClass("FAOtab",representation=representation(bestTab="matrix",
                                              tables="list",
                                              iters="integer",
                                              objective="numeric",
                                              call="call"))#,
                                              #args="list"))

setMethod("show",signature("FAOtab"),function(object){
	cat("Call:\n")
	print(object@call)
	cat("\nOptimal Table: ")
	print(object@bestTab)
	cat("Number of Iterations: ")
	cat(object@iters,"\n")
	cat("Objective Function: ")
	cat(object@objective,"\n")
	})

sampleTables <- function(n0,muTab, bounds,nIter=100,N=10000,controlCol=NULL,controlRow=NULL,sdev=5,verbose=TRUE,objFun=function(tab){colSums(tab)[1]}){

	call <- match.call()
	
	argz <- lapply(2:length(call),function(i)eval(call[[i]]))
	names(argz) <- names(call)[2:length(names(call))]

	sdTab <- abs((bounds[,,2]-bounds[,,1])/2)
	
	nr<-nrow(muTab)
	nc<-ncol(muTab)
	
	okTab <- list()
	
	if(is.null(controlRow)) controlRow <- do.call(rbind,lapply(1:nr,function(i){
		if(muTab[i,nc]==0) nc <- which.max(apply(bounds[i,,],1,function(x)diff(range(x))))
		
		range(bounds[i,nc,])
		
		}))
		
	# NB: questo é SOLO se non diamo control Col come input, quindi credo vada bene una normale (é il caso nel quale non abbiamo idea dei boundaries)
	if(is.null(controlCol)) controlCol <- do.call(cbind,lapply(1:nc,function(j)range(round(rnorm(N,colSums(muTab),sqrt(colSums(sdTab^2)))))))
	
	iter <- t <- 1L
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
		#browser()
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
			
			t <- t+1L
			if(verbose)print(t)
			}
		iter <- iter + 1L
        }

      okTab <- okTab[1:uniqueT]
      bestTab <- okTab[[which.max(unlist(lapply(okTab,objFun)))]]
      

      return(new("FAOtab",bestTab=bestTab,tables=okTab,iters=iter,objective=objFun(bestTab),call=call))
      
      }




muTab <- read.csv("ExpextedValues_140.csv",row.names=1,sep=";")[-141,-c(7,8)]
n0<-read.csv("ExpextedValues_140.csv",header=TRUE,row.names=1,sep=";")[-141,8]
muCol<-read.csv("ExpextedValues_140.csv",row.names=1,sep=";")[141,-c(7,8)]

bounds <- array(NA,c(dim(muTab),2))
dimnames(bounds) <- c(dimnames(muTab),list(c("Lower","Upper")))
bounds[,,1] <- as.matrix(read.csv("lowerBounds_140.csv",row.names=1,sep=';'))[-141,]
bounds[,,2] <- as.matrix(read.csv("upperBounds_140.csv",row.names=1,sep=';'))[-141,]

lowCol<-read.csv("lowerBounds_140.csv",row.names=1,sep=';')[141,]
uppCol<-read.csv("upperBounds_140.csv",row.names=1,sep=';')[141,]
controlCol1 <- rbind(lowCol,uppCol)


# l'ultima cella é il bad boy..... questa é un po' una pezza....
controlCol1[,ncol(controlCol1)]<- controlCol1[,ncol(controlCol1)]*c(5,1/5)


	
#tab <- sampleTables(n0,muTab,bounds,sdev=5)

system.time(tab <- sampleTables(n0,muTab,bounds,controlCol=controlCol1,sdev=5))
tab

#as.matrix(lapply(tab$tabs, length))

unlist(lapply(tab@tables,function(x)attr(x,"mult")))
length(tab@tables)

Scenario <- read.csv("Scenario_140.csv",sep=";",nrows=nrow(muTab))
rownames(Scenario) <- make.unique(as.character(Scenario[,1]))
Scenario <- Scenario[,-c(1,8:9)]

rmse <- unlist(lapply(tab@tables,function(t)sqrt((mean((t-Scenario)^2)))))
rmse
rrmse <- unlist(lapply(tab@tables,function(t)sqrt((mean(((t-Scenario)/(Scenario+2*sqrt(.Machine$double.eps)))^2)))))
rrmse

order(rrmse)
min(rrmse)


# check
aa<-order(rrmse)[1]
aa
sqrt((sum(((tab[[1]][,,aa]-Scenario)/(Scenario+0.0000001))^2))/length(tab[[1]][,,aa]))

--------------------------------------------------------------------------------------------------------------------

## SBAGLIATE
#rmse <- sqrt(sum(unlist(lapply(tab[[1]],function(t)sum((t-Scenario)^2)))))
#rmse <- sqrt(unlist(lapply(tab[[1]],function(t)sum((t-Scenario)^2))))

#rrmse <- sqrt(sum(unlist(lapply(tab[[1]],function(t)sum(((t-Scenario)^2)/abs(Scenario))))))

den <- 1:3
reject <- rmse <- time <- rep(NA,length(den))

for(i in 1:length(den)){
	print(i)
#	controlCol <- do.call(cbind,lapply(1:ncol(muTab),function(j)range(round(rnorm(10000,colSums(muTab),sqrt(colSums(sdTab^2))/den[i])))))
	time[i]<-system.time(res <- sampleTables(n0,muTab,bounds,controlCol=controlCol1))[3]
	rmse[i] <- unlist(lapply(res@tables,function(t)sqrt((mean((t-Scenario)^2)))))
	reject[i] <- res[[2]]
	gc()
	}
plot(1/den,time,type="b",col="blue",main="Time vs. Strenght of the Constraint",xlab="narrowness")
points(1/den,time,col="red",pch=16)

