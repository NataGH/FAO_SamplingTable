rm(list=ls())
gc()


require(msm)

#setwd("C:/Users/nata/Desktop/FAO/Sampling Contingency Tables/R_code/var_colonne")

muTab <- read.csv("ExpectedValues.csv",row.names=1,sep=";")[-8,-c(7,8)]
muTab

n0<-read.csv("ExpectedValues.csv",header=TRUE,row.names=1,sep=";")[-8,8]
n0


muCol<-read.csv("ExpectedValues.csv",row.names=1,sep=";")[8,-c(7,8)]
muCol


bounds <- array(NA,c(dim(muTab),2))
dimnames(bounds) <- c(dimnames(muTab),list(c("Lower","Upper")))
bounds[,,1] <- as.matrix(read.csv("lowerBounds_2.csv",row.names=1,sep=';'))[-8,]
bounds[,,2] <- as.matrix(read.csv("upperBounds_2.csv",row.names=1,sep=';'))[-8,]
bounds



lowCol<-read.csv("lowerBounds.csv",row.names=1,sep=';')[8,]
lowCol

uppCol<-read.csv("upperBounds.csv",row.names=1,sep=';')[8,]
uppCol
uppCol[6]=0
uppCol

controlCol1 <- matrix(c(lowCol,uppCol),byrow=T,nrow=2)

sdTab <- abs((bounds[,,2]-bounds[,,1])/2)

sampleTables <- function(n0,muTab, bounds,nIter=100,N=10000,controlCol=NULL,controlRow=NULL,sdev=100){
	
	sdTab <- abs((bounds[,,2]-bounds[,,1])/2)
	
	nr<-nrow(muTab)
	nc<-ncol(muTab)
	
	okTab <- array(NA,c(nr,nc,nIter))
	if(is.null(controlRow)) controlRow <- do.call(rbind,lapply(1:nr,function(i){
		if(muTab[i,nc]==0) nc <- which.max(apply(bounds[i,,],1,function(x)diff(range(x))))
		c(bounds[i,nc,1], bounds[i,nc,2])
		}))
    # qui dobbiamo modificare la distribuzione... ci vuole una normale troncata
	if(is.null(controlCol)) controlCol <- do.call(cbind,lapply(1:nc,function(j)range(round(rnorm(N,colSums(muTab),sqrt(colSums(sdTab^2)))))))
	iter <- t <- 1
	while(t <= nIter){
		
		sim <- do.call(rbind,lapply(1:nr,function(i){
			
			rrow <- rep(NA,nc)
			repeat{
				### VARSTOCK structural 0
				if(muTab[i,nc]==0){
					nc <- max(which(muTab[i,-nc]!=0))
					rrow[(1:length(rrow))>nc]<-0
					
					maxTol <- which.max(apply(bounds[i,,],1,function(x)diff(range(x))))
					
					rrow[-c(maxTol,(nc+1):length(rrow))] <- round(rtnorm(nc-1, mean=unlist(muTab[i,-c(maxTol,(nc+1):length(rrow))]), sd=sdev, lower=bounds[i,-c(maxTol,(nc+1):length(rrow)),1],upper=bounds[i,-c(maxTol,(nc+1):length(rrow)),2]),0)

					rrow[maxTol] <- (n0[i]-sum(rrow[-c(maxTol,nc+1:length(rrow))]))
					nc<-maxTol
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

ptm <- proc.time()	
tab <- sampleTables(n0,muTab,bounds,controlCol=controlCol1)
proc.time() - ptm

tab


# Confronto 

compaTab = function(tab){
  eqv.tab = list()                     # Inizializzo il listone contenente gli indici delle tabelle equivalenti...
  eqv.idx = 0                          # ...e l'indice che le scorre...
  act.tab = 1:dim(tab[[1]])[3]        # Inizializzo indice tabelle attive
  
  while ( length(act.tab) > 0 ){                            # Finché ho tabelle attive...
    ref.tab = act.tab[1]                                    # ...seleziono la prima attiva come tabella di riferimento
    eqv.idx = eqv.idx + 1                                   # ...aggiungo un elemento alla lista delle tabelle equivalenti...  
    eqv.tab[[eqv.idx]] = ref.tab                            # ...lo inizializzo...
    names(eqv.tab)[eqv.idx] <- paste("Table", ref.tab)      # ...lo battezzo...
    act.tab = act.tab[-1]                                   # ...lo rimuovo dalla lista della tabelle attive
    act.idx = act.tab                                       # Salva gli indici per il ciclo <for> successivo... 
    for (j in act.idx){                                     # Confronto la prima ancora attiva con tutte le successive ancora attive...
      if ( all(tab[[1]][,,ref.tab] == tab[[1]][,,j]) ){     # ...se uguale...
        eqv.tab[[eqv.idx]] = c(eqv.tab[[eqv.idx]], j)       # ...aggiungo l'indice attuale al mucchio...
        act.tab = setdiff(act.tab, j)                       # ...e lo rimuovo da quelle attive...
      }
    }  
  }
  return(eqv.tab)
}

out <- compaTab(tab)
out

# Per avere le frequenze...
as.matrix(lapply(out, length))

Scenario <- read.csv("ValoriItalia.csv",sep=";",row.names=1)[,-(7:9)]

rmse <- do.call(rbind,lapply(1:dim(tab[[1]])[3],function(t)sqrt((sum((tab[[1]][,,t]-Scenario)^2))/(length(tab[[1]][,,t])))))
rmse
rrmse <- do.call(rbind,lapply(1:dim(tab[[1]])[3],function(t)sqrt((sum(((tab[[1]][,,t]-Scenario)/(Scenario+0.0000000001))^2))/(length(tab[[1]][,,t])))))
rrmse

order(rrmse)
min(rrmse)


# check
aa<-order(rrmse)[1]
aa
sqrt((sum(((tab[[1]][,,aa]-Scenario)/(Scenario+0.0000001))^2))/length(tab[[1]][,,aa]))
tab[[1]][,,aa]

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
	rmse[i] <- sqrt(sum(unlist(lapply(res[[1]],function(t)sum((t-Scenario)^2)))))
	reject[i] <- res[[2]]
	gc()
	}
plot(1/den,time,type="b",col="blue",main="Time vs. Strenght of the Constraint",xlab="narrowness")
points(1/den,time,col="red",pch=16)

