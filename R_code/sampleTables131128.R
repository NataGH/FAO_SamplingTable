rm(list=ls())
gc()


source("FAOTables.R")



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

system.time(tab <- sampleTables(n0,muTab,bounds,controlCol=controlCol1))
tab


# now let's say, for instance, the 3rd, 56th and 101st rows of the previous best table were OK
system.time(tab2 <- sampleTables(n0,muTab,bounds,controlCol=controlCol1,fixed=c(3,56,101),fixedRows=tab@bestTab[c(3,56,101),]))
tab2

tab@bestTab[c(3,56,101),]
tab2@bestTab[c(3,56,101),]

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



#### Let's use the RRMSE as a loss function


system.time(tab <- sampleTables(n0,muTab,bounds,controlCol=controlCol1))