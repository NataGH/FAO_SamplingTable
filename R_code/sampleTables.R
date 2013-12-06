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



system.time(tab <- sampleTables(n0,muTab,bounds,controlCol=controlCol1))
tab


unlist(lapply(tab@tables,function(x)attr(x,"mult")))
length(tab@tables)



### now let's say, for instance, the 3rd, 56th and 101st rows of the previous best table were OK
system.time(tab2 <- sampleTables(n0,muTab,bounds,controlCol=controlCol1,fixed=c(3,56,101),fixedRows=tab@bestTab[c(3,56,101),]))

tab@bestTab[c(3,56,101),]
tab2@bestTab[c(3,56,101),]


### now let's fix the ()-th element of the table
system.time(tab2 <- sampleTables(n0,muTab,bounds,controlCol=controlCol1,fixed=c(3,56,101),fixedRows=tab@bestTab[c(3,56,101),]))

tab@bestTab[c(3,56,101),]
tab2@bestTab[c(3,56,101),]




### now let's use the MSE to pick the best table
muTab <- read.csv("ExpectedValues_t_Italy_2006.csv",row.names=1,sep=",")[-114,-(8:10)]
n0<-read.csv("ExpectedValues_t_Italy_2006.csv",header=TRUE,row.names=1,sep=",")[-114,8]

bounds <- array(NA,c(dim(muTab),2))
dimnames(bounds) <- c(dimnames(muTab),list(c("Lower","Upper")))
bounds[,,1] <- as.matrix(read.csv("LowerBounds_t_Italy_2006.csv",row.names=1,sep=','))[-114,colnames(muTab)]
bounds[,,2] <- as.matrix(read.csv("UpperBounds_t_Italy_2006.csv",row.names=1,sep=','))[-114,colnames(muTab)]


lowCol<-read.csv("LowerBounds_t_Italy_2006.csv",row.names=1,sep=',')[114,]
uppCol<-read.csv("UpperBounds_t_Italy_2006.csv",row.names=1,sep=',')[114,]
controlCol1 <- rbind(lowCol,uppCol)
controlCol1[,ncol(controlCol1)]<- controlCol1[,ncol(controlCol1)]*c(1/5,5)
#controlCol1[1,] <- -controlCol1[2,]


Scenario <- read.csv("Scenario_t_Italy_2006.csv",sep=",",nrows=nrow(muTab))
rownames(Scenario) <- make.unique(as.character(Scenario[,1]))
Scenario <- Scenario[,-c(1,9:11)]

system.time(tab3 <- sampleTables(n0,muTab,bounds,controlCol=controlCol1,objFun=rrmseObj(Scenario)))
# system.time(tab3 <- sampleTables(n0,muTab,bounds,controlCol=controlCol1,objFun=rrmseObj(Scenario))



### Now let's consider the case one row is all zero
muTab1 <- muTab
muTab1[4,] <- muTab1[4,]*0
system.time(tab4 <- sampleTables(n0,muTab1,bounds,controlCol=controlCol1))



### do we still care about this?
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



