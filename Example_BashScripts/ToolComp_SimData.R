library(ggplot2)


#### ---- Calculate Sensitivity & Specificity ---- 

toolNames = c("DiTector", "ViReMa")

# load in data
sumDataMBP <- read.table("~/Documents/Project/Summary_res_mbp.txt", sep = "\t", header = FALSE)
colnames(sumDataMBP) <- c("row_names", "count")
row.names(sumDataMBP) <- sumDataMBP$row_names
sumDataMBP <- subset(sumDataMBP, select = -row_names)

sumDataID <- read.table("~/Documents/Project/Summary_res_id.txt", sep = "\t", header = FALSE)
colnames(sumDataID) <- c("row_names", "count")
row.names(sumDataID) <- sumDataID$row_names
sumDataID <- subset(sumDataID, select = -row_names)

sumDataCB <- read.table("~/Documents/Project/Summary_res_cb.txt", sep = "\t", header = FALSE)
colnames(sumDataCB) <- c("row_names", "count")
row.names(sumDataCB) <- sumDataCB$row_names
sumDataCB <- subset(sumDataCB, select = -row_names)

sumDataMS <- read.table("~/Documents/Project/Summary_res_ms.txt", sep = "\t", header = FALSE)
colnames(sumDataMS) <- c("row_names", "count")
row.names(sumDataMS) <- sumDataMS$row_names
sumDataMS <- subset(sumDataMS, select = -row_names)

sumDataND <- read.table("~/Documents/Project/Summary_res_nd.txt", sep = "\t", header = FALSE)
colnames(sumDataND) <- c("row_names", "count")
row.names(sumDataND) <- sumDataND$row_names
sumDataND <- subset(sumDataND, select = -row_names)

tndit <- (sumDataND[1,]-sumDataND[2,])
tnvir <- (sumDataND[1,]-sumDataND[3,])


dataList <- list(sumDataMBP, sumDataID, sumDataCB, sumDataMS)
methodList <- list("MBP", "INDEL", "CopyBack", "MultiSeg")

calcList <- list()
x = 1

for(data in dataList){
  
  tpdit <- data[4,]
  tpvir <- data[5,]
  fndit <- data[1,] - data[2,]
  fnvir <- data[1,] - data[3,]
  fpdit <- data[6,]
  fpvir <- data[7,]
  
  st_dit <- tpdit/(tpdit + fndit) * 100
  sp_dit <- tndit/(fpdit + tndit) * 100
  st_vir <- tpvir/(tpvir + fnvir) * 100
  sp_vir <- tnvir/(fpvir + tnvir) * 100
  
  sensitivity <- c(st_dit, st_vir)
  specificity <- c(sp_dit, sp_vir)
  
  ppv_dit <- tpdit/(tpdit + tndit) * 100
  npv_dit <- tndit/(fndit + tndit) * 100
  ppv_vir <- tpvir/(tpvir + tnvir) * 100
  npv_vir <- tnvir/(fnvir + tnvir) * 100  
  
  ppv <- c(ppv_dit, ppv_vir)
  npv <- c(npv_dit, npv_vir)
  
  acc_dit <- (tpdit + tndit) / (tpdit + tndit + fpdit + fndit)
  acc_vir <- (tpvir + tnvir) / (tpvir + tnvir + fpvir + fnvir)
  
  acc <- c(acc_dit, acc_vir)
  
  method <- c(methodList[x], methodList[x])
  
  calcList[[x]] <- data.frame(sensitivity, specificity, ppv, npv, toolNames, method, acc)
  colnames(calcList[[x]]) <- c("Sensitivity", "Specificity", "PPV", "NPV", "Tool", "Method", "Accuracy")
  x = x + 1
}

# Use GGplot, with sensitivity on the x axis and specificity on the y axis 

plot_strict <- ggplot(calcList[[1]], aes(x=Sensitivity, y=Specificity, shape=Tool, colour=Method), size = 4) + 
  geom_point(aes(x=calcList[[1]]$Sensitivity, y=calcList[[1]]$Specificity, shape=calcList[[1]]$Tool, colour=calcList[[1]]$Method), size = 4) +
  geom_point(aes(x=calcList[[2]]$Sensitivity, y=calcList[[2]]$Specificity, shape=calcList[[2]]$Tool, colour=calcList[[2]]$Method), size = 4) +
  geom_point(aes(x=calcList[[3]]$Sensitivity, y=calcList[[3]]$Specificity, shape=calcList[[3]]$Tool, colour=calcList[[3]]$Method), size = 4) +
  geom_point(aes(x=calcList[[4]]$Sensitivity, y=calcList[[4]]$Specificity, shape=calcList[[4]]$Tool, colour=calcList[[4]]$Method), size = 4) +
  xlim(0,100) + ylim(0,100)

pplot_strict <- ggplot(calcList[[1]], aes(x=PPV, y=NPV, shape=Tool, colour=Method), size = 3) + 
  geom_point(aes(x=calcList[[1]]$PPV, y=calcList[[1]]$NPV, shape=calcList[[1]]$Tool, colour=calcList[[1]]$Method), size = 3) +
  geom_point(aes(x=calcList[[2]]$PPV, y=calcList[[2]]$NPV, shape=calcList[[2]]$Tool, colour=calcList[[2]]$Method), size = 3) +
  geom_point(aes(x=calcList[[3]]$PPV, y=calcList[[3]]$NPV, shape=calcList[[3]]$Tool, colour=calcList[[3]]$Method), size = 3) +
  geom_point(aes(x=calcList[[4]]$PPV, y=calcList[[4]]$NPV, shape=calcList[[4]]$Tool, colour=calcList[[4]]$Method), size = 3) 


#### ---- With standard deviation plot ---- 

# load in data
sumDataMBP <- read.table("~/Documents/Project/Summary_resSD_mbp.txt", sep = "\t", header = FALSE)
colnames(sumDataMBP) <- c("row_names", "count")
row.names(sumDataMBP) <- sumDataMBP$row_names
sumDataMBP <- subset(sumDataMBP, select = -row_names)

sumDataID <- read.table("~/Documents/Project/Summary_resSD_id.txt", sep = "\t", header = FALSE)
colnames(sumDataID) <- c("row_names", "count")
row.names(sumDataID) <- sumDataID$row_names
sumDataID <- subset(sumDataID, select = -row_names)

sumDataCB <- read.table("~/Documents/Project/Summary_resSD_cb.txt", sep = "\t", header = FALSE)
colnames(sumDataCB) <- c("row_names", "count")
row.names(sumDataCB) <- sumDataCB$row_names
sumDataCB <- subset(sumDataCB, select = -row_names)

sumDataMS <- read.table("~/Documents/Project/Summary_resSD_ms.txt", sep = "\t", header = FALSE)
colnames(sumDataMS) <- c("row_names", "count")
row.names(sumDataMS) <- sumDataMS$row_names
sumDataMS <- subset(sumDataMS, select = -row_names)

sumDataND <- read.table("~/Documents/Project/Summary_resSD_nd.txt", sep = "\t", header = FALSE)
colnames(sumDataND) <- c("row_names", "count")
row.names(sumDataND) <- sumDataND$row_names
sumDataND <- subset(sumDataND, select = -row_names)

tndit <- (sumDataND[1,]-sumDataND[2,])
tnvir <- (sumDataND[1,]-sumDataND[3,])

sddataList <- list(sumDataMBP, sumDataID, sumDataCB, sumDataMS)

sdcalcList <- list()
x = 1

for(data in sddataList){
  tpdit <- data[4,]
  tpvir <- data[5,]
  fndit <- data[1,] - data[2,]
  fnvir <- data[1,] - data[3,]
  fpdit <- data[6,]
  fpvir <- data[7,]
  
  st_dit <- tpdit/(tpdit + fndit) * 100
  sp_dit <- tndit/(fpdit + tndit) * 100
  st_vir <- tpvir/(tpvir + fnvir) * 100
  sp_vir <- tnvir/(fpvir + tnvir) * 100
  
  sensitivity <- c(st_dit, st_vir)
  specificity <- c(sp_dit, sp_vir)
  
  ppv_dit <- tpdit/(tpdit + tndit) * 100
  npv_dit <- tndit/(fndit + tndit) * 100
  ppv_vir <- tpvir/(tpvir + tnvir) * 100
  npv_vir <- tnvir/(fnvir + tnvir) * 100  
  
  ppv <- c(ppv_dit, ppv_vir)
  npv <- c(npv_dit, npv_vir)
  
  acc_dit <- (tpdit + tndit) / (tpdit + tndit + fpdit + fndit)
  acc_vir <- (tpvir + tnvir) / (tpvir + tnvir + fpvir + fnvir)
  
  acc <- c(acc_dit, acc_vir)
  
  method <- c(methodList[x], methodList[x])
  
  sdcalcList[[x]] <- data.frame(sensitivity, specificity, ppv, npv, toolNames, method, acc)
  colnames(sdcalcList[[x]]) <- c("Sensitivity", "Specificity", "PPV", "NPV", "Tool", "Method", "Accuracy")
  x = x + 1
}

# Use GGplot, with sensitivity on the x axis and specificity on the y axis 

plot_sd <- ggplot(sdcalcList[[1]], aes(x=Sensitivity, y=Specificity, shape=Tool, colour=Method), size = 4) + 
  geom_point(aes(x=sdcalcList[[1]]$Sensitivity, y=sdcalcList[[1]]$Specificity, shape=sdcalcList[[1]]$Tool, colour=sdcalcList[[1]]$Method), size = 4) +
  geom_point(aes(x=sdcalcList[[2]]$Sensitivity, y=sdcalcList[[2]]$Specificity, shape=sdcalcList[[2]]$Tool, colour=sdcalcList[[2]]$Method), size = 4) +
  geom_point(aes(x=sdcalcList[[3]]$Sensitivity, y=sdcalcList[[3]]$Specificity, shape=sdcalcList[[3]]$Tool, colour=sdcalcList[[3]]$Method), size = 4) +
  geom_point(aes(x=sdcalcList[[4]]$Sensitivity, y=sdcalcList[[4]]$Specificity, shape=sdcalcList[[4]]$Tool, colour=sdcalcList[[4]]$Method), size = 4) +
  xlim(0,100) + ylim(0,100)

pplot_sd <- ggplot(sdcalcList[[1]], aes(x=PPV, y=NPV, shape=Tool, colour=Method), size = 3) + 
  geom_point(aes(x=sdcalcList[[1]]$PPV, y=sdcalcList[[1]]$NPV, shape=sdcalcList[[1]]$Tool, colour=sdcalcList[[1]]$Method), size = 3) +
  geom_point(aes(x=sdcalcList[[2]]$PPV, y=sdcalcList[[2]]$NPV, shape=sdcalcList[[2]]$Tool, colour=sdcalcList[[2]]$Method), size = 3) +
  geom_point(aes(x=sdcalcList[[3]]$PPV, y=sdcalcList[[3]]$NPV, shape=sdcalcList[[3]]$Tool, colour=sdcalcList[[3]]$Method), size = 3) +
  geom_point(aes(x=sdcalcList[[4]]$PPV, y=sdcalcList[[4]]$NPV, shape=sdcalcList[[4]]$Tool, colour=sdcalcList[[4]]$Method), size = 3) 

