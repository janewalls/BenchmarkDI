
library(ggplot2)
library(dplyr)
library(BiocManager)
library(limma)
library(VennDiagram)

# Key: dit = DI-Tector, vir/vrm = ViReMa


#### ---- Calculate Sensitivity & Specificity ---- 

toolNames = c("DiTector", "ViReMa")

# load in data - summary_results.csv file from OutParse.py
sumDataMBP <- read.table("Summary_res_mbp.txt", sep = "\t", header = FALSE)
colnames(sumDataMBP) <- c("row_names", "count")
row.names(sumDataMBP) <- sumDataMBP$row_names
sumDataMBP <- subset(sumDataMBP, select = -row_names)

sumDataID <- read.table("Summary_res_id.txt", sep = "\t", header = FALSE)
colnames(sumDataID) <- c("row_names", "count")
row.names(sumDataID) <- sumDataID$row_names
sumDataID <- subset(sumDataID, select = -row_names)

sumDataCB <- read.table("Summary_res_cb.txt", sep = "\t", header = FALSE)
colnames(sumDataCB) <- c("row_names", "count")
row.names(sumDataCB) <- sumDataCB$row_names
sumDataCB <- subset(sumDataCB, select = -row_names)

sumDataMS <- read.table("Summary_res_ms.txt", sep = "\t", header = FALSE)
colnames(sumDataMS) <- c("row_names", "count")
row.names(sumDataMS) <- sumDataMS$row_names
sumDataMS <- subset(sumDataMS, select = -row_names)

sumDataND <- read.table("Summary_res_nd.txt", sep = "\t", header = FALSE)
colnames(sumDataND) <- c("row_names", "count")
row.names(sumDataND) <- sumDataND$row_names
sumDataND <- subset(sumDataND, select = -row_names)

# Calculate True Negative
tndit <- (sumDataND[1,]-sumDataND[2,])
tnvir <- (sumDataND[1,]-sumDataND[3,])


dataList <- list(sumDataMBP, sumDataID, sumDataCB, sumDataMS)
methodList <- list("MBP", "INDEL", "CopyBack", "MultiSeg")

calcList <- list()
x = 1 # List position
# tp = True Positive, fn = False Negative, fp = False Positive, ppv = Positive Predicited Value, npv = Negative Predictive Value
# Do calculations for all samples
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
  
  acc_dit <- (tpdit + tndit) / (tpdit + tndit + fpdit + fndit) * 100
  acc_vir <- (tpvir + tnvir) / (tpvir + tnvir + fpvir + fnvir) * 100
  
  acc <- c(acc_dit, acc_vir) # Accuracy
  
  method <- c(methodList[[x]], methodList[[x]])
  
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

# Use GGplot, with ppv on the x axis and npv on the y axis 
pplot_strict <- ggplot(calcList[[1]], aes(x=PPV, y=NPV, shape=Tool, colour=Method), size = 3) + 
  geom_point(aes(x=calcList[[1]]$PPV, y=calcList[[1]]$NPV, shape=calcList[[1]]$Tool, colour=calcList[[1]]$Method), size = 3) +
  geom_point(aes(x=calcList[[2]]$PPV, y=calcList[[2]]$NPV, shape=calcList[[2]]$Tool, colour=calcList[[2]]$Method), size = 3) +
  geom_point(aes(x=calcList[[3]]$PPV, y=calcList[[3]]$NPV, shape=calcList[[3]]$Tool, colour=calcList[[3]]$Method), size = 3) +
  geom_point(aes(x=calcList[[4]]$PPV, y=calcList[[4]]$NPV, shape=calcList[[4]]$Tool, colour=calcList[[4]]$Method), size = 3) 


#### ---- Same With standard deviation plot ---- 

# load in data - summary_results.csv file from OutParse.py
sumDataMBP <- read.table("Summary_resSD_mbp.txt", sep = "\t", header = FALSE)
colnames(sumDataMBP) <- c("row_names", "count")
row.names(sumDataMBP) <- sumDataMBP$row_names
sumDataMBP <- subset(sumDataMBP, select = -row_names)

sumDataID <- read.table("Summary_resSD_id.txt", sep = "\t", header = FALSE)
colnames(sumDataID) <- c("row_names", "count")
row.names(sumDataID) <- sumDataID$row_names
sumDataID <- subset(sumDataID, select = -row_names)

sumDataCB <- read.table("Summary_resSD_cb.txt", sep = "\t", header = FALSE)
colnames(sumDataCB) <- c("row_names", "count")
row.names(sumDataCB) <- sumDataCB$row_names
sumDataCB <- subset(sumDataCB, select = -row_names)

sumDataMS <- read.table("Summary_resSD_ms.txt", sep = "\t", header = FALSE)
colnames(sumDataMS) <- c("row_names", "count")
row.names(sumDataMS) <- sumDataMS$row_names
sumDataMS <- subset(sumDataMS, select = -row_names)

sumDataND <- read.table("Summary_resSD_nd.txt", sep = "\t", header = FALSE)
colnames(sumDataND) <- c("row_names", "count")
row.names(sumDataND) <- sumDataND$row_names
sumDataND <- subset(sumDataND, select = -row_names)

# Calculate True Negative for each tool
tndit <- (sumDataND[1,]-sumDataND[2,])
tnvir <- (sumDataND[1,]-sumDataND[3,])

sddataList <- list(sumDataMBP, sumDataID, sumDataCB, sumDataMS)

sdcalcList <- list()
x = 1 # List position

# tp = True Positive, fn = False Negative, fp = False Positive, ppv = Positive Predicited Value, npv = Negative Predictive Value
# Do calculations for all samples
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
  
  acc_dit <- (tpdit + tndit) / (tpdit + tndit + fpdit + fndit) * 100
  acc_vir <- (tpvir + tnvir) / (tpvir + tnvir + fpvir + fnvir) * 100
  
  acc <- c(acc_dit, acc_vir)
  
  method <- c(methodList[[x]], methodList[[x]])
  
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

# Use GGplot, with ppv on the x axis and npv on the y axis 
pplot_sd <- ggplot(sdcalcList[[1]], aes(x=PPV, y=NPV, shape=Tool, colour=Method), size = 3) + 
  geom_point(aes(x=sdcalcList[[1]]$PPV, y=sdcalcList[[1]]$NPV, shape=sdcalcList[[1]]$Tool, colour=sdcalcList[[1]]$Method), size = 3) +
  geom_point(aes(x=sdcalcList[[2]]$PPV, y=sdcalcList[[2]]$NPV, shape=sdcalcList[[2]]$Tool, colour=sdcalcList[[2]]$Method), size = 3) +
  geom_point(aes(x=sdcalcList[[3]]$PPV, y=sdcalcList[[3]]$NPV, shape=sdcalcList[[3]]$Tool, colour=sdcalcList[[3]]$Method), size = 3) +
  geom_point(aes(x=sdcalcList[[4]]$PPV, y=sdcalcList[[4]]$NPV, shape=sdcalcList[[4]]$Tool, colour=sdcalcList[[4]]$Method), size = 3) 




#### ---- Bar Chart for Segment distribution ---- 

segNames <- c("JX680447seg1","JX680448seg2", "JX680449seg3", "JX680450seg4", "JX680451seg5", "JX680452seg6", "JX680453seg7", "JX680454seg8", "JX680455seg9", "JX680456seg10", "Multi Segment", "JX680447seg1","JX680448seg2", "JX680449seg3", "JX680450seg4", "JX680451seg5", "JX680452seg6", "JX680453seg7", "JX680454seg8", "JX680455seg9", "JX680456seg10", "Multi Segment")
toolNames <- c("DITector", "DITector", "DITector", "DITector", "DITector", "DITector", "DITector", "DITector", "DITector", "DITector", "DITector", "ViReMa", "ViReMa", "ViReMa", "ViReMa", "ViReMa", "ViReMa", "ViReMa", "ViReMa", "ViReMa", "ViReMa", "ViReMa")
cols = c("darkblue", "cornflowerblue", "darkcyan", "aquamarine1", "darkgreen", "chartreuse", "darkolivegreen3", "chocolate3", "darkorange", "darkgoldenrod1", "darkorchid1")

# load in Bld_33_S1 - parser_output.csv from OutParse.py
data_33 <- read.table("Bld_33_S1_SDoutput.csv", sep = "\t", header = TRUE)

# Creates data frames for DIPs across single segments
seg1 <- subset(data_33, BP_seg == 1 & RI_seg == 1)
seg2 <- subset(data_33, BP_seg == 2 & RI_seg == 2)
seg3 <- subset(data_33, BP_seg == 3 & RI_seg == 3)
seg4 <- subset(data_33, BP_seg == 4 & RI_seg == 4)
seg5 <- subset(data_33, BP_seg == 5 & RI_seg == 5)
seg6 <- subset(data_33, BP_seg == 6 & RI_seg == 6)
seg7 <- subset(data_33, BP_seg == 7 & RI_seg == 7)
seg8 <- subset(data_33, BP_seg == 8 & RI_seg == 8)
seg9 <- subset(data_33, BP_seg == 9 & RI_seg == 9)
seg10 <- subset(data_33, BP_seg == 10 & RI_seg == 10)
multiSeg <- subset(data_33, BP_seg != RI_seg) # If Break point and Re-initiation points occur on different segments

# Create bar chart of DIP counts across individual segments or multiple segments
dat_comb <- c(length(which(seg1[5]>0)), length(which(seg2[5]>0)), length(which(seg3[5]>0)), length(which(seg4[5]>0)), length(which(seg5[5]>0)), length(which(seg6[5]>0)), length(which(seg7[5]>0)), length(which(seg8[5]>0)), length(which(seg9[5]>0)), length(which(seg10[5]>0)), length(which(multiSeg[5]>0)), length(which(seg1[6]>0)), length(which(seg2[6]>0)), length(which(seg3[6]>0)), length(which(seg4[6]>0)), length(which(seg5[6]>0)), length(which(seg6[6]>0)), length(which(seg7[6]>0)), length(which(seg8[6]>0)), length(which(seg9[6]>0)), length(which(seg10[6]>0)), length(which(multiSeg[6]>0)))
barFrame <- data.frame(dat_comb, segNames, toolNames)
colnames(barFrame) <- c("DIP_Count", "Segment", "Tool")
barPlot <- ggplot(barFrame, aes(x=Segment, y=DIP_Count)) + geom_bar(aes(fill=Tool), width=0.7, position=position_dodge(width=0.75), stat="identity") + theme(axis.text.x = element_text(angle = 75, vjust=0.7)) + ylab("Identified DIPs") + scale_y_log10() 

# Find proportion of DIPs found in relative segments
ditDIP <- length(which(data_33[5]>0))
ditVRM <- length(which(data_33[6]>0))

# Create stacked bar chart to illistrate proportions
barFrame$Total <- c(ditDIP, ditDIP, ditDIP, ditDIP, ditDIP, ditDIP, ditDIP, ditDIP, ditDIP, ditDIP, ditDIP, ditVRM, ditVRM, ditVRM, ditVRM, ditVRM, ditVRM, ditVRM, ditVRM, ditVRM, ditVRM, ditVRM)
barFrame$Proportion <- ((barFrame[1]/barFrame[4]) *100)
barFrame$Proportion <- unlist(barFrame$Proportion)
stackPlot <- ggplot(barFrame, aes(x=Tool, y=Proportion, fill=Segment)) + geom_bar(aes(fill=Segment), position='stack', stat="identity") + theme(axis.text.x = element_text(vjust=0.7)) + ylab("Proportion of Identified DIPs %") + ylim(0,100) + scale_fill_manual(values=cols)




# load in Bld_36_S2 - parser_output.csv from OutParse.py
data_36 <- read.table("Bld_36_S2_SDoutput.csv", sep = "\t", header = TRUE)

# Creates data frames for DIPs across single segments
seg1 <- subset(data_36, BP_seg == 1 & RI_seg == 1)
seg2 <- subset(data_36, BP_seg == 2 & RI_seg == 2)
seg3 <- subset(data_36, BP_seg == 3 & RI_seg == 3)
seg4 <- subset(data_36, BP_seg == 4 & RI_seg == 4)
seg5 <- subset(data_36, BP_seg == 5 & RI_seg == 5)
seg6 <- subset(data_36, BP_seg == 6 & RI_seg == 6)
seg7 <- subset(data_36, BP_seg == 7 & RI_seg == 7)
seg8 <- subset(data_36, BP_seg == 8 & RI_seg == 8)
seg9 <- subset(data_36, BP_seg == 9 & RI_seg == 9)
seg10 <- subset(data_36, BP_seg == 10 & RI_seg == 10)
multiSeg <- subset(data_36, BP_seg != RI_seg) # If Break point and Re-initiation points occur on different segments

# Create bar chart of DIP counts across individual segments or multiple segments
dat_comb <- c(length(which(seg1[5]>0)), length(which(seg2[5]>0)), length(which(seg3[5]>0)), length(which(seg4[5]>0)), length(which(seg5[5]>0)), length(which(seg6[5]>0)), length(which(seg7[5]>0)), length(which(seg8[5]>0)), length(which(seg9[5]>0)), length(which(seg10[5]>0)), length(which(multiSeg[5]>0)), length(which(seg1[6]>0)), length(which(seg2[6]>0)), length(which(seg3[6]>0)), length(which(seg4[6]>0)), length(which(seg5[6]>0)), length(which(seg6[6]>0)), length(which(seg7[6]>0)), length(which(seg8[6]>0)), length(which(seg9[6]>0)), length(which(seg10[6]>0)), length(which(multiSeg[6]>0)))
barFrame <- data.frame(dat_comb, segNames, toolNames)
colnames(barFrame) <- c("DIP_Count", "Segment", "Tool")
barPlot <- ggplot(barFrame, aes(x=Segment, y=DIP_Count)) + geom_bar(aes(fill=Tool), width=0.7, position=position_dodge(width=0.75), stat="identity") + theme(axis.text.x = element_text(angle = 75, vjust=0.7)) + ylab("Identified DIPs")  + scale_y_log10() 

# Find proportion of DIPs found in relative segments
ditDIP <- length(which(data_36[5]>0))
ditVRM <- length(which(data_36[6]>0))

# Create stacked bar chart to illistrate proportions
barFrame$Total <- c(ditDIP, ditDIP, ditDIP, ditDIP, ditDIP, ditDIP, ditDIP, ditDIP, ditDIP, ditDIP, ditDIP, ditVRM, ditVRM, ditVRM, ditVRM, ditVRM, ditVRM, ditVRM, ditVRM, ditVRM, ditVRM, ditVRM)
barFrame$Proportion <- ((barFrame[1]/barFrame[4]) *100)
barFrame$Proportion <- unlist(barFrame$Proportion)
stackPlot <- ggplot(barFrame, aes(x=Tool, y=Proportion, fill=Segment)) + geom_bar(aes(fill=Segment), position='stack', stat="identity") + theme(axis.text.x = element_text(vjust=0.7)) + ylab("Proportion of Identified DIPs %") + ylim(0,100) + scale_fill_manual(values=cols)




# load in Bld_38_S4 - parser_output.csv from OutParse.py
data_38 <- read.table("Bld_38_S4_SDoutput.csv", sep = "\t", header = TRUE)

# Creates data frames for DIPs across single segments
seg1 <- subset(data_38, BP_seg == 1 & RI_seg == 1)
seg2 <- subset(data_38, BP_seg == 2 & RI_seg == 2)
seg3 <- subset(data_38, BP_seg == 3 & RI_seg == 3)
seg4 <- subset(data_38, BP_seg == 4 & RI_seg == 4)
seg5 <- subset(data_38, BP_seg == 5 & RI_seg == 5)
seg6 <- subset(data_38, BP_seg == 6 & RI_seg == 6)
seg7 <- subset(data_38, BP_seg == 7 & RI_seg == 7)
seg8 <- subset(data_38, BP_seg == 8 & RI_seg == 8)
seg9 <- subset(data_38, BP_seg == 9 & RI_seg == 9)
seg10 <- subset(data_38, BP_seg == 10 & RI_seg == 10)
multiSeg <- subset(data_38, BP_seg != RI_seg) # If Break point and Re-initiation points occur on different segments

# Create bar chart of DIP counts across individual segments or multiple segments
dat_comb <- c(length(which(seg1[5]>0)), length(which(seg2[5]>0)), length(which(seg3[5]>0)), length(which(seg4[5]>0)), length(which(seg5[5]>0)), length(which(seg6[5]>0)), length(which(seg7[5]>0)), length(which(seg8[5]>0)), length(which(seg9[5]>0)), length(which(seg10[5]>0)), length(which(multiSeg[5]>0)), length(which(seg1[6]>0)), length(which(seg2[6]>0)), length(which(seg3[6]>0)), length(which(seg4[6]>0)), length(which(seg5[6]>0)), length(which(seg6[6]>0)), length(which(seg7[6]>0)), length(which(seg8[6]>0)), length(which(seg9[6]>0)), length(which(seg10[6]>0)), length(which(multiSeg[6]>0)))
barFrame <- data.frame(dat_comb, segNames, toolNames)
colnames(barFrame) <- c("DIP_Count", "Segment", "Tool")
barPlot <- ggplot(barFrame, aes(x=Segment, y=DIP_Count)) + geom_bar(aes(fill=Tool), width=0.7, position=position_dodge(width=0.75), stat="identity") + theme(axis.text.x = element_text(angle = 75, vjust=0.7)) + ylab("Identified DIPs") + scale_y_log10() 

# Find proportion of DIPs found in relative segments
ditDIP <- length(which(data_38[5]>0))
ditVRM <- length(which(data_38[6]>0))

# Create stacked bar chart to illistrate proportions
barFrame$Total <- c(ditDIP, ditDIP, ditDIP, ditDIP, ditDIP, ditDIP, ditDIP, ditDIP, ditDIP, ditDIP, ditDIP, ditVRM, ditVRM, ditVRM, ditVRM, ditVRM, ditVRM, ditVRM, ditVRM, ditVRM, ditVRM, ditVRM)
barFrame$Proportion <- ((barFrame[1]/barFrame[4]) *100)
barFrame$Proportion <- unlist(barFrame$Proportion)
stackPlot <- ggplot(barFrame, aes(x=Tool, y=Proportion, fill=Segment)) + geom_bar(aes(fill=Segment), position='stack', stat="identity") + theme(axis.text.x = element_text(vjust=0.7)) + ylab("Proportion of Identified DIPs %") + ylim(0,100) + scale_fill_manual(values=cols)




# load in Bld_40_S6 - parser_output.csv from OutParse.py
data_40 <- read.table("Bld_40_S6_SDoutput.csv", sep = "\t", header = TRUE)

# Creates data frames for DIPs across single segments
seg1 <- subset(data_40, BP_seg == 1 & RI_seg == 1)
seg2 <- subset(data_40, BP_seg == 2 & RI_seg == 2)
seg3 <- subset(data_40, BP_seg == 3 & RI_seg == 3)
seg4 <- subset(data_40, BP_seg == 4 & RI_seg == 4)
seg5 <- subset(data_40, BP_seg == 5 & RI_seg == 5)
seg6 <- subset(data_40, BP_seg == 6 & RI_seg == 6)
seg7 <- subset(data_40, BP_seg == 7 & RI_seg == 7)
seg8 <- subset(data_40, BP_seg == 8 & RI_seg == 8)
seg9 <- subset(data_40, BP_seg == 9 & RI_seg == 9)
seg10 <- subset(data_40, BP_seg == 10 & RI_seg == 10)
multiSeg <- subset(data_40, BP_seg != RI_seg) # If Break point and Re-initiation points occur on different segments

# Create bar chart of DIP counts across individual segments or multiple segments
dat_comb <- c(length(which(seg1[5]>0)), length(which(seg2[5]>0)), length(which(seg3[5]>0)), length(which(seg4[5]>0)), length(which(seg5[5]>0)), length(which(seg6[5]>0)), length(which(seg7[5]>0)), length(which(seg8[5]>0)), length(which(seg9[5]>0)), length(which(seg10[5]>0)), length(which(multiSeg[5]>0)), length(which(seg1[6]>0)), length(which(seg2[6]>0)), length(which(seg3[6]>0)), length(which(seg4[6]>0)), length(which(seg5[6]>0)), length(which(seg6[6]>0)), length(which(seg7[6]>0)), length(which(seg8[6]>0)), length(which(seg9[6]>0)), length(which(seg10[6]>0)), length(which(multiSeg[6]>0)))
barFrame <- data.frame(dat_comb, segNames, toolNames)
colnames(barFrame) <- c("DIP_Count", "Segment", "Tool")
barPlot <- ggplot(barFrame, aes(x=Segment, y=DIP_Count)) + geom_bar(aes(fill=Tool), width=0.7, position=position_dodge(width=0.75), stat="identity") + theme(axis.text.x = element_text(angle = 75, vjust=0.7)) + ylab("Identified DIPs") + scale_y_log10()

# Find proportion of DIPs found in relative segments
ditDIP <- length(which(data_40[5]>0))
ditVRM <- length(which(data_40[6]>0))

# Create stacked bar chart to illistrate proportions
barFrame$Total <- c(ditDIP, ditDIP, ditDIP, ditDIP, ditDIP, ditDIP, ditDIP, ditDIP, ditDIP, ditDIP, ditDIP, ditVRM, ditVRM, ditVRM, ditVRM, ditVRM, ditVRM, ditVRM, ditVRM, ditVRM, ditVRM, ditVRM)
barFrame$Proportion <- ((barFrame[1]/barFrame[4]) *100)
barFrame$Proportion <- unlist(barFrame$Proportion)
stackPlot <- ggplot(barFrame, aes(x=Tool, y=Proportion, fill=Segment)) + geom_bar(aes(fill=Segment), position='stack', stat="identity") + theme(axis.text.x = element_text(vjust=0.7)) + ylab("Proportion of Identified DIPs %") + ylim(0,100) + scale_fill_manual(values=cols)




# load in Isol - parser_output.csv from OutParse.py
data_isol <- read.table("Isol_output/parser_output_SD.csv", sep = "\t", header = TRUE)

# Creates data frames for DIPs across single segments
seg1 <- subset(data_isol, BP_seg == 1 & RI_seg == 1)
seg2 <- subset(data_isol, BP_seg == 2 & RI_seg == 2)
seg3 <- subset(data_isol, BP_seg == 3 & RI_seg == 3)
seg4 <- subset(data_isol, BP_seg == 4 & RI_seg == 4)
seg5 <- subset(data_isol, BP_seg == 5 & RI_seg == 5)
seg6 <- subset(data_isol, BP_seg == 6 & RI_seg == 6)
seg7 <- subset(data_isol, BP_seg == 7 & RI_seg == 7)
seg8 <- subset(data_isol, BP_seg == 8 & RI_seg == 8)
seg9 <- subset(data_isol, BP_seg == 9 & RI_seg == 9)
seg10 <- subset(data_isol, BP_seg == 10 & RI_seg == 10)
multiSeg <- subset(data_isol, BP_seg != RI_seg) # If Break point and Re-initiation points occur on different segments

# Create bar chart of DIP counts across individual segments or multiple segments
dat_comb <- c(length(which(seg1[5]>0)), length(which(seg2[5]>0)), length(which(seg3[5]>0)), length(which(seg4[5]>0)), length(which(seg5[5]>0)), length(which(seg6[5]>0)), length(which(seg7[5]>0)), length(which(seg8[5]>0)), length(which(seg9[5]>0)), length(which(seg10[5]>0)), length(which(multiSeg[5]>0)), length(which(seg1[6]>0)), length(which(seg2[6]>0)), length(which(seg3[6]>0)), length(which(seg4[6]>0)), length(which(seg5[6]>0)), length(which(seg6[6]>0)), length(which(seg7[6]>0)), length(which(seg8[6]>0)), length(which(seg9[6]>0)), length(which(seg10[6]>0)), length(which(multiSeg[6]>0)))
barFrame <- data.frame(dat_comb, segNames, toolNames)
colnames(barFrame) <- c("DIP_Count", "Segment", "Tool")
barPlot <- ggplot(barFrame, aes(x=Segment, y=DIP_Count)) + geom_bar(aes(fill=Tool), width=0.7, position=position_dodge(width=0.75), stat="identity") + theme(axis.text.x = element_text(angle = 75, vjust=0.7)) + ylab("Identified DIPs") + scale_y_log10() 


# Find proportion of DIPs found in relative segments
ditDIP <- length(which(data_isol[5]>0))
ditVRM <- length(which(data_isol[6]>0))

# Create stacked bar chart to illistrate proportions
barFrame$Total <- c(ditDIP, ditDIP, ditDIP, ditDIP, ditDIP, ditDIP, ditDIP, ditDIP, ditDIP, ditDIP, ditDIP, ditVRM, ditVRM, ditVRM, ditVRM, ditVRM, ditVRM, ditVRM, ditVRM, ditVRM, ditVRM, ditVRM)
barFrame$Proportion <- ((barFrame[1]/barFrame[4]) *100)
barFrame$Proportion <- unlist(barFrame$Proportion)
stackPlot <- ggplot(barFrame, aes(x=Tool, y=Proportion, fill=Segment)) + geom_bar(aes(fill=Segment), position='stack', stat="identity") + theme(axis.text.x = element_text(vjust=0.7)) + ylab("Proportion of Identified DIPs %") + ylim(0,100) + scale_fill_manual(values=cols)


# Find average deletion size, on same segments
dataSub <- subset(data_33, BP_seg == RI_seg)
dataSub$delSize <- with(dataSub, abs(BP - RI))
summary(dataSub$delSize)

dataSub <- subset(data_36, BP_seg == RI_seg)
dataSub$delSize <- with(dataSub, abs(BP - RI))
summary(dataSub$delSize)

dataSub <- subset(data_38, BP_seg == RI_seg)
dataSub$delSize <- with(dataSub, abs(BP - RI))
summary(dataSub$delSize)

dataSub <- subset(data_40, BP_seg == RI_seg)
dataSub$delSize <- with(dataSub, abs(BP - RI))
summary(dataSub$delSize)

dataSub <- subset(data_isol, BP_seg == RI_seg)
dataSub$delSize <- with(dataSub, abs(BP - RI))
summary(dataSub$delSize)

# How many found in both
data <- subset(data, DITector_count > 0 & ViReMa_count > 0)



#### ---- Filter ViReMa bed file  ----

bed_33 <- read.table("Bld_33_S1_VRM.bed", sep = "\t", skip = 1)
bed_33 <- subset(bed_33, bed_33[,5] > 5)
write.table(bed_33, "Bld_33_S1_VRM_filt.bed", sep="\t", row.names=FALSE, col.names=FALSE, quote = FALSE) 

bed_36 <- read.table("Bld_36_S2_VRM.bed", sep = "\t", skip = 1)
bed_36 <- subset(bed_36, bed_36[,5] > 5)
write.table(bed_36, "Bld_36_S2_VRM_filt.bed", sep="\t", row.names=FALSE, col.names=FALSE, quote = FALSE) 

bed_38 <- read.table("Bld_38_S4_VRM.bed", sep = "\t", skip = 1)
bed_38 <- subset(bed_38, bed_38[,5] > 5)
write.table(bed_38, "Bld_38_S4_VRM_filt.bed", sep="\t", row.names=FALSE, col.names=FALSE, quote = FALSE) 

bed_40 <- read.table("Bld_40_S6_VRM.bed", sep = "\t", skip = 1)
bed_40 <- subset(bed_40, bed_40[,5] > 5)
write.table(bed_40, "Bld_40_S6_VRM_filt.bed", sep="\t", row.names=FALSE, col.names=FALSE, quote = FALSE) 



#### ---- Venn Diagram comparison  ----

# Read in compiled output from CompDIP.py
compiled <- read.table("compiledOutput.txt", sep = "\t", header = TRUE)


# Create Venn Diagram for DI-Tector reads
dit.venn <- select(compiled, -contains('VRM')) # Filter ViReMa reads
dit.venn <- subset(dit.venn, (s1_DIT > 0 | s2_DIT > 0 | s3_DIT > 0 | s4_DIT > 0))
dit.venn.ls <- list()

# Create list of DIPs for Venn Diagram
for(i in 5:8){
  samV <- c()
  for(n in 1:nrow(dit.venn)){
    DIP <- ""
    DIP <- paste(as.character(dit.venn[n,1]), as.character(dit.venn[n,2]), as.character(dit.venn[n,3]), as.character(dit.venn[n,4]),sep=",")
    if(dit.venn[n,i] > 0){
        samV <- c(samV,DIP)
    }
  }
  len <- length(dit.venn.ls)
  dit.venn.ls[[len+1]] <- samV
}

venn.diagram(x = dit.venn.ls, filename = 'DIT_venn_diagramm.png', category.names = c("Bld_33_S1" , "Bld_36_S2" , "Bld_38_S4", "Bld_40_S6"), output=TRUE,
             lwd = 1, lty = 'blank', col=c("cyan3","darkorchid","brown2", "darkorange1"), fill = c(alpha("cyan3",0.3), alpha("darkorchid",0.3), alpha("brown2",0.3),alpha("darkorange1",0.3)))

# Total DIPs found in each sample 
nrow(dit.venn[dit.venn$s1_DIT > 0,])
nrow(dit.venn[dit.venn$s2_DIT > 0,])
nrow(dit.venn[dit.venn$s3_DIT > 0,])
nrow(dit.venn[dit.venn$s4_DIT > 0,])

# Create Venn Diagram for ViReMa reads
vrm.venn <- select(compiled, -contains('DIT'))
vrm.venn <- subset(vrm.venn, (s1_VRM > 0 | s2_VRM > 0 | s3_VRM > 0 | s4_VRM > 0))
vrm.venn.ls <- list()

for(i in 5:8){
  samV <- c()
  for(n in 1:nrow(vrm.venn)){
    DIP <- ""
    DIP <- paste(as.character(vrm.venn[n,1]), as.character(vrm.venn[n,2]), as.character(vrm.venn[n,3]), as.character(vrm.venn[n,4]),sep=",")
    if(vrm.venn[n,i] > 0){
        samV <- c(samV,DIP)
    }
  }
  len <- length(vrm.venn.ls)
  vrm.venn.ls[[len+1]] <- samV
}
venn.diagram(x = vrm.venn.ls, filename = 'VRM_venn_diagramm.png', category.names = c("Bld_33_S1" , "Bld_36_S2" , "Bld_38_S4", "Bld_40_S6"), output=TRUE,
             lwd = 1, col=c("cyan3","darkorchid","brown2", "darkorange1"), fill = c(alpha("cyan3",0.3), alpha("darkorchid",0.3), alpha("brown2",0.3),alpha("darkorange1",0.3)))

# Total DIPs found in each sample 
nrow(vrm.venn[vrm.venn$s1_VRM > 0,])
nrow(vrm.venn[vrm.venn$s2_VRM > 0,])
nrow(vrm.venn[vrm.venn$s3_VRM > 0,])
nrow(vrm.venn[vrm.venn$s4_VRM > 0,])


# Create Venn Diagram for all/combined reads

# Combine counts for ViReMa and DI-Tector 
compiled$s1_tot <- compiled$s1_DIT + compiled$s1_VRM
compiled$s2_tot <- compiled$s2_DIT + compiled$s2_VRM
compiled$s3_tot <- compiled$s3_DIT + compiled$s3_VRM
compiled$s4_tot <- compiled$s4_DIT + compiled$s4_VRM

combined.venn <- data.frame(compiled$BP_seg, compiled$BP, compiled$RI_seg, compiled$RI, compiled$s1_tot, compiled$s2_tot, compiled$s3_tot, compiled$s4_tot)
combined.venn.ls <- list()

for(i in 5:8){
  samV <- c()
  for(n in 1:nrow(combined.venn)){
    DIP <- ""
    DIP <- paste(as.character(combined.venn[n,1]), as.character(combined.venn[n,2]), as.character(combined.venn[n,3]), as.character(combined.venn[n,4]),sep=",")
    if(combined.venn[n,i] > 0){
      samV <- c(samV,DIP)
    }
  }
  len <- length(combined.venn.ls)
  combined.venn.ls[[len+1]] <- samV
}


venn.diagram(x = combined.venn.ls,  filename = 'DIP_venn_diagramm.png', category.names = c("Bld_33_S1" , "Bld_36_S2" , "Bld_38_S4", "Bld_40_S6"), output=TRUE,
             lwd = 1, col=c("cyan3","darkorchid","brown2", "darkorange1"), fill = c(alpha("cyan3",0.3), alpha("darkorchid",0.3), alpha("brown2",0.3),alpha("darkorange1",0.3)))

# Total DIPs found in each sample 
nrow(combined.venn[combined.venn$compiled.s1_tot > 0,])
nrow(combined.venn[combined.venn$compiled.s2_tot > 0,])
nrow(combined.venn[combined.venn$compiled.s3_tot > 0,])
nrow(combined.venn[combined.venn$compiled.s4_tot > 0,])


# Create Venn Diagram for reads comparing Isol and Bld_38 Sample

# Combine counts for between Bld_38 Isol
isol <- read.table("isolCompiledOutput.txt", sep = "\t", header = TRUE)

isol$s1_tot <- isol$s1_DIT + isol$s1_VRM
isol$s2_tot <- isol$s2_DIT + isol$s2_VRM

isol.venn <- data.frame(isol$BP_seg, isol$BP, isol$RI_seg, isol$RI, isol$s1_tot, isol$s2_tot)
isol.venn.ls <- list()

for(i in 5:6){
  samV <- c()
  for(n in 1:nrow(isol.venn)){
    DIP <- ""
    DIP <- paste(as.character(isol.venn[n,1]), as.character(isol.venn[n,2]), as.character(isol.venn[n,3]), as.character(isol.venn[n,4]),sep=",")
    if(isol.venn[n,i] > 0){
      samV <- c(samV,DIP)
    }
  }
  len <- length(isol.venn.ls)
  isol.venn.ls[[len+1]] <- samV
}


venn.diagram(x = isol.venn.ls,  filename = 'ISOLno_venn_diagramm.png', category.names = c("" , ""), output=TRUE,
             lwd = 1, col=c("cyan3","darkorange1"), fill = c(alpha("cyan3",0.3),alpha("darkorange1",0.3)))

# Total DIPs found in each sample 
nrow(combined.venn[combined.venn$compiled.s1_tot > 0,])
nrow(combined.venn[combined.venn$compiled.s2_tot > 0,])
nrow(combined.venn[combined.venn$compiled.s3_tot > 0,])
nrow(combined.venn[combined.venn$compiled.s4_tot > 0,])






#### ---- Create contingency table for relationship between segements ---- 
segNames <- c('Seg1', 'Seg2', 'Seg3', 'Seg4', 'Seg5', 'Seg6', 'Seg7', 'Seg8', 'Seg9', 'Seg10')

data_33[nrow(data_33)+1,] <- c(5,0,5,0,0,0) # No BP are in seg 5, but need row for table
table_33 <- table(data_33$BP_seg, data_33$RI_seg, dnn=c("BP Seg", "RI Seg"))
row.names(table_33) <- segNames
colnames(table_33) <- segNames
write.table(table_33, "Bld_33_S1_crc.txt", sep="\t", quote = FALSE) 

data_36[nrow(data_36)+1,] <- c(5,0,5,0,0,0) # No BP are in seg 5, but need row for table
table_36 <- table(data_36$BP_seg, data_36$RI_seg)
row.names(table_36) <- segNames
colnames(table_36) <- segNames
write.table(table_36, "Bld_36_S2_crc.txt", sep="\t", quote = FALSE) 

data_38[nrow(data_38)+1,] <- c(5,0,5,0,0,0) # No BP are in seg 5, but need row for table
table_38 <- table(data_38$BP_seg, data_38$RI_seg)
row.names(table_38) <- segNames
colnames(table_38) <- segNames
write.table(table_38, "Bld_38_S4_crc.txt", sep="\t", quote = FALSE) 

data_40[nrow(data_40)+1,] <- c(5,0,5,0,0,0) # No BP are in seg 5, but need row for table
table_40 <- table(data_40$BP_seg, data_40$RI_seg)
row.names(table_40) <- segNames
colnames(table_40) <- segNames
write.table(table_40, "Bld_40_S6_crc.txt", sep="\t", quote = FALSE) 




