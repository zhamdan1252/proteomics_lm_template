#Import packages
library(MSstats)
library(tidyverse)
library(limma)
library(ggplot2)

############ LOAD IN PSMs AND ANNOTATION FILE #################

setwd('..')
setwd(".\\input\\PSM_files")
#aggregate all PSM files to one data frame
files <- list.files(pattern = "\\.csv")
for (file in files) {
  if(!exists("raw_psm")) {
    raw_psm <- read.csv(file)
  }
  
  else {
    temp <- read.csv(file)
    raw_psm <- rbind(raw_psm, temp)
    rm(temp)
  }
}

setwd("..")
setwd("..")
#load in annotations.csv
annot <- read.csv("input\\annotations_file\\annotations.csv")
    
    #  the annotation file should be manually created as described in the 
    #  MSstats user Guide for Proteome Discoverer Data
    

############# COVERT FROM PD TO MSSTATS AND PREP FOR DATA PROCESS ##########

#part of MSstats package: converts Proteome Discoverer file format to MSstats compatible format
input_df <- PDtoMSstatsFormat(raw_psm,
                              annot,
                              which.proteinid = "Master.Protein.Accessions",
                              which.sequence = "Annotated.Sequence",
                              which.quantification = "Precursor.Area"
                              )

#remove any peptides that could be associated with 2 or more Proteins
input_df <- input_df[!grepl(';', input_df$ProteinName),]

#remove any rows with no signal
input_df <- input_df[!is.na(input_df$Intensity),]


############# DATA PROCESSING WITH MSSTATS PACKAGE ###################

#dataProcess is part of the MSstats package. See their reference manual for additional parameters
processed <- dataProcess(input_df,
                         logTrans = 2,
                         normalization = 'FALSE',
                         summaryMethod = 'linear'
                        )

    
######################## GET PROTEIN DATA FRAME FOR LIMMA #######################

#defining function to extract protein data from dataProcess output
get_prot_data <- function(processed_data) {
  #intialize dataframe with number of rows = to number of unique proteins found, 
  #                     and number of columns = number of animals in the study
  prot_data <- data.frame(matrix(ncol = length(unique(processed_data$ProteinLevelData$SUBJECT)), 
                                 nrow = length(unique(processed_data$ProteinLevelData$Protein))))
  
  #setting row and column names in respects to the same variables we set the dimensions to
  rownames(prot_data) <- unique(processed_data$ProteinLevelData$Protein)
  colnames(prot_data) <- unique(processed_data$ProteinLevelData$SUBJECT)
  
  
  prt_lvl <- subset(processed_data$ProteinLevelData, select = c("Protein", "LogIntensities", "SUBJECT"))
  
  #translate the data frame such that the rows = protein name, columns = animals in experiment, cells = log2Intensity
  for(i in 1:nrow(prt_lvl)) {
    prt_name = as.character(prt_lvl[i, "Protein"])
    subject = as.character(prt_lvl[i, "SUBJECT"])
    abundance = as.numeric(prt_lvl[i, "LogIntensities"])
    prot_data[prt_name, subject] = abundance
  }
  #removes animals that produced NA results for ALL proteins
  prot_data <- prot_data[colSums(!is.na(prot_data)) > 0]

  return(prot_data)
}

protein_data = get_prot_data(processed)


################## GET ANNOTATION DATA FRAME FOR LINEAR MODELLING ##############

make_lm_annot <- function(prot_data, annot) {
  
  #initialize dataframe with rows = animals in exp, column = "Condition" 
  lm_annot  <- data.frame(matrix(nrow = ncol(prot_data), ncol = 1))
  rownames(lm_annot) = colnames(prot_data)
  colnames(lm_annot) = c("Condition")
  
  for (i in 1:nrow(lm_annot)) {
    row_name = rownames(lm_annot)[i]
    cond_index = which(annot$BioReplicate == row_name)
    lm_annot$Condition[i] = annot$Condition[cond_index]
  }
  
  return(lm_annot)
}

lin_var_annot <- make_lm_annot(protein_data, annot)


################### LINEAR REGRESSION ###########################

lin_regression <- function(prot_data, condition_annot) {
  
  condition_annot$Condition <- gsub(" ", "_", condition_annot$Condition)
  
  #creates a matrix where the rows = animals in exp and columns = conditions
  #matrix is filled with 1s and 0s to say whether the animal has the condition or not
  design_matrix <- model.matrix(~0 + Condition, data = condition_annot)
  for ( col in 1:ncol(design_matrix)){
    colnames(design_matrix)[col] <-  sub("Condition", "", colnames(design_matrix)[col])
  }
  prot_data[] <- sapply(prot_data, as.numeric)
  
  fit_matrix <- limma::lmFit(prot_data, design_matrix)
  
  #define which conditions you want to compare results 
  contrast_matrix <- limma::makeContrasts(
                                          cond1_vs_cond2 = cond1 - cond2,
                                          cond1_vs_cond3 = cond1 - cond3,
                                          cond2_vs_cond3 = cond2 - cond3,
                                          levels = design_matrix
                                          )
  contrast_fit_matrix <- limma::contrasts.fit(fit_matrix, contrast_matrix)

  efit<-limma::eBayes(contrast_fit_matrix)
  return(efit)
}

fit_data <- lin_regression(protein_data, lin_var_annot)


############## PLOT AND SAVE DATA ###############

#save most top 50 most significant proteins from 1 comparison

setwd("results")
write.csv(topTable(fit_data, coef = "cond1_vs_cond2", number = 50), "cond1_vs_cond2.csv", row.names=TRUE)

pdf("volcanoPlot.pdf")
volcanoplot(fit_data, coef = "cond1_vs_cond2", style = "p-value", highlight = 0, names = fit_data$genes$ID, hl.col="blue",
            xlab = "Log2 Fold Change", ylab = NULL, pch=16, cex=0.35)
dev.off()


