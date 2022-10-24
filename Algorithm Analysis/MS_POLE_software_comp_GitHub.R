#!/usr/bin/env Rscript

print("================> Bioinformatic Tools Analysis MMR-POLE <================")

library(ggplot2)
library(ggpubr)
library(ggh4x)

setwd("mutational_signatures")
MMR_signatures <- c("Signature_6","Signature_14","Signature_15","Signature_20","Signature_21","Signature_26")
POLE_signatures <- c("Signature_10")
MSI_CRC <- read.table("MSI_Cell_Lines_230")
MSS_CRC <- read.table("MSS_Cell_Lines_230")
MSI_TCGA <- read.table("TCGA_MSI")
MSS_TCGA <- read.table("TCGA_MSS")
MSS_CRC <- as.vector(MSS_CRC$V1)
POLE_MUT_CRC <- c("HCA24", "HCC2998", "HDC114", "HROC69", "HT115", "SNU81", "VACO400")
POLE_WT_CRC <- setdiff(MSS_CRC, POLE_MUT_CRC)
POLE_MUT_TCGA <- c("TCGA-AA-3510", "TCGA-AA-3678", "TCGA-AA-3977", "TCGA-AA-3984", "TCGA-AZ-4315")
MSS_TCGA <- as.vector(MSS_TCGA$V1)
POLE_WT_TCGA <- setdiff(MSS_TCGA, POLE_MUT_TCGA)

DS_CRC <- read.table("deconstructSigs_Fitting_Norm_CRCBank_Cosmic2.txt", header = T, row.names = 1, sep = "\t")
DS_TCGA <- read.table("deconstructSigs_Fitting_Norm_TCGA_Cosmic2.txt", header = T, row.names = 1, sep = "\t")

MP_CRC <- read.table("MutationalPatterns_Fitting_Norm_CRCBank_Cosmic2.txt", header = T, row.names = 1, sep = "\t")
MP_TCGA <- read.table("MutationalPatterns_Fitting_Norm_TCGA_Cosmic2.txt", header = T, row.names = 1, sep = "\t")

STL_CRC <- read.table("signature.tools.lib_Fitting_Norm_CRCBank_Cosmic2.txt", header = T, row.names = 1, sep = "\t")
STL_TCGA <- read.table("signature.tools.lib_Fitting_Norm_TCGA_Cosmic2.txt", header = T, row.names = 1, sep = "\t")

SA_CRC <- read.table("signatureanalyzer_Fitting_Norm_CRCBank_Cosmic2.txt", header = T,row.names = 1, sep = "\t")
SA_TCGA <- read.table("signatureanalyzer_Fitting_Norm_TCGA_Cosmic2.txt", header = T, row.names = 1, sep = "\t")

SPA_CRC <- read.table("sigprofilerassignment_Fitting_Norm_CRCBank_Cosmic2.txt", header = T, row.names = 1, sep = "\t")
SPA_TCGA <- read.table("sigprofilerassignment_Fitting_Norm_TCGA_Cosmic2.txt", header = T, row.names = 1, sep = "\t")


###Preclinical Dataset Analysis

##MS Status

MSS_DS_CRC <- apply(DS_CRC[MMR_signatures, MSS_CRC$V1], 2, sum)
MSI_DS_CRC <- apply(DS_CRC[MMR_signatures, MSI_CRC$V1], 2, sum)

MSS_MP_CRC <- apply(MP_CRC[MMR_signatures, MSS_CRC$V1], 2, sum)
MSI_MP_CRC <- apply(MP_CRC[MMR_signatures, MSI_CRC$V1], 2, sum)

MSS_STL_CRC <- apply(STL_CRC[MMR_signatures, MSS_CRC$V1], 2, sum)
MSI_STL_CRC <- apply(STL_CRC[MMR_signatures, MSI_CRC$V1], 2, sum)

MSS_SA_CRC <- apply(SA_CRC[MMR_signatures, MSS_CRC$V1], 2, sum)
MSI_SA_CRC <- apply(SA_CRC[MMR_signatures, MSI_CRC$V1], 2, sum)

MSS_SPA_CRC <- apply(SPA_CRC[MMR_signatures, MSS_CRC$V1], 2, sum)
MSI_SPA_CRC <- apply(SPA_CRC[MMR_signatures, MSI_CRC$V1], 2, sum)


DF_DS_MSS <- data.frame(matrix(nrow = 152, ncol = 3))
colnames(DF_DS_MSS) <- c("MMR_signatures", "MS_Status", "Tool")
DF_DS_MSS$MMR_signatures <- MSS_DS_CRC
DF_DS_MSS$MS_Status <- "MSS"
DF_DS_MSS$Tool <- "deconstructSigs"

DF_DS_MSI <- data.frame(matrix(nrow = 78, ncol = 3))
colnames(DF_DS_MSI) <- c("MMR_signatures", "MS_Status", "Tool")
DF_DS_MSI$MMR_signatures <- MSI_DS_CRC
DF_DS_MSI$MS_Status <- "MSI"
DF_DS_MSI$Tool <- "deconstructSigs"

###

DF_MP_MSS <- data.frame(matrix(nrow = 152, ncol = 3))
colnames(DF_MP_MSS) <- c("MMR_signatures", "MS_Status", "Tool")
DF_MP_MSS$MMR_signatures <- MSS_MP_CRC
DF_MP_MSS$MS_Status <- "MSS"
DF_MP_MSS$Tool <- "MutationalPatterns"

DF_MP_MSI <- data.frame(matrix(nrow = 78, ncol = 3))
colnames(DF_MP_MSI) <- c("MMR_signatures", "MS_Status", "Tool")
DF_MP_MSI$MMR_signatures <- MSI_MP_CRC
DF_MP_MSI$MS_Status <- "MSI"
DF_MP_MSI$Tool <- "MutationalPatterns"

###

DF_STL_MSS <- data.frame(matrix(nrow = 152, ncol = 3))
colnames(DF_STL_MSS) <- c("MMR_signatures", "MS_Status", "Tool")
DF_STL_MSS$MMR_signatures <- MSS_STL_CRC
DF_STL_MSS$MS_Status <- "MSS"
DF_STL_MSS$Tool <- "signature.tools.lib"

DF_STL_MSI <- data.frame(matrix(nrow = 78, ncol = 3))
colnames(DF_STL_MSI) <- c("MMR_signatures", "MS_Status", "Tool")
DF_STL_MSI$MMR_signatures <- MSI_STL_CRC
DF_STL_MSI$MS_Status <- "MSI"
DF_STL_MSI$Tool <- "signature.tools.lib"

###

DF_SA_MSS <- data.frame(matrix(nrow = 152, ncol = 3))
colnames(DF_SA_MSS) <- c("MMR_signatures", "MS_Status", "Tool")
DF_SA_MSS$MMR_signatures <- MSS_SA_CRC
DF_SA_MSS$MS_Status <- "MSS"
DF_SA_MSS$Tool <- "SignatureAnalyzer"

DF_SA_MSI <- data.frame(matrix(nrow = 78, ncol = 3))
colnames(DF_SA_MSI) <- c("MMR_signatures", "MS_Status", "Tool")
DF_SA_MSI$MMR_signatures <- MSI_SA_CRC
DF_SA_MSI$MS_Status <- "MSI"
DF_SA_MSI$Tool <- "SignatureAnalyzer"

###

DF_SPA_MSS <- data.frame(matrix(nrow = 152, ncol = 3))
colnames(DF_SPA_MSS) <- c("MMR_signatures", "MS_Status", "Tool")
DF_SPA_MSS$MMR_signatures <- MSS_SPA_CRC
DF_SPA_MSS$MS_Status <- "MSS"
DF_SPA_MSS$Tool <- "SigProfilerAssignment"

DF_SPA_MSI <- data.frame(matrix(nrow = 78, ncol = 3))
colnames(DF_SPA_MSI) <- c("MMR_signatures", "MS_Status", "Tool")
DF_SPA_MSI$MMR_signatures <- MSI_SPA_CRC
DF_SPA_MSI$MS_Status <- "MSI"
DF_SPA_MSI$Tool <- "SigProfilerAssignment"


DF_MS_CRC <- rbind(DF_DS_MSI, DF_DS_MSS, DF_MP_MSI, DF_MP_MSS, DF_SA_MSI, DF_SA_MSS, DF_SPA_MSI, DF_SPA_MSS, DF_STL_MSI, DF_STL_MSS)

DF_MS_CRC$Tool <- factor(DF_MS_CRC$Tool, levels = c("MutationalPatterns", "deconstructSigs", "signature.tools.lib", "SigProfilerAssignment", "SignatureAnalyzer"))

ggplot(DF_MS_CRC, aes(x=Tool, y=MMR_signatures, fill=MS_Status)) + 
  geom_boxplot(aes(fill=MS_Status)) +
  scale_fill_manual(values=c("lightblue1", "salmon2")) +
  scale_y_continuous(breaks = seq(0, 1, 0.25),
                     limits = c(0,1),
                     expand = c(0.01, 0.01))+
  force_panelsizes(row=unit(6, "cm"), cols = unit(8, "cm")) +
  labs(x= "Software", y="MMR Signatures Contribution") + 
  theme(axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10), legend.text = element_text(size = 10), 
        axis.text.x = element_text(angle = -45, vjust = 0, hjust = 0.1),
        legend.title = element_text(size = 10), axis.text = element_text(size = 10),
        panel.background = element_rect(fill = "White", colour = "Black", linetype = 1, size = 0.4),
        panel.grid.major = element_line(colour = "Black", size = 0.1)) +
  stat_compare_means(method ="wilcox.test", label = "p.signif")


##POLE Status

POLE_WT_DS_CRC <- apply(DS_CRC[POLE_signatures, POLE_WT_CRC], 2, sum)
POLE_MUT_DS_CRC <- apply(DS_CRC[POLE_signatures, POLE_MUT_CRC], 2, sum)

POLE_WT_MP_CRC <- apply(MP_CRC[POLE_signatures, POLE_WT_CRC], 2, sum)
POLE_MUT_MP_CRC <- apply(MP_CRC[POLE_signatures, POLE_MUT_CRC], 2, sum)

POLE_WT_STL_CRC <- apply(STL_CRC[POLE_signatures, POLE_WT_CRC], 2, sum)
POLE_MUT_STL_CRC <- apply(STL_CRC[POLE_signatures, POLE_MUT_CRC], 2, sum)

POLE_WT_SA_CRC <- apply(SA_CRC[POLE_signatures, POLE_WT_CRC], 2, sum)
POLE_MUT_SA_CRC <- apply(SA_CRC[POLE_signatures, POLE_MUT_CRC], 2, sum)

POLE_WT_SPA_CRC <- apply(SPA_CRC[POLE_signatures, POLE_WT_CRC], 2, sum)
POLE_MUT_SPA_CRC <- apply(SPA_CRC[POLE_signatures, POLE_MUT_CRC], 2, sum)


DF_DS_POLE_WT <- data.frame(matrix(nrow = 145, ncol = 3))
colnames(DF_DS_POLE_WT) <- c("POLE_signatures", "POLE_Status", "Tool")
DF_DS_POLE_WT$POLE_signatures <- POLE_WT_DS_CRC
DF_DS_POLE_WT$POLE_Status <- "POLE_WT"
DF_DS_POLE_WT$Tool <- "deconstructSigs"

DF_DS_POLE_MUT <- data.frame(matrix(nrow = 7, ncol = 3))
colnames(DF_DS_POLE_MUT) <- c("POLE_signatures", "POLE_Status", "Tool")
DF_DS_POLE_MUT$POLE_signatures <- POLE_MUT_DS_CRC
DF_DS_POLE_MUT$POLE_Status <- "POLE_MUT"
DF_DS_POLE_MUT$Tool <- "deconstructSigs"

###

DF_MP_POLE_WT <- data.frame(matrix(nrow = 145, ncol = 3))
colnames(DF_MP_POLE_WT) <- c("POLE_signatures", "POLE_Status", "Tool")
DF_MP_POLE_WT$POLE_signatures <- POLE_WT_MP_CRC
DF_MP_POLE_WT$POLE_Status <- "POLE_WT"
DF_MP_POLE_WT$Tool <- "MutationalPatterns"

DF_MP_POLE_MUT <- data.frame(matrix(nrow = 7, ncol = 3))
colnames(DF_MP_POLE_MUT) <- c("POLE_signatures", "POLE_Status", "Tool")
DF_MP_POLE_MUT$POLE_signatures <- POLE_MUT_MP_CRC
DF_MP_POLE_MUT$POLE_Status <- "POLE_MUT"
DF_MP_POLE_MUT$Tool <- "MutationalPatterns"

###

DF_STL_POLE_WT <- data.frame(matrix(nrow = 145, ncol = 3))
colnames(DF_STL_POLE_WT) <- c("POLE_signatures", "POLE_Status", "Tool")
DF_STL_POLE_WT$POLE_signatures <- POLE_WT_STL_CRC
DF_STL_POLE_WT$POLE_Status <- "POLE_WT"
DF_STL_POLE_WT$Tool <- "signature.tools.lib"

DF_STL_POLE_MUT <- data.frame(matrix(nrow = 7, ncol = 3))
colnames(DF_STL_POLE_MUT) <- c("POLE_signatures", "POLE_Status", "Tool")
DF_STL_POLE_MUT$POLE_signatures <- POLE_MUT_STL_CRC
DF_STL_POLE_MUT$POLE_Status <- "POLE_MUT"
DF_STL_POLE_MUT$Tool <- "signature.tools.lib"

###

DF_SA_POLE_WT <- data.frame(matrix(nrow = 145, ncol = 3))
colnames(DF_SA_POLE_WT) <- c("POLE_signatures", "POLE_Status", "Tool")
DF_SA_POLE_WT$POLE_signatures <- POLE_WT_SA_CRC
DF_SA_POLE_WT$POLE_Status <- "POLE_WT"
DF_SA_POLE_WT$Tool <- "SignatureAnalyzer"

DF_SA_POLE_MUT <- data.frame(matrix(nrow = 7, ncol = 3))
colnames(DF_SA_POLE_MUT) <- c("POLE_signatures", "POLE_Status", "Tool")
DF_SA_POLE_MUT$POLE_signatures <- POLE_MUT_SA_CRC
DF_SA_POLE_MUT$POLE_Status <- "POLE_MUT"
DF_SA_POLE_MUT$Tool <- "SignatureAnalyzer"

###

DF_SPA_POLE_WT <- data.frame(matrix(nrow = 145, ncol = 3))
colnames(DF_SPA_POLE_WT) <- c("POLE_signatures", "POLE_Status", "Tool")
DF_SPA_POLE_WT$POLE_signatures <- POLE_WT_SPA_CRC
DF_SPA_POLE_WT$POLE_Status <- "POLE_WT"
DF_SPA_POLE_WT$Tool <- "SigProfilerAssignment"

DF_SPA_POLE_MUT <- data.frame(matrix(nrow = 7, ncol = 3))
colnames(DF_SPA_POLE_MUT) <- c("POLE_signatures", "POLE_Status", "Tool")
DF_SPA_POLE_MUT$POLE_signatures <- POLE_MUT_SPA_CRC
DF_SPA_POLE_MUT$POLE_Status <- "POLE_MUT"
DF_SPA_POLE_MUT$Tool <- "SigProfilerAssignment"


DF_POLE_CRC <- rbind(DF_DS_POLE_MUT, DF_DS_POLE_WT, DF_MP_POLE_MUT, DF_MP_POLE_WT, DF_SA_POLE_MUT, DF_SA_POLE_WT, DF_SPA_POLE_MUT, DF_SPA_POLE_WT, DF_STL_POLE_MUT, DF_STL_POLE_WT)

DF_POLE_CRC$Tool <- factor(DF_POLE_CRC$Tool, levels = c("MutationalPatterns", "deconstructSigs", "signature.tools.lib", "SigProfilerAssignment", "SignatureAnalyzer"))

ggplot(DF_POLE_CRC, aes(x=Tool, y=POLE_signatures, fill=POLE_Status)) + 
  geom_boxplot(aes(fill=POLE_Status)) +
  scale_fill_manual(values=c("darkgoldenrod1", "olivedrab4")) +
  scale_y_continuous(breaks = seq(0, 1, 0.25),
                     limits = c(0,1),
                     expand = c(0, 0)) +
  force_panelsizes(row=unit(6, "cm"), cols = unit(8, "cm")) +
  labs(x= "Software", y="POLE Signatures Contribution") + 
  theme_bw() +
  theme(axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10), legend.text = element_text(size = 10), 
        axis.text.x = element_text(angle = -45, vjust = 0, hjust = 0.1),
        legend.title = element_text(size = 10), axis.text = element_text(size = 10),
        panel.background = element_rect(fill = "White", colour = "Black", linetype = 1, size = 0.4),
        panel.grid.major = element_line(colour = "Black", size = 0.1)) +
  stat_compare_means(method ="wilcox.test", label = "p.signif")



###Clinical Dataset Analysis

#MS Status

colnames(DS_TCGA) <- gsub('\\.', '-', colnames(DS_TCGA))
colnames(DS_TCGA) <- substr(colnames(DS_TCGA), 1, nchar(colnames(DS_TCGA))-16)

colnames(MP_TCGA) <- gsub('\\.', '-', colnames(MP_TCGA))
colnames(MP_TCGA) <- substr(colnames(MP_TCGA), 1, nchar(colnames(MP_TCGA))-16)

colnames(STL_TCGA) <- gsub('\\.', '-', colnames(STL_TCGA))
colnames(STL_TCGA) <- substr(colnames(STL_TCGA), 1, nchar(colnames(STL_TCGA))-16)

colnames(SA_TCGA) <- gsub('\\.', '-', colnames(SA_TCGA))
colnames(SA_TCGA) <- substr(colnames(SA_TCGA), 1, nchar(colnames(SA_TCGA))-16)

colnames(SPA_TCGA) <- gsub('\\.', '-', colnames(SPA_TCGA))
colnames(SPA_TCGA) <- substr(colnames(SPA_TCGA), 1, nchar(colnames(SPA_TCGA))-16)



MSS_DS_TCGA <- apply(DS_TCGA[MMR_signatures, MSS_TCGA$V1], 2, sum)
MSI_DS_TCGA <- apply(DS_TCGA[MMR_signatures, MSI_TCGA$V1], 2, sum)

MSS_MP_TCGA <- apply(MP_TCGA[MMR_signatures, MSS_TCGA$V1], 2, sum)
MSI_MP_TCGA <- apply(MP_TCGA[MMR_signatures, MSI_TCGA$V1], 2, sum)

MSS_STL_TCGA <- apply(STL_TCGA[MMR_signatures, MSS_TCGA$V1], 2, sum)
MSI_STL_TCGA <- apply(STL_TCGA[MMR_signatures, MSI_TCGA$V1], 2, sum)

MSS_SA_TCGA <- apply(SA_TCGA[MMR_signatures, MSS_TCGA$V1], 2, sum)
MSI_SA_TCGA <- apply(SA_TCGA[MMR_signatures, MSI_TCGA$V1], 2, sum)

MSS_SPA_TCGA <- apply(SPA_TCGA[MMR_signatures, MSS_TCGA$V1], 2, sum)
MSI_SPA_TCGA <- apply(SPA_TCGA[MMR_signatures, MSI_TCGA$V1], 2, sum)


DF_DS_MSS_TCGA <- data.frame(matrix(nrow = 64, ncol = 3))
colnames(DF_DS_MSS_TCGA) <- c("MMR_signatures", "MS_Status", "Tool")
DF_DS_MSS_TCGA$MMR_signatures <- MSS_DS_TCGA
DF_DS_MSS_TCGA$MS_Status <- "MSS"
DF_DS_MSS_TCGA$Tool <- "deconstructSigs"

DF_DS_MSI_TCGA <- data.frame(matrix(nrow = 68, ncol = 3))
colnames(DF_DS_MSI_TCGA) <- c("MMR_signatures", "MS_Status", "Tool")
DF_DS_MSI_TCGA$MMR_signatures <- MSI_DS_TCGA
DF_DS_MSI_TCGA$MS_Status <- "MSI"
DF_DS_MSI_TCGA$Tool <- "deconstructSigs"

###

DF_MP_MSS_TCGA <- data.frame(matrix(nrow = 64, ncol = 3))
colnames(DF_MP_MSS_TCGA) <- c("MMR_signatures", "MS_Status", "Tool")
DF_MP_MSS_TCGA$MMR_signatures <- MSS_MP_TCGA
DF_MP_MSS_TCGA$MS_Status <- "MSS"
DF_MP_MSS_TCGA$Tool <- "MutationalPatterns"

DF_MP_MSI_TCGA <- data.frame(matrix(nrow = 68, ncol = 3))
colnames(DF_MP_MSI_TCGA) <- c("MMR_signatures", "MS_Status", "Tool")
DF_MP_MSI_TCGA$MMR_signatures <- MSI_MP_TCGA
DF_MP_MSI_TCGA$MS_Status <- "MSI"
DF_MP_MSI_TCGA$Tool <- "MutationalPatterns"

###

DF_STL_MSS_TCGA <- data.frame(matrix(nrow = 64, ncol = 3))
colnames(DF_STL_MSS_TCGA) <- c("MMR_signatures", "MS_Status", "Tool")
DF_STL_MSS_TCGA$MMR_signatures <- MSS_STL_TCGA
DF_STL_MSS_TCGA$MS_Status <- "MSS"
DF_STL_MSS_TCGA$Tool <- "signature.tools.lib"

DF_STL_MSI_TCGA <- data.frame(matrix(nrow = 68, ncol = 3))
colnames(DF_STL_MSI_TCGA) <- c("MMR_signatures", "MS_Status", "Tool")
DF_STL_MSI_TCGA$MMR_signatures <- MSI_STL_TCGA
DF_STL_MSI_TCGA$MS_Status <- "MSI"
DF_STL_MSI_TCGA$Tool <- "signature.tools.lib"

###

DF_SA_MSS_TCGA <- data.frame(matrix(nrow = 64, ncol = 3))
colnames(DF_SA_MSS_TCGA) <- c("MMR_signatures", "MS_Status", "Tool")
DF_SA_MSS_TCGA$MMR_signatures <- MSS_SA_TCGA
DF_SA_MSS_TCGA$MS_Status <- "MSS"
DF_SA_MSS_TCGA$Tool <- "SignatureAnalyzer"

DF_SA_MSI_TCGA <- data.frame(matrix(nrow = 68, ncol = 3))
colnames(DF_SA_MSI_TCGA) <- c("MMR_signatures", "MS_Status", "Tool")
DF_SA_MSI_TCGA$MMR_signatures <- MSI_SA_TCGA
DF_SA_MSI_TCGA$MS_Status <- "MSI"
DF_SA_MSI_TCGA$Tool <- "SignatureAnalyzer"

###

DF_SPA_MSS_TCGA <- data.frame(matrix(nrow = 64, ncol = 3))
colnames(DF_SPA_MSS_TCGA) <- c("MMR_signatures", "MS_Status", "Tool")
DF_SPA_MSS_TCGA$MMR_signatures <- MSS_SPA_TCGA
DF_SPA_MSS_TCGA$MS_Status <- "MSS"
DF_SPA_MSS_TCGA$Tool <- "SigProfilerAssignment"

DF_SPA_MSI_TCGA <- data.frame(matrix(nrow = 68, ncol = 3))
colnames(DF_SPA_MSI_TCGA) <- c("MMR_signatures", "MS_Status", "Tool")
DF_SPA_MSI_TCGA$MMR_signatures <- MSI_SPA_TCGA
DF_SPA_MSI_TCGA$MS_Status <- "MSI"
DF_SPA_MSI_TCGA$Tool <- "SigProfilerAssignment"


DF_MS_TCGA <- rbind(DF_DS_MSI_TCGA, DF_DS_MSS_TCGA, DF_MP_MSI_TCGA, DF_MP_MSS_TCGA, DF_SA_MSI_TCGA, DF_SA_MSS_TCGA, DF_SPA_MSI_TCGA, DF_SPA_MSS_TCGA, DF_STL_MSI_TCGA, DF_STL_MSS_TCGA)

DF_MS_TCGA$Tool <- factor(DF_MS_TCGA$Tool, levels = c("MutationalPatterns", "deconstructSigs", "signature.tools.lib", "SigProfilerAssignment", "SignatureAnalyzer"))

ggplot(DF_MS_TCGA, aes(x=Tool, y=MMR_signatures, fill=MS_Status)) + 
  geom_boxplot(aes(fill=MS_Status)) +
  labs(title = "MS classification software comparison TCGA", x= "Software", y="MMR Signatures Contribution") + 
  scale_fill_manual(values=c("lightblue1", "salmon2")) +
  scale_y_continuous(breaks = seq(0, 1, 0.25),
                     limits = c(0,1),
                     expand = c(0.01, 0.01))+
  force_panelsizes(row=unit(6, "cm"), cols = unit(8, "cm")) +
  labs(x= "Software", y="MMR Signatures Contribution") + 
  theme(axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10), legend.text = element_text(size = 10), 
        axis.text.x = element_text(angle = -45, vjust = 0, hjust = 0.1),
        legend.title = element_text(size = 10), axis.text = element_text(size = 10),
        panel.background = element_rect(fill = "White", colour = "Black", linetype = 1, size = 0.4),
        panel.grid.major = element_line(colour = "Black", size = 0.1)) +
  stat_compare_means(method ="wilcox.test", label = "p.signif")

#POLE Status

POLE_WT_DS_TCGA <- apply(DS_TCGA[POLE_signatures, POLE_WT_TCGA], 2, sum)
POLE_MUT_DS_TCGA <- apply(DS_TCGA[POLE_signatures, POLE_MUT_TCGA], 2, sum)

POLE_WT_MP_TCGA <- apply(MP_TCGA[POLE_signatures, POLE_WT_TCGA], 2, sum)
POLE_MUT_MP_TCGA <- apply(MP_TCGA[POLE_signatures, POLE_MUT_TCGA], 2, sum)

POLE_WT_STL_TCGA <- apply(STL_TCGA[POLE_signatures, POLE_WT_TCGA], 2, sum)
POLE_MUT_STL_TCGA <- apply(STL_TCGA[POLE_signatures, POLE_MUT_TCGA], 2, sum)

POLE_WT_SA_TCGA <- apply(SA_TCGA[POLE_signatures, POLE_WT_TCGA], 2, sum)
POLE_MUT_SA_TCGA <- apply(SA_TCGA[POLE_signatures, POLE_MUT_TCGA], 2, sum)

POLE_WT_SPA_TCGA <- apply(SPA_TCGA[POLE_signatures, POLE_WT_TCGA], 2, sum)
POLE_MUT_SPA_TCGA <- apply(SPA_TCGA[POLE_signatures, POLE_MUT_TCGA], 2, sum)


DF_DS_POLE_WT_TCGA <- data.frame(matrix(nrow = 59, ncol = 3))
colnames(DF_DS_POLE_WT_TCGA) <- c("POLE_signatures", "POLE_Status", "Tool")
DF_DS_POLE_WT_TCGA$POLE_signatures <- POLE_WT_DS_TCGA
DF_DS_POLE_WT_TCGA$POLE_Status <- "POLE_WT"
DF_DS_POLE_WT_TCGA$Tool <- "deconstructSigs"

DF_DS_POLE_MUT_TCGA <- data.frame(matrix(nrow = 5, ncol = 3))
colnames(DF_DS_POLE_MUT_TCGA) <- c("POLE_signatures", "POLE_Status", "Tool")
DF_DS_POLE_MUT_TCGA$POLE_signatures <- POLE_MUT_DS_TCGA
DF_DS_POLE_MUT_TCGA$POLE_Status <- "POLE_MUT"
DF_DS_POLE_MUT_TCGA$Tool <- "deconstructSigs"

###

DF_MP_POLE_WT_TCGA <- data.frame(matrix(nrow = 59, ncol = 3))
colnames(DF_MP_POLE_WT_TCGA) <- c("POLE_signatures", "POLE_Status", "Tool")
DF_MP_POLE_WT_TCGA$POLE_signatures <- POLE_WT_MP_TCGA
DF_MP_POLE_WT_TCGA$POLE_Status <- "POLE_WT"
DF_MP_POLE_WT_TCGA$Tool <- "MutationalPatterns"

DF_MP_POLE_MUT_TCGA <- data.frame(matrix(nrow = 5, ncol = 3))
colnames(DF_MP_POLE_MUT_TCGA) <- c("POLE_signatures", "POLE_Status", "Tool")
DF_MP_POLE_MUT_TCGA$POLE_signatures <- POLE_MUT_MP_TCGA
DF_MP_POLE_MUT_TCGA$POLE_Status <- "POLE_MUT"
DF_MP_POLE_MUT_TCGA$Tool <- "MutationalPatterns"

###

DF_STL_POLE_WT_TCGA <- data.frame(matrix(nrow = 59, ncol = 3))
colnames(DF_STL_POLE_WT_TCGA) <- c("POLE_signatures", "POLE_Status", "Tool")
DF_STL_POLE_WT_TCGA$POLE_signatures <- POLE_WT_STL_TCGA
DF_STL_POLE_WT_TCGA$POLE_Status <- "POLE_WT"
DF_STL_POLE_WT_TCGA$Tool <- "signature.tools.lib"

DF_STL_POLE_MUT_TCGA <- data.frame(matrix(nrow = 5, ncol = 3))
colnames(DF_STL_POLE_MUT_TCGA) <- c("POLE_signatures", "POLE_Status", "Tool")
DF_STL_POLE_MUT_TCGA$POLE_signatures <- POLE_MUT_STL_TCGA
DF_STL_POLE_MUT_TCGA$POLE_Status <- "POLE_MUT"
DF_STL_POLE_MUT_TCGA$Tool <- "signature.tools.lib"

###

DF_SA_POLE_WT_TCGA <- data.frame(matrix(nrow = 59, ncol = 3))
colnames(DF_SA_POLE_WT_TCGA) <- c("POLE_signatures", "POLE_Status", "Tool")
DF_SA_POLE_WT_TCGA$POLE_signatures <- POLE_WT_SA_TCGA
DF_SA_POLE_WT_TCGA$POLE_Status <- "POLE_WT"
DF_SA_POLE_WT_TCGA$Tool <- "SignatureAnalyzer"

DF_SA_POLE_MUT_TCGA <- data.frame(matrix(nrow = 5, ncol = 3))
colnames(DF_SA_POLE_MUT_TCGA) <- c("POLE_signatures", "POLE_Status", "Tool")
DF_SA_POLE_MUT_TCGA$POLE_signatures <- POLE_MUT_SA_TCGA
DF_SA_POLE_MUT_TCGA$POLE_Status <- "POLE_MUT"
DF_SA_POLE_MUT_TCGA$Tool <- "SignatureAnalyzer"

###

DF_SPA_POLE_WT_TCGA <- data.frame(matrix(nrow = 59, ncol = 3))
colnames(DF_SPA_POLE_WT_TCGA) <- c("POLE_signatures", "POLE_Status", "Tool")
DF_SPA_POLE_WT_TCGA$POLE_signatures <- POLE_WT_SPA_TCGA
DF_SPA_POLE_WT_TCGA$POLE_Status <- "POLE_WT"
DF_SPA_POLE_WT_TCGA$Tool <- "SigProfilerAssignment"

DF_SPA_POLE_MUT_TCGA <- data.frame(matrix(nrow = 5, ncol = 3))
colnames(DF_SPA_POLE_MUT_TCGA) <- c("POLE_signatures", "POLE_Status", "Tool")
DF_SPA_POLE_MUT_TCGA$POLE_signatures <- POLE_MUT_SPA_TCGA
DF_SPA_POLE_MUT_TCGA$POLE_Status <- "POLE_MUT"
DF_SPA_POLE_MUT_TCGA$Tool <- "SigProfilerAssignment"


DF_POLE_TCGA <- rbind(DF_DS_POLE_MUT_TCGA, DF_DS_POLE_WT_TCGA, DF_MP_POLE_MUT_TCGA, DF_MP_POLE_WT_TCGA, DF_SA_POLE_MUT_TCGA, DF_SA_POLE_WT_TCGA, DF_SPA_POLE_MUT_TCGA, DF_SPA_POLE_WT_TCGA, DF_STL_POLE_MUT_TCGA, DF_STL_POLE_WT_TCGA)

DF_POLE_TCGA$Tool <- factor(DF_POLE_TCGA$Tool, levels = c("MutationalPatterns", "deconstructSigs", "signature.tools.lib", "SigProfilerAssignment", "SignatureAnalyzer"))

ggplot(DF_POLE_TCGA, aes(x=Tool, y=POLE_signatures, fill=POLE_Status)) + 
  geom_boxplot(aes(fill=POLE_Status)) +
  scale_fill_manual(values=c("darkgoldenrod1", "olivedrab4")) +
  scale_y_continuous(breaks = seq(0, 1, 0.25),
                     limits = c(0,1),
                     expand = c(0, 0)) +
  force_panelsizes(row=unit(6, "cm"), cols = unit(8, "cm")) +
  labs(x= "Software", y="POLE Signatures Contribution") + 
  theme_bw() +
  theme(axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10), legend.text = element_text(size = 10), 
        axis.text.x = element_text(angle = -45, vjust = 0, hjust = 0.1),
        legend.title = element_text(size = 10), axis.text = element_text(size = 10),
        panel.background = element_rect(fill = "White", colour = "Black", linetype = 1, size = 0.4),
        panel.grid.major = element_line(colour = "Black", size = 0.1)) +
  stat_compare_means(method ="wilcox.test", label = "p.signif")



#Stat MMR Preclinical Dataset

delta_MSI_MP <- median(DF_MP_MSI$MMR_signatures)
delta_MSS_MP <- median(DF_MP_MSS$MMR_signatures)

distr_delta_MP <- DF_MP_MSI$MMR_signatures-delta_MSS_MP
distr_delta_MP <- append(distr_delta_MP, abs(DF_MP_MSS$MMR_signatures-delta_MSI_MP))

delta_MSI_DS <- median(DF_DS_MSI$MMR_signatures)
delta_MSS_DS <- median(DF_DS_MSS$MMR_signatures)

distr_delta_DS <- DF_DS_MSI$MMR_signatures-delta_MSS_DS
distr_delta_DS <- append(distr_delta_DS, abs(DF_DS_MSS$MMR_signatures-delta_MSI_DS))

delta_MSI_STL <- median(DF_STL_MSI$MMR_signatures)
delta_MSS_STL <- median(DF_STL_MSS$MMR_signatures)

distr_delta_STL <- DF_STL_MSI$MMR_signatures-delta_MSS_STL
distr_delta_STL <- append(distr_delta_STL, abs(DF_STL_MSS$MMR_signatures-delta_MSI_STL))

delta_MSI_SPA <- median(DF_SPA_MSI$MMR_signatures)
delta_MSS_SPA <- median(DF_SPA_MSS$MMR_signatures)

distr_delta_SPA <- DF_SPA_MSI$MMR_signatures-delta_MSS_SPA
distr_delta_SPA <- append(distr_delta_SPA, abs(DF_SPA_MSS$MMR_signatures-delta_MSI_SPA))

delta_MSI_SA <- median(DF_SA_MSI$MMR_signatures)
delta_MSS_SA <- median(DF_SA_MSS$MMR_signatures)

distr_delta_SA <- DF_SA_MSI$MMR_signatures-delta_MSS_SA
distr_delta_SA <- append(distr_delta_SA, abs(DF_SA_MSS$MMR_signatures-delta_MSI_SA))



DF_violin_MP <- data.frame(matrix(nrow = 230, ncol = 2))
DF_violin_MP$X1 <- abs(distr_delta_MP)
DF_violin_MP$X2 <- "MutationalPatterns"

DF_violin_DS <- data.frame(matrix(nrow = 230, ncol = 2))
DF_violin_DS$X1 <- abs(distr_delta_DS)
DF_violin_DS$X2 <- "deconstructSigs"

DF_violin_STL <- data.frame(matrix(nrow = 230, ncol = 2))
DF_violin_STL$X1 <- abs(distr_delta_STL)
DF_violin_STL$X2 <- "signature.tools.lib"

DF_violin_SPA <- data.frame(matrix(nrow = 230, ncol = 2))
DF_violin_SPA$X1 <- abs(distr_delta_SPA)
DF_violin_SPA$X2 <- "SigProfilerAssignment"

DF_violin_SA <- data.frame(matrix(nrow = 230, ncol = 2))
DF_violin_SA$X1 <- abs(distr_delta_SA)
DF_violin_SA$X2 <- "SignatureAnalyzer"


DF_violin_CRC <- rbind(DF_violin_MP, DF_violin_DS, DF_violin_STL, DF_violin_SPA, DF_violin_SA)
colnames(DF_violin_CRC) <- c("Delta MMR", "Tools")

DF_violin_CRC$Color <- rep("MS", 1150)

couples <- list(c("MutationalPatterns", "deconstructSigs"), c("MutationalPatterns", "signature.tools.lib"), c("MutationalPatterns", "SigProfilerAssignment"), c("MutationalPatterns", "SignatureAnalyzer"), c("signature.tools.lib", "deconstructSigs"), c("signature.tools.lib", "SigProfilerAssignment"), c("signature.tools.lib", "SignatureAnalyzer"), c("SigProfilerAssignment", "deconstructSigs"), c("SignatureAnalyzer", "deconstructSigs"), c("SigProfilerAssignment", "SignatureAnalyzer"))
DF_violin_CRC$Tool <- factor(DF_violin_CRC$Tool, levels = c("MutationalPatterns", "deconstructSigs", "signature.tools.lib", "SigProfilerAssignment", "SignatureAnalyzer"))


ggplot(DF_violin_CRC, aes(x=Tool, y=`Delta MMR`, color=Color, fill=Color)) +
  geom_violin(adjust = 4) +
  scale_color_manual(values=c("brown3")) +
  scale_fill_manual(values=c("coral")) +
  scale_y_continuous(breaks = seq(0, 1, 0.25),
                     limits = c(0,1),
                     expand = c(0, 0)) +
  force_panelsizes(row=unit(6, "cm"), cols = unit(8, "cm")) +
  theme(plot.title = element_text(size = 10, hjust = 0.5, face = "bold"), 
      axis.title.x = element_text(size = 10), 
      axis.title.y = element_text(size = 10),
      axis.text.x = element_text(angle = -45, vjust = 0, hjust = 0.1),
      legend.text = element_text(size = 16), 
      legend.position = "None" ,
      panel.background = element_rect(fill = "White", colour = "Black", linetype = 1, size = 0.8),
      panel.grid.major = element_line(colour = "Black", size = 0.1)) +
  geom_boxplot(width=0.1, colour = "black", fill = "white", outlier.colour = NA)+
  stat_compare_means(method ="wilcox", label = "p.label", ref.group = Tool, comparisons = couples)

#Stat POLE Preclinical Dataset

delta_POLE_MUT_MP <- median(DF_MP_POLE_MUT$POLE_signatures)
delta_POLE_WT_MP <- median(DF_MP_POLE_WT$POLE_signatures)

distr_delta_MP_POLE <- DF_MP_POLE_MUT$POLE_signatures-delta_POLE_WT_MP
distr_delta_MP_POLE <- append(distr_delta_MP_POLE, abs(DF_MP_POLE_WT$POLE_signatures-delta_POLE_MUT_MP))

delta_POLE_MUT_DS <- median(DF_DS_POLE_MUT$POLE_signatures)
delta_POLE_WT_DS <- median(DF_DS_POLE_WT$POLE_signatures)

distr_delta_DS_POLE <- DF_DS_POLE_MUT$POLE_signatures-delta_POLE_WT_DS
distr_delta_DS_POLE <- append(distr_delta_DS_POLE, abs(DF_DS_POLE_WT$POLE_signatures-delta_POLE_MUT_DS))

delta_POLE_MUT_STL <- median(DF_STL_POLE_MUT$POLE_signatures)
delta_POLE_WT_STL <- median(DF_STL_POLE_WT$POLE_signatures)

distr_delta_STL_POLE <- DF_STL_POLE_MUT$POLE_signatures-delta_POLE_WT_STL
distr_delta_STL_POLE <- append(distr_delta_STL_POLE, abs(DF_STL_POLE_WT$POLE_signatures-delta_POLE_MUT_STL))

delta_POLE_MUT_SPA <- median(DF_SPA_POLE_MUT$POLE_signatures)
delta_POLE_WT_SPA <- median(DF_SPA_POLE_WT$POLE_signatures)

distr_delta_SPA_POLE <- DF_SPA_POLE_MUT$POLE_signatures-delta_POLE_WT_SPA
distr_delta_SPA_POLE <- append(distr_delta_SPA_POLE, abs(DF_SPA_POLE_WT$POLE_signatures-delta_POLE_MUT_SPA))

delta_POLE_MUT_SA <- median(DF_SA_POLE_MUT$POLE_signatures)
delta_POLE_WT_SA <- median(DF_SA_POLE_WT$POLE_signatures)

distr_delta_SA_POLE <- DF_SA_POLE_MUT$POLE_signatures-delta_POLE_WT_SA
distr_delta_SA_POLE <- append(distr_delta_SA_POLE, abs(DF_SA_POLE_WT$POLE_signatures-delta_POLE_MUT_SA))


DF_violin_MP_POLE <- data.frame(abs(distr_delta_MP_POLE), "MutationalPatterns")
colnames(DF_violin_MP_POLE) <- c("Delta", "Tool")
DF_violin_DS_POLE <- data.frame(abs(distr_delta_DS_POLE), "deconstructSigs")
colnames(DF_violin_DS_POLE) <- c("Delta", "Tool")
DF_violin_STL_POLE <- data.frame(abs(distr_delta_STL_POLE), "signature.tools.lib")
colnames(DF_violin_STL_POLE) <- c("Delta", "Tool")
DF_violin_SPA_POLE <- data.frame(abs(distr_delta_SPA_POLE), "SigProfilerAssignment")
colnames(DF_violin_SPA_POLE) <- c("Delta", "Tool")
DF_violin_SA_POLE <- data.frame(abs(distr_delta_SA_POLE), "SignatureAnalyzer")
colnames(DF_violin_SA_POLE) <- c("Delta", "Tool")


DF_violin_CRC_POLE <- rbind(DF_violin_MP_POLE, DF_violin_DS_POLE, DF_violin_STL_POLE, DF_violin_SPA_POLE, DF_violin_SA_POLE)
colnames(DF_violin_CRC_POLE) <- c("Delta POLE", "Tool")

DF_violin_CRC_POLE$Color <- rep("MS", 760)


couples <- list(c("MutationalPatterns", "deconstructSigs"), c("MutationalPatterns", "signature.tools.lib"), c("MutationalPatterns", "SigProfilerAssignment"), c("MutationalPatterns", "SignatureAnalyzer"), c("signature.tools.lib", "deconstructSigs"), c("signature.tools.lib", "SigProfilerAssignment"), c("signature.tools.lib", "SignatureAnalyzer"), c("SigProfilerAssignment", "deconstructSigs"), c("SignatureAnalyzer", "deconstructSigs"), c("SigProfilerAssignment", "SignatureAnalyzer"))
DF_violin_CRC_POLE$Tool <- factor(DF_violin_CRC_POLE$Tool, levels = c("MutationalPatterns", "deconstructSigs", "signature.tools.lib", "SigProfilerAssignment", "SignatureAnalyzer"))


ggplot(DF_violin_CRC_POLE, aes(x=Tool, y=`Delta POLE`, color=Color, fill=Color)) +
  geom_violin(adjust = 8) +
  scale_color_manual(values=c("Dark Blue")) +
  scale_fill_manual(values=c("alice blue")) +
  scale_y_continuous(breaks = seq(0, 1, 0.25),
                     limits = c(0,1),
                     expand = c(0, 0)) +
  force_panelsizes(row=unit(6, "cm"), cols = unit(8, "cm")) +
  theme(plot.title = element_text(size = 10, hjust = 0.5, face = "bold"), 
        axis.title.x = element_text(size = 10), 
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(angle = -45, vjust = 0, hjust = 0.1),
        legend.text = element_text(size = 16), 
        legend.position = "None" ,
        panel.background = element_rect(fill = "White", colour = "Black", linetype = 1, size = 0.8),
        panel.grid.major = element_line(colour = "Black", size = 0.1)) +
  geom_boxplot(width=0.1, colour = "black", fill = "white", outlier.colour = NA)+
  stat_compare_means(method ="wilcox", label = "p.label", ref.group = Tool, comparisons = couples)


#Stat MMR Clinical Dataset


delta_MSI_TCGA_MP <- median(DF_MP_MSI_TCGA$MMR_signatures)
delta_MSS_TCGA_MP <- median(DF_MP_MSS_TCGA$MMR_signatures)

distr_delta_MP_TCGA <- DF_MP_MSI_TCGA$MMR_signatures-delta_MSS_TCGA_MP
distr_delta_MP_TCGA <- append(distr_delta_MP_TCGA, abs(DF_MP_MSS_TCGA$MMR_signatures-delta_MSI_TCGA_MP))

delta_MSI_TCGA_DS <- median(DF_DS_MSI_TCGA$MMR_signatures)
delta_MSS_TCGA_DS <- median(DF_DS_MSS_TCGA$MMR_signatures)

distr_delta_DS_TCGA <- DF_DS_MSI_TCGA$MMR_signatures-delta_MSS_TCGA_DS
distr_delta_DS_TCGA <- append(distr_delta_DS_TCGA, abs(DF_DS_MSS_TCGA$MMR_signatures-delta_MSI_TCGA_DS))

delta_MSI_TCGA_STL <- median(DF_STL_MSI_TCGA$MMR_signatures)
delta_MSS_TCGA_STL <- median(DF_STL_MSS_TCGA$MMR_signatures)

distr_delta_STL_TCGA <- DF_STL_MSI_TCGA$MMR_signatures-delta_MSS_TCGA_STL
distr_delta_STL_TCGA <- append(distr_delta_STL_TCGA, abs(DF_STL_MSS_TCGA$MMR_signatures-delta_MSI_TCGA_STL))

delta_MSI_TCGA_SPA <- median(DF_SPA_MSI_TCGA$MMR_signatures)
delta_MSS_TCGA_SPA <- median(DF_SPA_MSS_TCGA$MMR_signatures)

distr_delta_SPA_TCGA <- DF_SPA_MSI_TCGA$MMR_signatures-delta_MSS_TCGA_SPA
distr_delta_SPA_TCGA <- append(distr_delta_SPA_TCGA, abs(DF_SPA_MSS_TCGA$MMR_signatures-delta_MSI_TCGA_SPA))

delta_MSI_TCGA_SA <- median(DF_SA_MSI_TCGA$MMR_signatures)
delta_MSS_TCGA_SA <- median(DF_SA_MSS_TCGA$MMR_signatures)

distr_delta_SA_TCGA <- DF_SA_MSI_TCGA$MMR_signatures-delta_MSS_TCGA_SA
distr_delta_SA_TCGA <- append(distr_delta_SA_TCGA, abs(DF_SA_MSS_TCGA$MMR_signatures-delta_MSI_TCGA_SA))



DF_violin_MP_TCGA <- data.frame(matrix(nrow = 132, ncol = 2))
DF_violin_MP_TCGA$X1 <- abs(distr_delta_MP_TCGA)
DF_violin_MP_TCGA$X2 <- "MutationalPatterns"

DF_violin_DS_TCGA <- data.frame(matrix(nrow = 132, ncol = 2))
DF_violin_DS_TCGA$X1 <- abs(distr_delta_DS_TCGA)
DF_violin_DS_TCGA$X2 <- "deconstructSigs"

DF_violin_STL_TCGA <- data.frame(matrix(nrow = 132, ncol = 2))
DF_violin_STL_TCGA$X1 <- abs(distr_delta_STL_TCGA)
DF_violin_STL_TCGA$X2 <- "signature.tools.lib"

DF_violin_SPA_TCGA <- data.frame(matrix(nrow = 132, ncol = 2))
DF_violin_SPA_TCGA$X1 <- abs(distr_delta_SPA_TCGA)
DF_violin_SPA_TCGA$X2 <- "SigProfilerAssignment"

DF_violin_SA_TCGA <- data.frame(matrix(nrow = 132, ncol = 2))
DF_violin_SA_TCGA$X1 <- abs(distr_delta_SA_TCGA)
DF_violin_SA_TCGA$X2 <- "SignatureAnalyzer"


DF_violin_TCGA <- rbind(DF_violin_MP_TCGA, DF_violin_DS_TCGA, DF_violin_STL_TCGA, DF_violin_SPA_TCGA, DF_violin_SA_TCGA)
colnames(DF_violin_TCGA) <- c("Delta MMR", "Tool")

couples <- list(c("MutationalPatterns", "deconstructSigs"), c("MutationalPatterns", "signature.tools.lib"), c("MutationalPatterns", "SigProfilerAssignment"), c("MutationalPatterns", "SignatureAnalyzer"), c("signature.tools.lib", "deconstructSigs"), c("signature.tools.lib", "SigProfilerAssignment"), c("signature.tools.lib", "SignatureAnalyzer"), c("SigProfilerAssignment", "deconstructSigs"), c("SignatureAnalyzer", "deconstructSigs"), c("SigProfilerAssignment", "SignatureAnalyzer"))
DF_violin_TCGA$Tool <- factor(DF_violin_TCGA$Tool, levels = c("MutationalPatterns", "deconstructSigs", "signature.tools.lib", "SigProfilerAssignment", "SignatureAnalyzer"))

DF_violin_TCGA$Color <- rep("MS", 660)

ggplot(DF_violin_TCGA, aes(x=Tool, y=`Delta MMR`, fill=Color)) +
  geom_violin(adjust = 4) +
  scale_color_manual(values=c("brown3")) +
  scale_fill_manual(values=c("coral")) +
  scale_y_continuous(breaks = seq(0, 1, 0.25),
                     limits = c(0,1),
                     expand = c(0, 0)) +
  force_panelsizes(row=unit(6, "cm"), cols = unit(8, "cm")) +
  theme(plot.title = element_text(size = 10, hjust = 0.5, face = "bold"), 
        axis.title.x = element_text(size = 10), 
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(angle = -45, vjust = 0, hjust = 0.1),
        legend.text = element_text(size = 16), 
        legend.position = "None" ,
        panel.background = element_rect(fill = "White", colour = "Black", linetype = 1, size = 0.8),
        panel.grid.major = element_line(colour = "Black", size = 0.1)) +
  geom_boxplot(width=0.1, colour = "black", fill = "white", outlier.colour = NA)+
  stat_compare_means(method ="wilcox", label = "p.label", ref.group = Tool, comparisons = couples)


#Stat POLE Clinical Dataset

delta_POLE_MUT_TCGA_MP <- median(DF_MP_POLE_MUT_TCGA$POLE_signatures)
delta_POLE_WT_TCGA_MP <- median(DF_MP_POLE_WT_TCGA$POLE_signatures)

distr_delta_MP_TCGA_POLE <- DF_MP_POLE_MUT_TCGA$POLE_signatures-delta_POLE_WT_TCGA_MP
distr_delta_MP_TCGA_POLE <- append(distr_delta_MP_TCGA_POLE, abs(DF_MP_POLE_WT_TCGA$POLE_signatures-delta_POLE_MUT_TCGA_MP))

delta_POLE_MUT_TCGA_DS <- median(DF_DS_POLE_MUT_TCGA$POLE_signatures)
delta_POLE_WT_TCGA_DS <- median(DF_DS_POLE_WT_TCGA$POLE_signatures)

distr_delta_DS_TCGA_POLE <- DF_DS_POLE_MUT_TCGA$POLE_signatures-delta_POLE_WT_TCGA_DS
distr_delta_DS_TCGA_POLE <- append(distr_delta_DS_TCGA_POLE, abs(DF_DS_POLE_WT_TCGA$POLE_signatures-delta_POLE_MUT_TCGA_DS))

delta_POLE_MUT_TCGA_STL <- median(DF_STL_POLE_MUT_TCGA$POLE_signatures)
delta_POLE_WT_TCGA_STL <- median(DF_STL_POLE_WT_TCGA$POLE_signatures)

distr_delta_STL_TCGA_POLE <- DF_STL_POLE_MUT_TCGA$POLE_signatures-delta_POLE_WT_TCGA_STL
distr_delta_STL_TCGA_POLE <- append(distr_delta_STL_TCGA_POLE, abs(DF_STL_POLE_WT_TCGA$POLE_signatures-delta_POLE_MUT_TCGA_STL))

delta_POLE_MUT_TCGA_SPA <- median(DF_SPA_POLE_MUT_TCGA$POLE_signatures)
delta_POLE_WT_TCGA_SPA <- median(DF_SPA_POLE_WT_TCGA$POLE_signatures)

distr_delta_SPA_TCGA_POLE <- DF_SPA_POLE_MUT_TCGA$POLE_signatures-delta_POLE_WT_TCGA_SPA
distr_delta_SPA_TCGA_POLE <- append(distr_delta_SPA_TCGA_POLE, abs(DF_SPA_POLE_WT_TCGA$POLE_signatures-delta_POLE_MUT_TCGA_SPA))

delta_POLE_MUT_TCGA_SA <- median(DF_SA_POLE_MUT_TCGA$POLE_signatures)
delta_POLE_WT_TCGA_SA <- median(DF_SA_POLE_WT_TCGA$POLE_signatures)

distr_delta_SA_TCGA_POLE <- DF_SA_POLE_MUT_TCGA$POLE_signatures-delta_POLE_WT_TCGA_SA
distr_delta_SA_TCGA_POLE <- append(distr_delta_SA_TCGA_POLE, abs(DF_SA_POLE_WT_TCGA$POLE_signatures-delta_POLE_MUT_TCGA_SA))


DF_violin_MP_POLE_TCGA <- data.frame(abs(distr_delta_MP_TCGA_POLE), "MutationalPatterns")
colnames(DF_violin_MP_POLE_TCGA) <- c("Delta", "Tool")
DF_violin_DS_POLE_TCGA <- data.frame(abs(distr_delta_DS_TCGA_POLE), "deconstructSigs")
colnames(DF_violin_DS_POLE_TCGA) <- c("Delta", "Tool")
DF_violin_STL_POLE_TCGA <- data.frame(abs(distr_delta_STL_TCGA_POLE), "signature.tools.lib")
colnames(DF_violin_STL_POLE_TCGA) <- c("Delta", "Tool")
DF_violin_SPA_POLE_TCGA <- data.frame(abs(distr_delta_SPA_TCGA_POLE), "SigProfilerAssignment")
colnames(DF_violin_SPA_POLE_TCGA) <- c("Delta", "Tool")
DF_violin_SA_POLE_TCGA <- data.frame(abs(distr_delta_SA_TCGA_POLE), "SignatureAnalyzer")
colnames(DF_violin_SA_POLE_TCGA) <- c("Delta", "Tool")


DF_violin_TCGA_POLE <- rbind(DF_violin_MP_POLE_TCGA, DF_violin_DS_POLE_TCGA, DF_violin_STL_POLE_TCGA, DF_violin_SPA_POLE_TCGA, DF_violin_SA_POLE_TCGA)
colnames(DF_violin_TCGA_POLE) <- c("Delta POLE", "Tool")

DF_violin_TCGA_POLE$Color <- rep("MS", 320)


couples <- list(c("MutationalPatterns", "deconstructSigs"), c("MutationalPatterns", "signature.tools.lib"), c("MutationalPatterns", "SigProfilerAssignment"), c("MutationalPatterns", "SignatureAnalyzer"), c("signature.tools.lib", "deconstructSigs"), c("signature.tools.lib", "SigProfilerAssignment"), c("signature.tools.lib", "SignatureAnalyzer"), c("SigProfilerAssignment", "deconstructSigs"), c("SignatureAnalyzer", "deconstructSigs"), c("SigProfilerAssignment", "SignatureAnalyzer"))
DF_violin_TCGA_POLE$Tool <- factor(DF_violin_TCGA_POLE$Tool, levels = c("MutationalPatterns", "deconstructSigs", "signature.tools.lib", "SigProfilerAssignment", "SignatureAnalyzer"))


ggplot(DF_violin_TCGA_POLE, aes(x=Tool, y=`Delta POLE`, color=Color, fill=Color)) +
  geom_violin(adjust = 8) +
  scale_color_manual(values=c("Dark Blue")) +
  scale_fill_manual(values=c("alice blue")) +
  scale_y_continuous(breaks = seq(0, 1, 0.25),
                     limits = c(0,1),
                     expand = c(0, 0)) +
  force_panelsizes(row=unit(6, "cm"), cols = unit(8, "cm")) +
  theme(plot.title = element_text(size = 10, hjust = 0.5, face = "bold"), 
        axis.title.x = element_text(size = 10), 
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(angle = -45, vjust = 0, hjust = 0.1),
        legend.text = element_text(size = 16), 
        legend.position = "None" ,
        panel.background = element_rect(fill = "White", colour = "Black", linetype = 1, size = 0.8),
        panel.grid.major = element_line(colour = "Black", size = 0.1)) +
  geom_boxplot(width=0.1, colour = "black", fill = "white", outlier.colour = NA)+
  stat_compare_means(method ="wilcox", label = "p.label", ref.group = Tool, comparisons = couples)
