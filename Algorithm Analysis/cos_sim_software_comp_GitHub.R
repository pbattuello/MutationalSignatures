#!/usr/bin/env Rscript

print("================> Bioinformatic Tools Analysis <================")

setwd("mutational_signatures")

library(ggplot2)
library(ggpubr)
library(ggh4x)

MP_CRC <- read.table("cos_sim_CRC_COSMIC2_MP.txt", sep="\t")
MP_TCGA <- read.table("cos_sim_TCGA_COSMIC2_MP.txt", sep="\t")
SA_CRC <- read.table("cos_sim_CRC_COSMIC2_SA.txt", sep="\t")
SA_TCGA <- read.table("cos_sim_TCGA_COSMIC2_SA.txt", sep="\t")
SPA_CRC <- read.table("cos_sim_CRC_COSMIC2_SPA.txt", sep="\t")
SPA_TCGA <- read.table("cos_sim_TCGA_COSMIC2_SPA.txt", sep="\t")
DS_CRC <- read.table("cos_sim_CRC_COSMIC2_DS.txt", sep="\t")
DS_TCGA <- read.table("cos_sim_TCGA_COSMIC2_DS.txt", sep="\t")
STL_CRC <- read.table("cos_sim_CRC_COSMIC2_STL.txt", sep="\t")
STL_TCGA <- read.table("cos_sim_TCGA_COSMIC2_STL.txt", sep="\t")

TCGA_name <- rep("CRC Clinical Dataset", 152)
CRC_name <- rep("CRC Preclinical Dataset", 230)
MP_CRC_name <- rep("MutationalPatterns", 230)
MP_TCGA_name <- rep("MutationalPatterns", 152)
SA_CRC_name <- rep("SignatureAnalyzer", 230)
SA_TCGA_name <- rep("SignatureAnalyzer", 152)
STL_CRC_name <- rep("Signature.tool.lib", 230)
STL_TCGA_name <- rep("Signature.tool.lib", 152)
DS_CRC_name <- rep("DeconstructSigs", 230)
DS_TCGA_name <- rep("DeconstructSigs", 152)
SPA_CRC_name <- rep("SigProfilerAssignment", 230)
SPA_TCGA_name <- rep("SigProfilerAssignment", 152)

#ggplot format

MP_CRC$V1 <- CRC_name
MP_TCGA$V1 <- TCGA_name
SA_CRC$V1 <- CRC_name
SA_TCGA$V1 <- TCGA_name
SPA_CRC$V1 <- CRC_name
SPA_TCGA$V1 <- TCGA_name
DS_CRC$V1 <- CRC_name
DS_TCGA$V1 <- TCGA_name
STL_CRC$V1 <- CRC_name
STL_TCGA$V1 <- TCGA_name

MP_CRC$V3 <- MP_CRC_name
MP_TCGA$V3 <- MP_TCGA_name
SA_CRC$V3 <- SA_CRC_name
SA_TCGA$V3 <- SA_TCGA_name
SPA_CRC$V3 <- SPA_CRC_name
SPA_TCGA$V3 <- SPA_TCGA_name
DS_CRC$V3 <- DS_CRC_name
DS_TCGA$V3 <- DS_TCGA_name
STL_CRC$V3 <- STL_CRC_name
STL_TCGA$V3 <- STL_TCGA_name

#unisco i dataset

cosine_sim_df <- rbind(MP_CRC, MP_TCGA, SA_CRC, SA_TCGA, STL_CRC, STL_TCGA, SPA_CRC, SPA_TCGA, DS_CRC, DS_TCGA)
colnames(cosine_sim_df) <- c("Dataset", "Cosine_Similarity", "Tool")

#grafico con ggplot

cosine_sim_df$Tool <- factor(cosine_sim_df$Tool, levels = c("MutationalPatterns", "DeconstructSigs", "Signature.tool.lib", "SigProfilerAssignment", "SignatureAnalyzer"))
cosine_sim_df$Dataset <- factor(cosine_sim_df$Dataset, levels = c("CRC Preclinical Dataset", "CRC Clinical Dataset"))

print("Figure 2A")

ggplot(cosine_sim_df, aes(x=Dataset, Cosine_Similarity, y=Cosine_Similarity, fill=Tool)) +
	geom_boxplot(aes(fill=Tool)) +
	labs(title = "Cosine Similarity Comparison", x= "Dataset", y="Cosine Similarity") +
      	scale_y_continuous(breaks = seq(0, 1, 0.25),
				 limits = c(0.75,1),
	    			 expand = c(0, 0))+
      	force_panelsizes(row=unit(10, "cm"), cols = unit(14, "cm")) +
	theme_bw() +
      	theme(axis.title.x = element_text(size = 10), 	
		    axis.title.y = element_text(size = 10), legend.text = element_text(size = 10), 
		    legend.title = element_text(size = 10), 
	    	    axis.text = element_text(size = 10),
	    	    panel.background = element_rect(fill = "White", colour = "Black", linetype = 1, size = 0.4),
	    	    panel.grid.major = element_line(colour = "Black", size = 0.1)) +
	geom_hline(yintercept=0.9, linetype="dashed")

ggsave("Cosine_Similarity_Tools.svg")

print("END")


