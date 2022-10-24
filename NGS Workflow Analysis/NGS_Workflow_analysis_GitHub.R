#!/usr/bin/env Rscript

print("================> NGS Workflow Analysis <================")

setwd("mutational_signatures")

library(MutationalPatterns)
library(ggplot2)
library(ggh4x)
library(ggpubr)

print("--> Reference Datasets Loading...")

COSMIC3.2 <- read.table(paste(getwd(), "COSMIC_v3.2_SBS_GRCh38.txt", sep = "/") , header = TRUE, row.names = 1)
COSMIC3.2 <- as.matrix(sapply(COSMIC3.2, as.numeric))

CRC_Cell_Bank <- read.table(paste(getwd(), "CRCBank.SBS96.all", sep = "/"), header = T, row.names = 1)
CRC_WGS <- read.table(paste(getwd(), "CRC_WGS.SBS96.all", sep = "/"), header = T, row.names = 1)
CRC_Cell_Bank_TSO_500 <- read.table(paste(getwd(), "TSO_500_230.SBS96.all", sep = "/"), header = T, row.names = 1)

MSI <- scan(file=paste(getwd(), "MSI_Cell_Lines_230", sep = "/"), what = "character")
MSS <- scan(file=paste(getwd(), "MSS_Cell_Lines_230", sep = "/"), what = "character")
MMR_signatures <- c("SBS6", "SBS14", "SBS15", "SBS20", "SBS21", "SBS26", "SBS44")
POLE_signatures <- c("SBS10a", "SBS10b")
Artefacts_signatures <- c("SBS27", "SBS43", "SBS45", "SBS46", "SBS47", "SBS48", 
			  "SBS49", "SBS50", "SBS51", "SBS52", "SBS53", "SBS54", 
			  "SBS55", "SBS56", "SBS57", "SBS58", "SBS59", "SBS60")
Flat_signatures <- c("SBS3", "SBS5", "SBS8", "SBS40", "SBS89")
POLE <- scan(file=paste(getwd(), "POLEmut_MSS_Cell_Lines", sep = "/"), what = "character")
MSS_POLE_WT <- setdiff(MSS, POLE)

print("--> Mutational Signatures Analysis")

print("WES")
Fitting_CRC_Cosmic3.2_MP <- fit_to_signatures(CRC_Cell_Bank, COSMIC3.2)
CRC_COSMIC3.2_MP <- Fitting_CRC_Cosmic3.2_MP$contribution
CRC_COSMIC3.2_MP_Norm <- CRC_COSMIC3.2_MP/colSums(CRC_COSMIC3.2_MP)[col(CRC_COSMIC3.2_MP)]
CRC_COSMIC3.2_MP_Norm <- as.data.frame(CRC_COSMIC3.2_MP_Norm)

print("WGS")
Fitting_CRC_WGS_Cosmic3.2_MP <- fit_to_signatures(CRC_WGS, COSMIC3.2)
CRC_WGS_COSMIC3.2_MP <- Fitting_CRC_WGS_Cosmic3.2_MP$contribution
CRC_WGS_COSMIC3.2_MP_Norm <- CRC_WGS_COSMIC3.2_MP/colSums(CRC_WGS_COSMIC3.2_MP)[col(CRC_WGS_COSMIC3.2_MP)]
CRC_WGS_COSMIC3.2_MP_Norm <- as.data.frame(CRC_WGS_COSMIC3.2_MP_Norm)

print("Pan-Cancer Panel")
Fitting_CRC_TSO_Cosmic3.2_MP <- fit_to_signatures(CRC_Cell_Bank_TSO_500, COSMIC3.2)
CRC_TSO_COSMIC3.2_MP <- Fitting_CRC_TSO_Cosmic3.2_MP$contribution
CRC_TSO_COSMIC3.2_MP_Norm <- CRC_TSO_COSMIC3.2_MP/colSums(CRC_TSO_COSMIC3.2_MP)[col(CRC_TSO_COSMIC3.2_MP)]
CRC_TSO_COSMIC3.2_MP_Norm <- as.data.frame(CRC_TSO_COSMIC3.2_MP_Norm)

print("--> Cosine similarity Analysis")

cosine_sim_CRC_Cosmic3.2_MP <- diag(cos_sim_matrix(Fitting_CRC_Cosmic3.2_MP$reconstructed, CRC_Cell_Bank))
cosine_sim_CRC_TSO500_Cosmic3.2_MP <- diag(cos_sim_matrix(Fitting_CRC_TSO_Cosmic3.2_MP$reconstructed, CRC_Cell_Bank_TSO_500))
cosine_sim_CRC_WGS_Cosmic3.2_MP <- diag(cos_sim_matrix(Fitting_CRC_WGS_Cosmic3.2_MP$reconstructed, CRC_WGS))

cosine_data_type_WES <- data.frame(cosine_sim_CRC_Cosmic3.2_MP, "WES")
colnames(cosine_data_type_WES) <- c("Cosine_Similarity", "Data_Type")

cosine_data_type_TSO500 <- data.frame(cosine_sim_CRC_TSO500_Cosmic3.2_MP, "TSO-500")
colnames(cosine_data_type_TSO500) <- c("Cosine_Similarity", "Data_Type")

cosine_data_type_WGS <- data.frame(cosine_sim_CRC_WGS_Cosmic3.2_MP, "WGS")
colnames(cosine_data_type_WGS) <- c("Cosine_Similarity", "Data_Type")

cosine_data_type <- rbind(cosine_data_type_WES, cosine_data_type_TSO500, cosine_data_type_WGS)
colnames(cosine_data_type) <- c("Cosine_Similarity", "Data_Type")

print("--> Figure 2A")

ggplot(cosine_data_type, aes(x= factor(Data_Type, level=c("WGS", "WES", "TSO-500")) , y= Cosine_Similarity)) + 
	geom_boxplot(color="black", fill="grey") + 
        coord_flip() +
	labs(x= "NGS Workflow", y="Cosine Similarity") + 
	scale_y_continuous(breaks = seq(0.9, 1, 0.025),
			   limits = c(0.89,1), 
			   expand = c(0, 0))+
	geom_hline(yintercept=0.9, linetype="dashed") +
	force_panelsizes(row=unit(6, "cm"), cols = unit(8, "cm")) +
	theme_bw() +
	theme(axis.title.x = element_text(size = 10), 
	      axis.title.y = element_text(size = 10),
	      axis.text = element_text(size = 10),
	      panel.border = element_rect(color = "black", size = 0.4))

ggsave("Figure_2A.svg")

print("--> MMR Analysis")

MSS_WES <- data.frame(apply(CRC_COSMIC3.2_MP_Norm[MMR_signatures, MSS], 2, sum), "MSS", "WES")
colnames(MSS_WES) <- c("Overall_MMR_Signatures", "MS", "Data_Type")

MSI_WES <- data.frame(apply(CRC_COSMIC3.2_MP_Norm[MMR_signatures, MSI], 2, sum), "MSI", "WES")
colnames(MSI_WES) <- c("Overall_MMR_Signatures", "MS", "Data_Type")

MSS_TSO <- data.frame(apply(CRC_TSO_COSMIC3.2_MP_Norm[MMR_signatures, MSS], 2, sum), "MSS", "TSO-500")
colnames(MSS_TSO) <- c("Overall_MMR_Signatures", "MS", "Data_Type")

MSI_TSO <- data.frame(apply(CRC_TSO_COSMIC3.2_MP_Norm[MMR_signatures, MSI], 2, sum), "MSI", "TSO-500")
colnames(MSI_TSO) <- c("Overall_MMR_Signatures", "MS", "Data_Type")

MSS_WGS <- data.frame(apply(CRC_WGS_COSMIC3.2_MP_Norm[MMR_signatures, intersect(colnames(CRC_WGS_COSMIC3.2_MP_Norm),MSS)], 2, sum), "MSS", "WGS")
colnames(MSS_WGS) <- c("Overall_MMR_Signatures", "MS", "Data_Type")

MSI_WGS <- data.frame(apply(CRC_WGS_COSMIC3.2_MP_Norm[MMR_signatures, intersect(colnames(CRC_WGS_COSMIC3.2_MP_Norm),MSI)], 2, sum), "MSI", "WGS")
colnames(MSI_WGS) <- c("Overall_MMR_Signatures", "MS", "Data_Type")


MS_data_type <- rbind(MSS_WES, MSI_WES, MSS_TSO, MSI_TSO, MSS_WGS, MSI_WGS)


MS_data_type$MS <- as.factor(MS_data_type$MS)
MS_data_type$Data_Type <- as.factor(MS_data_type$Data_Type)


print("--> Figure 2B")

ggplot(MS_data_type, aes(x=Data_Type, y=Overall_MMR_Signatures, fill=MS)) + 
	geom_boxplot(aes(fill=MS)) +
	scale_fill_manual(values=c("lightblue1", "salmon2")) +
	labs(x = "NGS Workflow", 
	     y = "Overall MMR Signatures Contribution",
	     fill = "Microsatellite Status") +
	scale_y_continuous(breaks = seq(0, 1, 0.25),
			   limits = c(0,1),
			   expand = c(0, 0))+
	stat_compare_means(method ="wilcox.test", label = "p.signif") +
	force_panelsizes(row=unit(6, "cm"), cols = unit(8, "cm")) +
	theme_bw() +
	theme(axis.title.x = element_text(size = 10),
	      axis.title.y = element_text(size = 10),
	      legend.text = element_text(size = 10), 
	      legend.title = element_text(size = 10),
	      legend.position = "top",
	      axis.text = element_text(size = 10),
	      element_line(color = "black", size = 0.4))

ggsave("Figure_2B.svg")

print("--> POLE Analysis")

POLE_WES <- data.frame(apply(CRC_COSMIC3.2_MP_Norm[POLE_signatures, POLE], 2, sum), "POLE_MUT", "WES")
colnames(POLE_WES) <- c("Overall_POLE_Signatures", "POLE", "Data_Type")

POLE_WT_WES <- data.frame(apply(CRC_COSMIC3.2_MP_Norm[POLE_signatures, MSS_POLE_WT], 2, sum), "POLE_WT", "WES")
colnames(POLE_WT_WES) <- c("Overall_POLE_Signatures", "POLE", "Data_Type")

POLE_TSO <- data.frame(apply(CRC_TSO_COSMIC3.2_MP_Norm[POLE_signatures, POLE], 2, sum), "POLE_MUT", "TSO-500")
colnames(POLE_TSO) <- c("Overall_POLE_Signatures", "POLE", "Data_Type")

POLE_WT_TSO <- data.frame(apply(CRC_TSO_COSMIC3.2_MP_Norm[POLE_signatures, MSS_POLE_WT], 2, sum), "POLE_WT", "TSO-500")
colnames(POLE_WT_TSO) <- c("Overall_POLE_Signatures", "POLE", "Data_Type")

POLE_WGS <- data.frame(apply(CRC_WGS_COSMIC3.2_MP_Norm[POLE_signatures, intersect(colnames(CRC_WGS_COSMIC3.2_MP_Norm),POLE)], 2, sum), "POLE_MUT", "WGS")
colnames(POLE_WGS) <- c("Overall_POLE_Signatures", "POLE", "Data_Type")

POLE_WT_WGS <- data.frame(apply(CRC_WGS_COSMIC3.2_MP_Norm[POLE_signatures, intersect(colnames(CRC_WGS_COSMIC3.2_MP_Norm),MSS_POLE_WT)], 2, sum), "POLE_WT", "WGS")
colnames(POLE_WT_WGS) <- c("Overall_POLE_Signatures", "POLE", "Data_Type")


POLE_data_type <- rbind(POLE_WES, POLE_WT_WES, POLE_TSO, POLE_WT_TSO, POLE_WGS, POLE_WT_WGS)

POLE_data_type$POLE <- as.factor(POLE_data_type$POLE)
POLE_data_type$Data_Type <- as.factor(POLE_data_type$Data_Type)

print("--> Figure 2C")

ggplot(POLE_data_type, aes(x=Data_Type, y=Overall_POLE_Signatures, fill=POLE)) +
	geom_boxplot(aes(fill=POLE)) +
	scale_y_continuous(breaks = seq(0, 1, 0.25),
			   limits = c(0,1),
			   expand = c(0, 0)) +
	stat_compare_means(method ="wilcox.test", label = "p.signif") +
	force_panelsizes(row=unit(6, "cm"), cols = unit(8, "cm")) +
	scale_fill_manual(values=c("darkgoldenrod1", "olivedrab4")) +
	labs(x= "NGS Workflow", 
	     y= "Overall POLE Signatures Contribution",
	     fill="POLE status") +
	theme_bw() +
	theme(axis.title.x = element_text(size = 10), 
	      axis.title.y = element_text(size = 10), 
	      legend.text = element_text(size = 10),
	      legend.title = element_text(size = 10),
	      legend.position = "top",
	      axis.text = element_text(size = 10),
	      element_line(color = "black", size = 0.4))  
		
ggsave("Figure_2C.svg")

print("--> Signatures Contribution Analysis")

WES_MSI <- apply(CRC_COSMIC3.2_MP_Norm[, MSI], 1, median)
WGS_MSI <- apply(CRC_WGS_COSMIC3.2_MP_Norm[, intersect(colnames(CRC_WGS_COSMIC3.2_MP_Norm),MSI)], 1, median)
TSO_MSI <- apply(CRC_TSO_COSMIC3.2_MP_Norm[, MSI], 1, median)

Signature_contribution <- data.frame(WES_MSI, WGS_MSI, TSO_MSI)
Signature_contribution <- Signature_contribution/colSums(Signature_contribution)[col(Signature_contribution)]

WES_contributions <- c(sum(Signature_contribution[MMR_signatures, 1]), sum(Signature_contribution[Artefacts_signatures, 1]), sum(Signature_contribution[Flat_signatures, 1]), 1-sum(c(sum(Signature_contribution[MMR_signatures, 1]), sum(Signature_contribution[Artefacts_signatures, 1]), sum(Signature_contribution[Flat_signatures, 1]))))

WGS_contributions <- c(sum(Signature_contribution[MMR_signatures, 2]), sum(Signature_contribution[Artefacts_signatures, 2]), sum(Signature_contribution[Flat_signatures, 2]), 1-sum(c(sum(Signature_contribution[MMR_signatures, 2]), sum(Signature_contribution[Artefacts_signatures, 2]), sum(Signature_contribution[Flat_signatures, 2]))))

TSO_contributions <- c(sum(Signature_contribution[MMR_signatures, 3]), sum(Signature_contribution[Artefacts_signatures, 3]), sum(Signature_contribution[Flat_signatures, 3]), 1-sum(c(sum(Signature_contribution[MMR_signatures, 3]), sum(Signature_contribution[Artefacts_signatures, 3]), sum(Signature_contribution[Flat_signatures, 3]))))

Signature_contribution_WES <- data.frame(WES_contributions, rep("WES", times=length(WES_contributions)), c("MMR signatures", "Artefacts signatures", "Flat signatures", "Other signatures"))
colnames(Signature_contribution_WES) <- c("Signature_Contributions", "NGS_Workflow", "Signatures_type")

Signature_contribution_WGS <- data.frame(WGS_contributions, rep("WGS", times=length(WGS_contributions)), c("MMR signatures", "Artefacts signatures", "Flat signatures", "Other signatures"))
colnames(Signature_contribution_WGS) <- c("Signature_Contributions", "NGS_Workflow", "Signatures_type")

Signature_contribution_TSO <- data.frame(TSO_contributions, rep("TSO-500", times=length(TSO_contributions)), c("MMR signatures", "Artefacts signatures", "Flat signatures", "Other signatures"))
colnames(Signature_contribution_TSO) <- c("Signature_Contributions", "NGS_Workflow", "Signatures_type")

Signature_contribution_figure <- rbind(Signature_contribution_WES, Signature_contribution_WGS, Signature_contribution_TSO)

Signature_contribution_figure$Signatures_type <- factor(Signature_contribution_figure$Signatures_type, levels = c("Other signatures", "Artefacts signatures", "Flat signatures", "MMR signatures"))

print("--> Figure 2D")

ggplot(Signature_contribution_figure, aes(x=NGS_Workflow, y=Signature_Contributions, fill=Signatures_type)) +
	geom_bar(stat="identity", width = 0.7) +
	scale_y_continuous(breaks = seq(0, 1, 0.25),
			   limits = c(0,1), 
			   expand = c(0, 0)) +
	force_panelsizes(row=unit(6, "cm"), cols = unit(8, "cm")) +
	scale_fill_manual(values=c("gray80", "deepskyblue3", "goldenrod1", "firebrick3")) +
	labs(x= "NGS Workflow",
	     y= "Mutational Signatures Contribution",
	     fill="Mutational Signatures Types") +
        theme_bw() +
	theme(axis.title.x = element_text(size = 10),
	      axis.title.y = element_text(size = 10),
	      legend.text = element_text(size = 10),
	      legend.title = element_text(size = 10),
	      legend.position = "top",
	      axis.text = element_text(size = 10),
	      element_line(color = "black", size = 0.4)) +
	guides(fill=guide_legend(ncol=2,nrow=3,byrow=TRUE, title.position = "top"))

ggsave("Figure_2D.svg")

print("DONE")
