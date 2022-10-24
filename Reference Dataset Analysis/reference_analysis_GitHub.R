#!/usr/bin/env Rscript

print("================> Reference Dataset Analysis <================")

setwd("mutational_signatures")

library(MutationalPatterns)
library(ggplot2)
library(ggh4x)
library(ggpubr)

print("--> Reference Loading...")

COSMIC2 <- read.table(paste(getwd(), "COSMIC_v2_SBS_GRCh38.txt", sep = "/") , header = TRUE, row.names = 1)
COSMIC2 <- as.matrix(sapply(COSMIC2, as.numeric))

COSMIC3.2 <- read.table(paste(getwd(), "COSMIC_v3.2_SBS_GRCh38.txt", sep = "/") , header = TRUE, row.names = 1)
COSMIC3.2 <- as.matrix(sapply(COSMIC3.2, as.numeric))

CRC_Specific <- c("SBS1", "SBS2", "SBS3", "SBS5", "SBS6", "SBS9", "SBS10a", "SBS10b", "SBS12", "SBS13", "SBS15", "SBS17a", "SBS17b", "SBS18", "SBS20", "SBS21", "SBS26", "SBS28", "SBS30", "SBS37", "SBS39", "SBS40", "SBS41", "SBS44")
Dataset_CRC_Tissue_Specific <- COSMIC3.2[,CRC_Specific]

CRC_Cell_Bank <- read.table(paste(getwd(), "CRCBank.SBS96.all", sep = "/"), header = T, row.names = 1)
TCGA_Dataset <-  read.table(paste(getwd(), "TCGA.SBS96.all", sep = "/"), header = T, row.names = 1)
colnames(TCGA_Dataset) <- gsub('\\.', "", substr(colnames(TCGA_Dataset), 6, 12))
names(TCGA_Dataset)[names(TCGA_Dataset) == '5MAAT6'] <- 'X5MAAT6'
CRC_Cell_Bank_TSO_500 <- read.table(paste(getwd(), "TCGA_TSO.SBS96.all", sep = "/"), header = T, row.names = 1)
CRC_WGS <- read.table(paste(getwd(), "CRC_WGS.SBS96.all", sep = "/"), header = T, row.names = 1)
TCGA_TSO500 <- read.table(paste(getwd(), "TCGA_TSO.SBS96.all", sep = "/"), header = T, row.names = 1)

#MSI-MSS
MSI_CRC_Bank <- read.table(paste(getwd(), "MSI_Cell_Lines_230", sep = "/"))
MSI_CRC_Bank <- as.vector(MSI_CRC_Bank$V1)

MSS_CRC_Bank <- read.table(paste(getwd(), "MSS_Cell_Lines_230", sep = "/"))
MSS_CRC_Bank <- as.vector(MSS_CRC_Bank$V1)

TCGA_MSS <- c("A62671","A62674","A62677","A62681","A63810","A64107","AA3495","AA3506","AA3510","AA3530","AA3666","AA3673","AA3675","AA3678","AA3679","AA3681","AA3684","AA3685","AA3693","AA3695","AA3696","AA3812","AA3814","AA3818","AA3831","AA3837","AA3841","AA3846","AA3848","AA3850","AA3851","AA3867","AA3875","AA3956","AA3968","AA3971","AA3975","AA3976","AA3977","AA3980","AA3984","AA3986","AA3989","AA3994","AAA01I","AAA01T","AAA01V","AAA01X","AAA01Z","AAA02H","AAA02K","AAA02O","AAA02Y","AAA03F","AAA03J","AY4071","AZ4308","AZ4315","AZ4682","CA5256","CM4744","CM4748","CM4752","CM5341")
	      
TCGA_MSI <- c("X5MAAT6", "A62672", "A62686", "A63809", "A65661", "A65665", "A66653", "A6A565", "AA3492", "AA3663", "AA3672", "AA3713", "AA3715", "AA3811", "AA3815", "AA3821", "AA3833", "AA3845", "AA3864", "AA3877", "AA3947", "AA3949", "AA3950", "AA3966", "AAA01P", "AAA01R", "AAA022", "AAA02R", "AD5900", "AD6889", "AD6895", "AD6964", "ADA5EJ", "AM5821", "AU6004", "AY6197", "AZ4313", "AZ4615", "AZ6598", "CK4951", "CK5913", "CK5916", "CK6746", "CM4743", "CM4746", "CM5861", "CM6162", "CM6171", "CM6674", "D56530", "D56540", "D56927", "D56928", "D56930", "DMA1HB","F46570", "F46703", "F46856", "G46302", "G46304", "G46309", "G46320", "G46586", "G46588", "G46628", "NHA5IV", "QGA5Z2", "WSAB45")

#POLE
MSS_POLEmut_CRC_Bank  <- read.table(paste(getwd(), "POLEmut_MSS_Cell_Lines", sep = "/"))
MSS_POLEmut_CRC_Bank  <- as.vector(MSS_POLEmut_CRC_Bank$V1)

MSS_POLEwt_CRC_Bank <- setdiff(MSS_CRC_Bank, MSS_POLEmut_CRC_Bank)

MSS_POLE_Mut_TCGA <- c("AA3510", "AA3678", "AA3977", "AA3984", "AZ4315")  
MSS_POLE_WT_TCGA <- setdiff(TCGA_MSS, MSS_POLE_Mut_TCGA)

#MMR Signatures
MMR_signature_cosmic2 <- c("Signature_6", "Signature_14", "Signature_15", "Signature_20", "Signature_21", "Signature_26")
MMR_signatures <- c("SBS6", "SBS14", "SBS15", "SBS20", "SBS21", "SBS26", "SBS44")
MMR_signatures_TS <- c("SBS6", "SBS15", "SBS20", "SBS21", "SBS26", "SBS44")

#Analisi di fitting con reference diversi
print("--> Mutational Signatures Analysis preclinical dataset")

Fitting_CRC_Cosmic2_MP <- fit_to_signatures(CRC_Cell_Bank, COSMIC2)
CRC_COSMIC2_MP <- Fitting_CRC_Cosmic2_MP$contribution
CRC_COSMIC2_MP_Norm <- CRC_COSMIC2_MP/colSums(CRC_COSMIC2_MP)[col(CRC_COSMIC2_MP)]
CRC_COSMIC2_MP_Norm <- as.data.frame(CRC_COSMIC2_MP_Norm)

Fitting_CRC_Cosmic3.2_MP <- fit_to_signatures(CRC_Cell_Bank, COSMIC3.2)
CRC_COSMIC3.2_MP <- Fitting_CRC_Cosmic3.2_MP$contribution
CRC_COSMIC3.2_MP_Norm <- CRC_COSMIC3.2_MP/colSums(CRC_COSMIC3.2_MP)[col(CRC_COSMIC3.2_MP)]
CRC_COSMIC3.2_MP_Norm <- as.data.frame(CRC_COSMIC3.2_MP_Norm)

Fitting_CRC_Tissue_Specific_MP <- fit_to_signatures(CRC_Cell_Bank, Dataset_CRC_Tissue_Specific)
CRC_Tissue_Specific_MP <- Fitting_CRC_Tissue_Specific_MP$contribution
CRC_Tissue_Specific_MP_Norm <- CRC_Tissue_Specific_MP/colSums(CRC_Tissue_Specific_MP)[col(CRC_Tissue_Specific_MP)]
CRC_Tissue_Specific_MP_Norm <- as.data.frame(CRC_Tissue_Specific_MP_Norm)


print("--> Mutational Signatures Analysis clinical dataset")

Fitting_TCGA_Cosmic2_MP <- fit_to_signatures(TCGA_Dataset, COSMIC2)
TCGA_COSMIC2_MP <- Fitting_TCGA_Cosmic2_MP$contribution
TCGA_COSMIC2_MP_Norm <- TCGA_COSMIC2_MP/colSums(TCGA_COSMIC2_MP)[col(TCGA_COSMIC2_MP)]
TCGA_COSMIC2_MP_Norm <- as.data.frame(TCGA_COSMIC2_MP_Norm)

Fitting_TCGA_Cosmic3.2_MP <- fit_to_signatures(TCGA_Dataset, COSMIC3.2)
TCGA_COSMIC3.2_MP <- Fitting_TCGA_Cosmic3.2_MP$contribution
TCGA_COSMIC3.2_MP_Norm <- TCGA_COSMIC3.2_MP/colSums(TCGA_COSMIC3.2_MP)[col(TCGA_COSMIC3.2_MP)]
TCGA_COSMIC3.2_MP_Norm <- as.data.frame(TCGA_COSMIC3.2_MP_Norm)

Fitting_TCGA_Tissue_Specific_MP <- fit_to_signatures(TCGA_Dataset, Dataset_CRC_Tissue_Specific)
TCGA_Tissue_Specific_MP <- Fitting_TCGA_Tissue_Specific_MP$contribution
TCGA_Tissue_Specific_MP_Norm <- TCGA_Tissue_Specific_MP/colSums(TCGA_Tissue_Specific_MP)[col(TCGA_Tissue_Specific_MP)]
TCGA_Tissue_Specific_MP_Norm <- as.data.frame(TCGA_Tissue_Specific_MP_Norm)

print("--> Cosine similarity")

cos_sim_CRC_COSMIC2_MP <- diag(cos_sim_matrix(CRC_Cell_Bank, Fitting_CRC_Cosmic2_MP$reconstructed))
cos_sim_TCGA_COSMIC2_MP <- diag(cos_sim_matrix(TCGA_Dataset, Fitting_TCGA_Cosmic2_MP$reconstructed))

cos_sim_CRC_COSMIC3.2_MP <- diag(cos_sim_matrix(CRC_Cell_Bank, Fitting_CRC_Cosmic3.2_MP$reconstructed))
cos_sim_TCGA_COSMIC3.2_MP <- diag(cos_sim_matrix(TCGA_Dataset, Fitting_TCGA_Cosmic3.2_MP$reconstructed))

cos_sim_CRC_Tissue_Specific_MP <- diag(cos_sim_matrix(CRC_Cell_Bank, Fitting_CRC_Tissue_Specific_MP$reconstructed))
cos_sim_TCGA_Tissue_Specific_MP <- diag(cos_sim_matrix(TCGA_Dataset, Fitting_TCGA_Tissue_Specific_MP$reconstructed))


df_cos_sim_CRC_COSMIC2_MP <- data.frame(cos_sim_CRC_COSMIC2_MP, "Cosmic2", "Preclinical Dataset")
colnames(df_cos_sim_CRC_COSMIC2_MP) <- c("Cosine_Similarity", "Reference_Dataset", "Dataset")

df_cos_sim_CRC_COSMIC3.2_MP <- data.frame(cos_sim_CRC_COSMIC3.2_MP, "Cosmic3.2", "Preclinical Dataset")
colnames(df_cos_sim_CRC_COSMIC3.2_MP) <- c("Cosine_Similarity", "Reference_Dataset", "Dataset")

df_cos_sim_CRC_Tissue_Specific_MP <- data.frame(cos_sim_CRC_Tissue_Specific_MP, "Tissue Specific", "Preclinical Dataset")
colnames(df_cos_sim_CRC_Tissue_Specific_MP) <- c("Cosine_Similarity", "Reference_Dataset", "Dataset")


df_cos_sim_TCGA_COSMIC2_MP <- data.frame(cos_sim_TCGA_COSMIC2_MP, "Cosmic2", "Clinical Dataset")
colnames(df_cos_sim_TCGA_COSMIC2_MP) <- c("Cosine_Similarity", "Reference_Dataset", "Dataset")

df_cos_sim_TCGA_COSMIC3.2_MP <- data.frame(cos_sim_TCGA_COSMIC3.2_MP, "Cosmic3.2", "Clinical Dataset")
colnames(df_cos_sim_TCGA_COSMIC3.2_MP) <- c("Cosine_Similarity", "Reference_Dataset", "Dataset")

df_cos_sim_TCGA_Tissue_Specific_MP <- data.frame(cos_sim_TCGA_Tissue_Specific_MP, "Tissue Specific", "Clinical Dataset")
colnames(df_cos_sim_TCGA_Tissue_Specific_MP) <- c("Cosine_Similarity", "Reference_Dataset", "Dataset")


df_Figure_3A <- rbind(df_cos_sim_CRC_COSMIC2_MP, df_cos_sim_CRC_COSMIC3.2_MP, df_cos_sim_CRC_Tissue_Specific_MP,
		      df_cos_sim_TCGA_COSMIC2_MP, df_cos_sim_TCGA_COSMIC3.2_MP, df_cos_sim_TCGA_Tissue_Specific_MP)

print("--> Figure 3A")

ggplot(df_Figure_3A, aes(x=Dataset, y=Cosine_Similarity, color=Reference_Dataset)) + 
	geom_boxplot(lwd=0.9) + 
	scale_fill_manual(aesthetics = "color", values=c("#1B8148", "Orange", "#3852BD","#1B8148", "Orange", "#3852BD")) +
	labs(y="Cosine Similarity",
	     color="Reference Dataset")+
	theme_bw() +
	theme(axis.title.y = element_blank(), 
	      legend.text = element_text(size = 10), 
	      legend.title = element_text(size = 10), 
	      axis.text = element_text(size = 10)) +
	scale_y_continuous(breaks = seq(0.825, 1, 0.05),
			   limits = c(0.825,1),
			   expand = c(0, 0)) +
	coord_flip() +
	force_panelsizes(row=unit(6, "cm"), cols = unit(8, "cm")) +
	geom_hline(yintercept=0.9, linetype="dashed")

ggsave("Figure_3A.svg")

print("--> MMR Analysis")

#Preclinical dataset
MSI_COSMIC2_MP <- as.data.frame(apply(CRC_COSMIC2_MP_Norm[MMR_signature_cosmic2,MSI_CRC_Bank], 2, sum))
MSI_COSMIC3.2_MP <- as.data.frame(apply(CRC_COSMIC3.2_MP_Norm[MMR_signatures,MSI_CRC_Bank], 2, sum))
MSI_Tissue_Specific_MP <- as.data.frame(apply(CRC_Tissue_Specific_MP_Norm[MMR_signatures_TS, MSI_CRC_Bank], 2, sum))

MSS_COSMIC2_MP <- as.data.frame(apply(CRC_COSMIC2_MP_Norm[MMR_signature_cosmic2,MSS_CRC_Bank], 2, sum))
MSS_COSMIC3.2_MP <- as.data.frame(apply(CRC_COSMIC3.2_MP_Norm[MMR_signatures,MSS_CRC_Bank], 2, sum))
MSS_Tissue_Specific_MP <- as.data.frame(apply(CRC_Tissue_Specific_MP_Norm[MMR_signatures_TS, MSS_CRC_Bank], 2, sum))

#Clinical dataset

MSI_TCGA_COSMIC2_MP <- as.data.frame(apply(TCGA_COSMIC2_MP_Norm[MMR_signature_cosmic2,TCGA_MSI], 2, sum))
MSI_TCGA_COSMIC3.2_MP <- as.data.frame(apply(TCGA_COSMIC3.2_MP_Norm[MMR_signatures,TCGA_MSI], 2, sum))
MSI_TCGA_Tissue_Specific_MP <- as.data.frame(apply(TCGA_Tissue_Specific_MP_Norm[MMR_signatures_TS, TCGA_MSI], 2, sum))

MSS_TCGA_COSMIC2_MP <- as.data.frame(apply(TCGA_COSMIC2_MP_Norm[MMR_signature_cosmic2,TCGA_MSS], 2, sum))
MSS_TCGA_COSMIC3.2_MP <- as.data.frame(apply(TCGA_COSMIC3.2_MP_Norm[MMR_signatures,TCGA_MSS], 2, sum))
MSS_TCGA_Tissue_Specific_MP <- as.data.frame(apply(TCGA_Tissue_Specific_MP_Norm[MMR_signatures_TS, TCGA_MSS], 2, sum))


df_MSS_COSMIC2_CRC <- data.frame(MSS_COSMIC2_MP, "Cosmic2", "MSS")
colnames(df_MSS_COSMIC2_CRC) <- c("MMR_Contribution", "Dataset", "MS_Status")

df_MSS_COSMIC3.2_CRC <- data.frame(MSS_COSMIC3.2_MP, "Cosmic3.2", "MSS")
colnames(df_MSS_COSMIC3.2_CRC) <- c("MMR_Contribution", "Dataset", "MS_Status")

df_MSS_Tissue_Specific_CRC <- data.frame(MSS_Tissue_Specific_MP, "Tissue_Specific", "MSS")
colnames(df_MSS_Tissue_Specific_CRC) <- c("MMR_Contribution", "Dataset", "MS_Status")

df_MSI_COSMIC2_CRC <- data.frame(MSI_COSMIC2_MP, "Cosmic2", "MSI")
colnames(df_MSI_COSMIC2_CRC) <- c("MMR_Contribution", "Dataset", "MS_Status")

df_MSI_COSMIC3.2_CRC <- data.frame(MSI_COSMIC3.2_MP, "Cosmic3.2", "MSI")
colnames(df_MSI_COSMIC3.2_CRC) <- c("MMR_Contribution", "Dataset", "MS_Status")

df_MSI_Tissue_Specific_CRC <- data.frame(MSI_Tissue_Specific_MP, "Tissue_Specific", "MSI")
colnames(df_MSI_Tissue_Specific_CRC) <- c("MMR_Contribution", "Dataset", "MS_Status")

df_Figure_3B <- rbind(df_MSS_COSMIC2_CRC, df_MSS_COSMIC3.2_CRC, df_MSS_Tissue_Specific_CRC, 
		      df_MSI_COSMIC2_CRC, df_MSI_COSMIC3.2_CRC, df_MSI_Tissue_Specific_CRC)


print("--> Figure 3B")

ggplot(df_Figure_3B, aes(x=Dataset, y=MMR_Contribution, fill=MS_Status)) +
	geom_boxplot(aes(fill=MS_Status)) +
	scale_fill_manual(values=c("lightblue1", "salmon2")) +
	labs(x = "Reference Dataset",
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

	
ggsave("Figure_3B.svg")


df_MSS_TCGA_COSMIC2_CRC <- data.frame(MSS_TCGA_COSMIC2_MP, "Cosmic2", "MSS")
colnames(df_MSS_TCGA_COSMIC2_CRC) <- c("MMR_Contribution", "Dataset", "MS_Status")

df_MSS_TCGA_COSMIC3.2_CRC <- data.frame(MSS_TCGA_COSMIC3.2_MP, "Cosmic3.2", "MSS")
colnames(df_MSS_TCGA_COSMIC3.2_CRC) <- c("MMR_Contribution", "Dataset", "MS_Status")

df_MSS_TCGA_Tissue_Specific_CRC <- data.frame(MSS_TCGA_Tissue_Specific_MP, "Tissue_Specific", "MSS")
colnames(df_MSS_TCGA_Tissue_Specific_CRC) <- c("MMR_Contribution", "Dataset", "MS_Status")

df_MSI_TCGA_COSMIC2_CRC <- data.frame(MSI_TCGA_COSMIC2_MP, "Cosmic2", "MSI")
colnames(df_MSI_TCGA_COSMIC2_CRC) <- c("MMR_Contribution", "Dataset", "MS_Status")

df_MSI_TCGA_COSMIC3.2_CRC <- data.frame(MSI_TCGA_COSMIC3.2_MP, "Cosmic3.2", "MSI")
colnames(df_MSI_TCGA_COSMIC3.2_CRC) <- c("MMR_Contribution", "Dataset", "MS_Status")

df_MSI_TCGA_Tissue_Specific_CRC <- data.frame(MSI_TCGA_Tissue_Specific_MP, "Tissue_Specific", "MSI")
colnames(df_MSI_TCGA_Tissue_Specific_CRC) <- c("MMR_Contribution", "Dataset", "MS_Status")

df_Figure_3C <- rbind(df_MSS_TCGA_COSMIC2_CRC, df_MSS_TCGA_COSMIC3.2_CRC, df_MSS_TCGA_Tissue_Specific_CRC, df_MSI_TCGA_COSMIC2_CRC, df_MSI_TCGA_COSMIC3.2_CRC, df_MSI_TCGA_Tissue_Specific_CRC )


print("--> Figure 3C")

ggplot(df_Figure_3C, aes(x=Dataset, y=MMR_Contribution, fill=MS_Status)) +
        geom_boxplot(aes(fill=MS_Status)) +
        scale_fill_manual(values=c("lightblue1", "salmon2")) +
        labs(x = "Reference Dataset",
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

	
ggsave("Figure_3C.svg")

print("--> POLE Analysis")

#Preclinical dataset
POLE_MUT_COSMIC2_MP <- as.data.frame(apply(CRC_COSMIC2_MP_Norm["Signature_10", MSS_POLEmut_CRC_Bank], 2, sum))
POLE_MUT_COSMIC3.2_MP <- as.data.frame(apply(CRC_COSMIC3.2_MP_Norm[c("SBS10a", "SBS10b"), MSS_POLEmut_CRC_Bank], 2, sum))
POLE_MUT_Tissue_Specific_MP <- as.data.frame(apply(CRC_Tissue_Specific_MP_Norm[c("SBS10a", "SBS10b"), MSS_POLEmut_CRC_Bank], 2, sum))

POLE_WT_COSMIC2_MP <- as.data.frame(apply(CRC_COSMIC2_MP_Norm["Signature_10", MSS_POLEwt_CRC_Bank], 2, sum))
POLE_WT_COSMIC3.2_MP <- as.data.frame(apply(CRC_COSMIC3.2_MP_Norm[c("SBS10a", "SBS10b"), MSS_POLEwt_CRC_Bank], 2, sum))
POLE_WT_Tissue_Specific_MP <- as.data.frame(apply(CRC_Tissue_Specific_MP_Norm[c("SBS10a", "SBS10b"), MSS_POLEwt_CRC_Bank], 2, sum))

#Clinical dataset

POLE_MUT_TCGA_COSMIC2_MP <- as.data.frame(apply(TCGA_COSMIC2_MP_Norm["Signature_10", MSS_POLE_Mut_TCGA], 2, sum))
POLE_MUT_TCGA_COSMIC3.2_MP <- as.data.frame(apply(TCGA_COSMIC3.2_MP_Norm[c("SBS10a", "SBS10b"), MSS_POLE_Mut_TCGA], 2, sum))
POLE_MUT_TCGA_Tissue_Specific_MP <- as.data.frame(apply(TCGA_Tissue_Specific_MP_Norm[c("SBS10a", "SBS10b"), MSS_POLE_Mut_TCGA], 2, sum))

POLE_WT_TCGA_COSMIC2_MP <- as.data.frame(apply(TCGA_COSMIC2_MP_Norm["Signature_10", MSS_POLE_WT_TCGA], 2, sum))
POLE_WT_TCGA_COSMIC3.2_MP <- as.data.frame(apply(TCGA_COSMIC3.2_MP_Norm[c("SBS10a", "SBS10b"), MSS_POLE_WT_TCGA], 2, sum))
POLE_WT_TCGA_Tissue_Specific_MP <- as.data.frame(apply(TCGA_Tissue_Specific_MP_Norm[c("SBS10a", "SBS10b"), MSS_POLE_WT_TCGA], 2, sum))



df_POLE_MUT_COSMIC2_CRC <- data.frame(POLE_MUT_COSMIC2_MP, "Cosmic2", "POLE Mutated")
colnames(df_POLE_MUT_COSMIC2_CRC) <- c("POLE_Contribution", "Dataset", "POLE_Status")

df_POLE_MUT_COSMIC3.2_CRC <- data.frame(POLE_MUT_COSMIC3.2_MP, "Cosmic3.2", "POLE Mutated")
colnames(df_POLE_MUT_COSMIC3.2_CRC) <- c("POLE_Contribution", "Dataset", "POLE_Status")

df_POLE_MUT_Tissue_Specific_CRC <- data.frame(POLE_MUT_COSMIC3.2_MP, "Tissue_Specific", "POLE Mutated")
colnames(df_POLE_MUT_Tissue_Specific_CRC) <- c("POLE_Contribution", "Dataset", "POLE_Status")


df_POLE_WT_COSMIC2_CRC <- data.frame(POLE_WT_COSMIC2_MP, "Cosmic2", "POLE Wild Type")
colnames(df_POLE_WT_COSMIC2_CRC) <- c("POLE_Contribution", "Dataset", "POLE_Status")

df_POLE_WT_COSMIC3.2_CRC <- data.frame(POLE_WT_COSMIC3.2_MP, "Cosmic3.2", "POLE Wild Type")
colnames(df_POLE_WT_COSMIC3.2_CRC) <- c("POLE_Contribution", "Dataset", "POLE_Status")

df_POLE_WT_Tissue_Specific_CRC <- data.frame(POLE_WT_Tissue_Specific_MP, "Tissue_Specific", "POLE Wild Type")
colnames(df_POLE_WT_Tissue_Specific_CRC) <- c("POLE_Contribution", "Dataset", "POLE_Status")

df_Figure_3D <- rbind(df_POLE_MUT_COSMIC2_CRC, df_POLE_MUT_COSMIC3.2_CRC, df_POLE_MUT_Tissue_Specific_CRC,
		      df_POLE_WT_COSMIC2_CRC, df_POLE_WT_COSMIC3.2_CRC, df_POLE_WT_Tissue_Specific_CRC)


print("--> Figure 3D")

ggplot(df_Figure_3D, aes(x=Dataset, y=POLE_Contribution, fill=POLE_Status))+
	geom_boxplot(aes(fill=POLE_Status)) +
       	scale_fill_manual(values=c("darkgoldenrod1", "olivedrab4")) + 
	labs(x = "Reference Dataset",		
	     y = "Overall POLE Signatures Contribution",
	     fill = "POLE Status") +
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



ggsave("Figure_3D.svg")

df_POLE_MUT_TCGA_COSMIC2_CRC <- data.frame(POLE_MUT_TCGA_COSMIC2_MP, "Cosmic2", "POLE Mutated")
colnames(df_POLE_MUT_TCGA_COSMIC2_CRC) <- c("POLE_Contribution", "Dataset", "POLE_Status")

df_POLE_MUT_TCGA_COSMIC3.2_CRC <- data.frame(POLE_MUT_TCGA_COSMIC3.2_MP, "Cosmic3.2", "POLE Mutated")
colnames(df_POLE_MUT_TCGA_COSMIC3.2_CRC) <- c("POLE_Contribution", "Dataset", "POLE_Status")

df_POLE_MUT_TCGA_Tissue_Specific_CRC <- data.frame(POLE_MUT_TCGA_COSMIC3.2_MP, "Tissue_Specific", "POLE Mutated")
colnames(df_POLE_MUT_TCGA_Tissue_Specific_CRC) <- c("POLE_Contribution", "Dataset", "POLE_Status")


df_POLE_WT_TCGA_COSMIC2_CRC <- data.frame(POLE_WT_TCGA_COSMIC2_MP, "Cosmic2", "POLE Wild Type")
colnames(df_POLE_WT_TCGA_COSMIC2_CRC) <- c("POLE_Contribution", "Dataset", "POLE_Status")

df_POLE_WT_TCGA_COSMIC3.2_CRC <- data.frame(POLE_WT_TCGA_COSMIC3.2_MP, "Cosmic3.2", "POLE Wild Type")
colnames(df_POLE_WT_TCGA_COSMIC3.2_CRC) <- c("POLE_Contribution", "Dataset", "POLE_Status")

df_POLE_WT_TCGA_Tissue_Specific_CRC <- data.frame(POLE_WT_TCGA_Tissue_Specific_MP, "Tissue_Specific", "POLE Wild Type")
colnames(df_POLE_WT_TCGA_Tissue_Specific_CRC) <- c("POLE_Contribution", "Dataset", "POLE_Status")

df_Figure_3E <- rbind(df_POLE_MUT_TCGA_COSMIC2_CRC, df_POLE_MUT_TCGA_COSMIC3.2_CRC, df_POLE_MUT_TCGA_Tissue_Specific_CRC,
		      df_POLE_WT_TCGA_COSMIC2_CRC, df_POLE_WT_TCGA_COSMIC3.2_CRC, df_POLE_WT_TCGA_Tissue_Specific_CRC)


print("--> Figure 3E")

ggplot(df_Figure_3E, aes(x=Dataset, y=POLE_Contribution, fill=POLE_Status))+
        geom_boxplot(aes(fill=POLE_Status)) +
        scale_fill_manual(values=c("darkgoldenrod1", "olivedrab4")) +
	labs(x = "Reference Dataset",
	     y = "Overall POLE Signatures Contribution",
	     fill = "POLE Status") +
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
	      element_line(color = "black", size = 0,4)) 
	      
ggsave("Figure_3E.svg")

print("MMR Contribution")

Cosmic2_MSI <- apply(CRC_COSMIC2_MP_Norm[, MSI_CRC_Bank], 1, median)
Cosmic3.2_MSI <- apply(CRC_COSMIC3.2_MP_Norm[, MSI_CRC_Bank], 1, median)
Tissue_Specific_MSI <- apply(CRC_Tissue_Specific_MP_Norm[, MSI_CRC_Bank], 1, median)

Cosmic2_MSI_Norm <- Cosmic2_MSI/sum(Cosmic2_MSI)
Cosmic3.2_MSI_Norm <- Cosmic3.2_MSI/sum(Cosmic3.2_MSI)
Tissue_Specific_MSI_Norm <- Tissue_Specific_MSI/sum(Tissue_Specific_MSI)

#COSMIC2

MSI_MMR_Cosmic2 <- Cosmic2_MSI_Norm[MMR_signature_cosmic2]
MSI_Other_Cosmic2 <- Cosmic2_MSI_Norm[setdiff(names(Cosmic2_MSI), names(MSI_MMR_Cosmic2))]

Cosmic2_MSI_MMR_df <- data.frame(c(MSI_MMR_Cosmic2, sum(MSI_Other_Cosmic2)), c("SBS6", "SBS14", "SBS15", "SBS20", "SBS21", "SBS26", "Other_SBS"), c(rep("Cosmic2", 7)))
colnames(Cosmic2_MSI_MMR_df) <- c("Contribution", "Signature", "Reference_Dataset")

#COSMIC3.2

MSI_MMR_Cosmic3.2 <- Cosmic3.2_MSI_Norm[MMR_signatures]
MSI_Other_Cosmic3.2 <- Cosmic3.2_MSI_Norm[setdiff(names(Cosmic3.2_MSI), names(MSI_MMR_Cosmic3.2))]

Cosmic3.2_MSI_MMR_df <- data.frame(c(MSI_MMR_Cosmic3.2, sum(MSI_Other_Cosmic3.2)), c("SBS6", "SBS14", "SBS15", "SBS20", "SBS21", "SBS26", "SBS44", "Other_SBS"), c(rep("Cosmic3.2", 8)))
colnames(Cosmic3.2_MSI_MMR_df) <- c("Contribution", "Signature", "Reference_Dataset")

#TISSUE SPECIFIC

MSI_MMR_TS <- Tissue_Specific_MSI_Norm[MMR_signatures_TS]
MSI_Other_TS <- Tissue_Specific_MSI_Norm[setdiff(names(Tissue_Specific_MSI), names(MSI_MMR_TS))]

Tissue_Specific_MSI_MMR_df <- data.frame(c(MSI_MMR_TS, sum(MSI_Other_TS)), c("SBS6", "SBS15", "SBS20", "SBS21", "SBS26", "SBS44", "Other_SBS"), c(rep("Tissue Specific", 7)))
colnames(Tissue_Specific_MSI_MMR_df) <- c("Contribution", "Signature", "Reference_Dataset")

MMR_contribution_MSI_reference <- rbind(Cosmic2_MSI_MMR_df, Cosmic3.2_MSI_MMR_df, Tissue_Specific_MSI_MMR_df)


ggplot(MMR_contribution_MSI_reference, aes(x=Reference_Dataset, y=Contribution, fill=Signature)) +
	geom_bar(stat="identity", width = 0.7) +
	scale_y_continuous(breaks = seq(0, 1, 0.25),
			   limits = c(0,1), 
			   expand = c(0, 0)) +
	force_panelsizes(row=unit(6, "cm"), cols = unit(8, "cm")) +
        scale_fill_manual(values=c("gray80", "seagreen1", "royalblue3", "orchid3", "turquoise2", "goldenrod1", "springgreen4", "orangered2")) + 
	labs(x= "Reference Dataset",
	     y= "MMR Mutational Signatures Contribution",
	     fill="MMR Mutational Signatures") +
	theme_bw() +
	theme(axis.title.x = element_text(size = 10),
	      axis.title.y = element_text(size = 10),
	      legend.text = element_text(size = 10),
	      legend.title = element_text(size = 10),
	      legend.position = "top",
	      panel.background = element_rect(fill="white"),
	      axis.text = element_text(size = 10),
	      element_line(color = "black", size = 0.4))

	ggsave("Figure_3F.svg")


print("FINE")
