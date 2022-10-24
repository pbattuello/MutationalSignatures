#!/usr/bin/env Rscript

library(MutationalPatterns)

CRC_Cell_Bank <- read.table("CRCBank.SBS96.txt" , header = T, row.names = 1)
TCGA <- read.table("TCGA.SBS96.txt" , header = T, row.names = 1)

COSMIC2 <- read.table("COSMIC_v2_SBS_GRCh38.txt", header = T, row.names = 1)
COSMIC2 <- as.matrix(sapply(COSMIC2, as.numeric))

Fitting_CRC_Cosmic2 <- read.table("Fitting_CRCBank_Cosmic2.txt" , header = T, row.names = 1)
Fitting_TCGA_Cosmic2 <- read.table("Fitting_TCGA_Cosmic2.txt" , header = T, row.names = 1)

#Reconstructed CRC

reconstructed_Cosmic2 <- data.frame(matrix(NA, nrow = length(seq(1:96)), ncol = length(seq(1:230))))
row.names(reconstructed_Cosmic2) <- rownames(CRC_Cell_Bank)
colnames(reconstructed_Cosmic2) <- colnames(CRC_Cell_Bank)

results_values_Cosmic2 <- data.frame(matrix(NA, nrow = length(seq(1:96)), ncol = length(seq(1:30))))
row.names(reconstructed_Cosmic2) <- rownames(CRC_Cell_Bank)
colnames(reconstructed_Cosmic2) <- colnames(COSMIC2)

for (j in 1:230){
	for (i in 1:30){
		results_values_Cosmic2[,i] <- Fitting_CRC_Cosmic2[i,j]*(COSMIC2[,i])
        }
        reconstructed_Cosmic2[,j] <- apply(results_values_Cosmic2, 1, sum)
}

#Cosine Similarity CRC

cos_sim_CRC_COSMIC2 <- diag(cos_sim_matrix(CRC_Cell_Bank, reconstructed_Cosmic2))
names(cos_sim_CRC_COSMIC2) <- colnames(CRC_Cell_Bank)
write.table(cos_sim_CRC_COSMIC2, "cos_sim_CRC_COSMIC2.txt", sep = "\t", quote = F)

#Reconstructed CRC

reconstructed_Cosmic2 <- data.frame(matrix(NA, nrow = length(seq(1:96)), ncol = length(seq(1:152))))
row.names(reconstructed_Cosmic2) <- rownames(TCGA)
colnames(reconstructed_Cosmic2) <- colnames(TCGA)

results_values_Cosmic2 <- data.frame(matrix(NA, nrow = length(seq(1:96)), ncol = length(seq(1:30))))
row.names(reconstructed_Cosmic2) <- rownames(TCGA)
colnames(reconstructed_Cosmic2) <- colnames(COSMIC2)

for (j in 1:230){
	for (i in 1:30){
	results_values_Cosmic2[,i] <- Fitting_TCGA_Cosmic2[i,j]*(COSMIC2[,i])
        }
        reconstructed_Cosmic2[,j] <- apply(results_values_Cosmic2, 1, sum)
}


#Cosine Similarity CRC

cos_sim_TCGA_COSMIC2 <- diag(cos_sim_matrix(TCGA, reconstructed_Cosmic2))
names(cos_sim_TCGA_COSMIC2) <- colnames(TCGA)
write.table(cos_sim_TCGA_COSMIC2, "cos_sim_TCGA_COSMIC2.txt", sep = "\t", quote = F)


