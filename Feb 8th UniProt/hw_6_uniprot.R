library(Biostrings)
library(UniprotR)
library(protti)

hasbroucki <- readDNAStringSet("C:/Users/emmag/OneDrive/School Work/2024-2025/SPRING 2024/ST GR Bioinformatics - BSC 6932/Bioinformatics_Spring2024/Megascops asio cytb MSA/m_asio_hasbroucki_cytb.fasta")
Biostrings::translate(hasbroucki)
m_asio_prot <- Biostrings::translate(hasbroucki)
writeXStringSet(m_asio_prot, "hasbroucki_prot.fasta", format="fasta")
read.csv("feb_8th_accessions.txt")
access_list <- read.table("feb_8th_accessions.txt")
access_list$V1
access_string <- access_list$V1
GetProteinGOInfo(access_string)
