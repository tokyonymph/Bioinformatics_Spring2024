# Open required packages to create a multiple sequence alignment (MSA) in R.
library(Biostrings)
library(msa)
library(seqinr)
library(ape)

# We will first do the MSA for Megascops asio. Start by creating vectors for each subspecies.
m_a_asio <- readDNAStringSet("C:/Users/emmag/OneDrive/School Work/2024-2025/SPRING 2024/ST GR Bioinformatics - BSC 6932/Bioinformatics_Spring2024/Independent Project/Screech Owl Cytb Data/Megascops asio/m_asio_asio_cytb.fasta")
m_a_hasbroucki <- readDNAStringSet("C:/Users/emmag/OneDrive/School Work/2024-2025/SPRING 2024/ST GR Bioinformatics - BSC 6932/Bioinformatics_Spring2024/Independent Project/Screech Owl Cytb Data/Megascops asio/m_asio_hasbroucki_cytb.fasta")
m_a_mccallii <- readDNAStringSet("C:/Users/emmag/OneDrive/School Work/2024-2025/SPRING 2024/ST GR Bioinformatics - BSC 6932/Bioinformatics_Spring2024/Independent Project/Screech Owl Cytb Data/Megascops asio/m_asio_mccallii_cytb.fasta")
m_a_naevius <- readDNAStringSet("C:/Users/emmag/OneDrive/School Work/2024-2025/SPRING 2024/ST GR Bioinformatics - BSC 6932/Bioinformatics_Spring2024/Independent Project/Screech Owl Cytb Data/Megascops asio/m_asio_naevius_cytb.fasta")

# Create vector with all four sequences at once.
m_asio_seqs <- c(m_a_asio, m_a_hasbroucki, m_a_mccallii, m_a_naevius)

# Generate multiple sequence alignment with new vector and export to fasta file.
m_asio_msa <- msa(m_asio_seqs, method=c("ClustalW"))
alignment_m_asio <- msaConvert(fm_asio_msa, type="phangorn::phyDat")
write.phyDat(alignment_m_asio, "meg_asio_msa.fasta", format="fasta")

