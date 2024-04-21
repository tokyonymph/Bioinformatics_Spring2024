# Open required packages to create a multiple sequence alignment (MSA) in R.
library(Biostrings)
library(msa)
library(seqinr)
library(ape)
library(phangorn)

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

# Repeat same process for other species; next, Megascops kennicottii.
m_k_aikeni <- readDNAStringSet("C:/Users/emmag/OneDrive/School Work/2024-2025/SPRING 2024/ST GR Bioinformatics - BSC 6932/Bioinformatics_Spring2024/Independent Project/Screech Owl Cytb Data/Megascops kennicottii/m kennicottii aikeni.fasta")
m_k_benderi <- readDNAStringSet("C:/Users/emmag/OneDrive/School Work/2024-2025/SPRING 2024/ST GR Bioinformatics - BSC 6932/Bioinformatics_Spring2024/Independent Project/Screech Owl Cytb Data/Megascops kennicottii/m kennicottii benderi.fasta")
m_k_kennicottii <- readDNAStringSet("C:/Users/emmag/OneDrive/School Work/2024-2025/SPRING 2024/ST GR Bioinformatics - BSC 6932/Bioinformatics_Spring2024/Independent Project/Screech Owl Cytb Data/Megascops kennicottii/m kennicottii kennicottii.fasta")
m_k_mcfarlandii <- readDNAStringSet("C:/Users/emmag/OneDrive/School Work/2024-2025/SPRING 2024/ST GR Bioinformatics - BSC 6932/Bioinformatics_Spring2024/Independent Project/Screech Owl Cytb Data/Megascops kennicottii/m kennicottii mcfarlandii.fasta")
m_k_suttoni <- readDNAStringSet("C:/Users/emmag/OneDrive/School Work/2024-2025/SPRING 2024/ST GR Bioinformatics - BSC 6932/Bioinformatics_Spring2024/Independent Project/Screech Owl Cytb Data/Megascops kennicottii/m kennicottii suttoni.fasta")

m_kenn_seqs <- c(m_k_aikeni, m_k_benderi, m_k_kennicottii, m_k_mcfarlandii, m_k_suttoni)
m_kenn_msa <- msa(m_kenn_seqs, method=c("ClustalW"))
align_m_kenn <- msaConvert(m_kenn_msa, type="phangorn::phyDat")
write.phyDat(align_m_kenn, "meg_kenn_msa.fasta", format="fasta")

m_wc_set_path <- "C:/Users/emmag/OneDrive/School Work/2024-2025/SPRING 2024/ST GR Bioinformatics - BSC 6932/Bioinformatics_Spring2024/Independent Project/Screech Owl Cytb Data/Megascops watsoni + choliba/Megascops_watsonii_and_choliba_cytb.fas"
m_wach_set <- readDNAStringSet(m_wc_set_path, format = "fasta")
m_wach_msa <- msa(m_wach_set, method = c("ClustalW"))
m_wach_align <- msaConvert(m_wach_msa, type="phangorn::phyDat")
write.phyDat(m_wach_align, "m_wat_cho_msa.fasta", format="fasta")
