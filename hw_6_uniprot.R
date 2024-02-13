library(Biostrings)
library(UniprotR)
library(protti)
library(r3dmol)

hasbroucki <- readDNAStringSet("C:/Users/emmag/OneDrive/School Work/2024-2025/SPRING 2024/ST GR Bioinformatics - BSC 6932/Bioinformatics_Spring2024/Megascops asio cytb MSA/m_asio_hasbroucki_cytb.fasta")
Biostrings::translate(hasbroucki)
m_asio_prot <- Biostrings::translate(hasbroucki)
writeXStringSet(m_asio_prot, "hasbroucki_prot.fasta", format="fasta")
read.csv("feb_8th_accessions.txt")
access_list <- read.table("feb_8th_accessions.txt")
access_list$V1
access_string <- access_list$V1
GetProteinGOInfo(access_string)
prot_functions <- GetProteinGOInfo(access_string)
PlotGoInfo(prot_functions) 
PlotGOAll(GOObj = prot_functions, Top = 10, directorypath = getwd(), width = 8, height = 5)
GetPathology_Biotech(access_string)
Get.diseases(access_list, directorypath = NULL)
fetch_uniprot(access_list)
uniprot_table <- fetch_uniprot(access_list)
prot_pdbs <- read.table("feb_8th_pdbs.txt")
prot_structure_table <- fetch_pdb(prot_pdbs)
prot_structure_table$pdb_ids
?fetch_alphafold_prediction()
fetch_alphafold_prediction
alphafold <- fetch_alphafold_prediction(uniprot_ids = c("P0A799", "P08839"), return_data_frame = TRUE)
prot_r3d <- r3dmol(alphafold)
r3dmol::m_add_models(id = prot_r3d, data = alphafold, format = c("pdb"))
