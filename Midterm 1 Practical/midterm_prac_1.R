# Load in required libraries for multiple sequence alignment, protein structure prediction, etc.
library(Biostrings)
library(msa)
library(seqinr)
library(phangorn)
library(UniprotR)
library(protti)
library(r3dmol)

# Set working directory for Midterm Practical 1 Folder; 
# I did this first to make the process easier.
setwd("C:/Users/emmag/OneDrive/School Work/2024-2025/SPRING 2024/ST GR Bioinformatics - BSC 6932/Bioinformatics_Spring2024/Midterm 1 Practical")

# QUESTION 1. 
# First, I created a character vector containing the file path to the FASTA File.
# Then I read in the DNA sequences and converted them to an XStringSet Object.
seq_path <- "C:/Users/emmag/OneDrive/School Work/2024-2025/SPRING 2024/ST GR Bioinformatics - BSC 6932/Bioinformatics_Spring2024/Midterm 1 Practical/mid1_sequences.fasta"
seq_set <- readDNAStringSet(seq_path, format = "fasta")
# With my DNAStringSet object now created, I generated my multiple sequence alignment.
# I used the "msa" function in the "msa" package and create another object.
mid_1_msa <- msa(seq_set, method = c("ClustalW"))

# QUESTION 2. 
# My MSA was now generated, but I still needed to compare the samples to each other.
# To do so, I used the print and frequency functions to observe any gaps, mismatches, and mutations.
print(mid_1_msa, show = "complete")
# Below, I limited the number of rows to make my evaluation of the mutations easier.
print(mid_1_msa, show=c("complete"), showNames=TRUE, showConsensus=TRUE, halfNrow = 2)
# I also converted my alphabet frequency to a table for convenience.
mid_1_freq <- alphabetFrequency(mid_1_msa)
# According to the frequency I calculated, there are three sequences that differ from the rest.
# These are sequences Homo_sapiens_4, Homo_sapiens_10, and Homo_sapiens_6. 
# They are listed on the printed MSA as [1], [2], and [20] respectively.
# Every other sequence has 144 As, 149 Cs, 182 Gs, and 167 Cs.
# Homo_sapiens_4 or [1] has 145 As and 148 Cs.
# Homo_sapiens_10 or [2] has 148 Cs, 183 Gs and 168 Ts.
# Homo_sapiens_6 or [20] has 141 As and 184 Gs. It has 641 bp total, unlike the others.
# I noted the following mutations in Homo_sapiens_6 or [20] (note: comparison is made to the consensus sequence):
# - Gap at position 1 instead of A, which is most likely a deletion.
# - A point substitution at position 3 instead of C.
# - C point substitution at position 134 instead of T.
# - T point substitution at position 145 instead of A.
# - G point substitution at position 152 instead of C.
# - C point substitution at position 586 instead of G.
# - G point substitution at position 623 instead of A.
# Homo_sapiens_4 or [1] had this mutation:
# - A point substitution at position 39 instead of C.
# Homo_sapiens_10 or [2] had these mutations:
# - G point substitution at position 39 instead of C.
# - T point substitution at position 45 instead of A.

# QUESTION 3. 
# For this question, I used my initial FASTA sequences file and plugged it into BLAST on the GenBank website, using the blastn version.
# According to my BLAST query, the most likely identity of my gene is:
# "Homo sapiens hbb gene for beta globin, partial cds."
# The accession number of this best match is GenBank LC121775.

# QUESTION 4.
# Based on the information from question 2, the individual that was the most different from the others in my sequence was Homo_sapiens_6 or [20].
# I made a copy of my original FASTA file and edited the copy so that it contained only the Homo_sapiens_6 sequence.
# I then read the sequence in and made a new object.
human_6_path <- "C:/Users/emmag/OneDrive/School Work/2024-2025/SPRING 2024/ST GR Bioinformatics - BSC 6932/Bioinformatics_Spring2024/Midterm 1 Practical/h_sap_6_seq.fasta"
human_6_dna <- readDNAStringSet(human_6_path, format = "fasta")
# I translated my new DNAStringSet sequence to protein using the Biostrings "translate" function.
human_6_prot <- Biostrings::translate(human_6_dna)
print(human_6_prot)
# I then wrote a new FASTA file with my translated protein sequence.
writeXStringSet(human_6_prot, "human_6_protein.fasta", format="fasta")

# QUESTION 5.
# I used the BLAST program on the UniProt database to figure out my protein's match.
# According to my BLAST query, the most likely identity of my protein is:
# "Hemoglobin subunit beta, HBB, Homo sapiens (human)."
# The accession number of this entry on UniProt is A0A0J9YWK4.

# QUESTION 6.
# I began by searching the UniProt database for my accession number, A0A0J9YWK4.
# This led to a UniProtKB entry for the Hemoglobin subunit beta protein in Homo sapiens.
# I clicked on the "Variant Viewer" for the entry and found out that
# HBB is associated with the following diseases:
# Beta thalassemia HBB/LCRB; beta-zero thalassemia; beta-plus thalassemia;
# dominant-beta thalassemia; sickle cell disease; methemeglobinemia, beta-globin type;
# and hemoglobin C disease.
# I did another Blast query on Ensembl with my translated FASTA file to identify which variant my sequece had.
# Based on the results, the best match to my protein sequence was the Ensembl transcript ID ENST00000633227.1.
# Also, according to Ensembl, this transcript is associated with the Ensembl variant rs33930165,
# which on UniProt has a synonym of VAR_002864.
# By reading the variant's entry on UniProt via Expasy, I figured out that
# it is associated with HbC or Hemoglobin C disease, which also confers resistance to severe malaria
# based on two linked publications on PubMed.
# Since the best match to my protein sequence is equivalent to this variant entry on UniProt,
# I conclude that my patient most likely has Hemoglobin C disease, which would explain
# why their DNA sequence differs greatly from the other sequences in my data.

# QUESTION 7.
# I have included a screenshot image of the predicted 3D structure
# of the Hemoglobin subunit beta protein in my "Midterm 1 Practical" folder
# within my GitHub Repository. It has the following file title:
# "Hemoglobin Subunit Beta 3D Structure.jpg".
