library(seqinr)

rm(list=ls())
WGS<-read.csv("Cosmo_WGS/Outputs/Cosmo_WGS_sequence_data.csv")
N<-read.csv("Cosmo_N_all/Outputs/new_designation.csv")
alignment<-read.alignment("Cosmo_N_all/Alignment/Cosmo_N_all_aligned.fasta", format = "fasta")

clades<-c("AF1a", "AF1b", "AF1c", "AF4", "AM1", "AM2a", "AM2b", "AM3a", "AM3b", "AM4",
          "CA1", "CA2", "CA3", "CE", "EE", "ME1a", "ME1b", "ME2", "NEE", "WE")

for (j in 1:length(clades)) {
  ref<-N[grep(clades[j], N$lineage),]
  for (i in 1:length(ref$ID)) {
    N$cluster[which(N$ID %in% ref$ID)]<-ref$lineage[i]
    data<-N[grep(clades[j], N$lineage),]
    data<-data[c(1,8)]
    names(data)<-c("ID", "cluster")
  }
  
  dir.create(clades[j])
  dir.create(paste(clades[j], "/reference", sep = ""))
  write.csv(data, file = paste(clades[j], "/updated_reference_clusters.csv", sep = ""), row.names = F)
  
  test<-which(alignment$nam %in% N$ID[grep(clades[j], N$cluster)])
  sequences<-alignment$seq[test]
  names<-alignment$nam[test]
  
  write.fasta(sequences, names, file.out = paste(clades[j], "/all_N.fasta", sep = ""))
  
  ref_align<-which(alignment$nam %in% N$ID[grep(clades[j], N$lineage)])
  sequences<-alignment$seq[ref_align]
  names<-alignment$nam[ref_align]
  
  write.fasta(sequences, names, file.out = paste(clades[j], "/reference/reference_aligned.fasta", sep = ""))
}




