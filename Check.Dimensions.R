#Count the number of scripts per chromosome for proliferation purposes
#Cluster Paths
dimension_path <- "/home/pmd-01/chemenya/CHS/txtDosage/dimensions/"
output_path <- "/home/pmd-01/chemenya/CHS/Split_Imputed_Results/"
dosage_path <- "/home/pmd-01/chemenya/CHS/Split_Imputed/"

#Loop through all 22 chromosomes
files <- do.call(rbind,lapply(1:22,function(i){
  
  #Set Chromosome
  chr=i
  
  #Read in the number of SNPs to be read in the dosage file
  dimension <- as.numeric(read.table(paste0(dimension_path,"chr",chr,".row"))[1])
  
  #Which files represent the chromosome?
  results <- list.files(output_path,paste0("chr",chr,".0"))
  results.dim <- lapply(1:length(results),function(x){
    file.dim <- dim(readRDS(results[x]))[1]
  })
  dimension.sum <- sum(as.numeric(results.dim))
  
  #Put all together to save
  cbind(chr,dimension,dimension.sum)

}))
write.table(files,"Results.dimensions.txt")
files