#Count the number of scripts per chromosome for proliferation purposes
#Cluster Paths
dimension_path <- "/home/pmd-01/chemenya/CHS/txtDosage/dimensions/"
output_path <- "/home/pmd-01/chemenya/CHS/Split_Imputed_Results/"
dosage_path <- "/home/pmd-01/chemenya/CHS/Split_Imputed/"

#Loop through all 22 chromosomes
first.last <- do.call(rbind,lapply(1:22,function(i){
  
  #Set Chromosome
  chr=i
  
  #Read in the number of SNPs to be read in the dosage file
  dimension <- as.numeric(read.table(paste0(dimension_path,"chr",chr,".row"))[1])
  
  #Get first and last SNPs from Map file
  first.map <- read.table(paste0("/home/pmd-01/chemenya/CHS/txtDosage/chr",chr,".map"),nrow=1)[2] 
  last.map <- read.table(paste0("/home/pmd-01/chemenya/CHS/txtDosage/chr",chr,".map"),skip=(dimension-1),nrow=1)[2]
  
  #Get first and last SNPs from results files
  files <- list.files(output_path,paste0("chr",chr,".0"))
  
  #Read in last file for the chromosome
  results.last <- readRDS(files[length(files)])
  
  #Assign first and last SNP names
  first.results <- as.character( readRDS(files[1])[1,1] )
  last.results <- as.character( results.last[nrow(results.last),1] )
  
  #Put all together to save
  SNPS <- as.data.frame( cbind(chr,dimension,first.map,first.results,last.map,last.results) )
  names(SNPS) <- c("chr","Dim","First_Map","First_Res","Last_Map","Last_Res")
  SNPS
}))
write.table(first.last,"First.Last.Match.txt")
first.last