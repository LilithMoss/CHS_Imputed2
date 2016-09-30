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
  
  #Read in how many files there are for this chromosome
  num.files <- ceiling(dimension/5000)
  
  #Which files to read
  scripts <- ceiling(num.files/20)
  
  #Put all together to save
  cbind(chr,num.files,scripts)

}))
files