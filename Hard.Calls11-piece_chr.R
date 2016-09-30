###########################
#HARD CALLS
###########################
library(foreach)
library(parallel)
library(doParallel)
library(doMC)
library(data.table)
library(Rcpp)
#setwd("C:/Users/Lilith Moss/Documents/MATH/RESEARCH/BMA_Project/#Real_Data/GEN2-Code/Imputed Data/")

chr = 1
script = 18

#Cluster Paths
dimension_path <- "/home/pmd-01/chemenya/CHS/txtDosage/dimensions/"
output_path <- "/home/pmd-01/chemenya/CHS/Split_Imputed_Results/"
dosage_path <- "/home/pmd-01/chemenya/CHS/Split_Imputed/"

#Read in the number of SNPs to be read in the dosage file
dimension <- as.numeric(read.table(paste0(dimension_path,"chr",chr,".row"))[1])

#Read in how many files there are for this chromosome
num.files <- ceiling(dimension/5000)

#Which files to read
scripts <- ceiling(num.files/20)
read.these <- c((20*(script-1)):(20+20*(script-1)-1))
files <- sprintf(".%0004d", read.these)

#Define dominant G function
cppFunction('NumericVector personloop2(NumericVector d) {
  int n = d.size();
            NumericVector G(n/3);
            for(int i=0; i<2999; ++i){
              if (d[(3*i)]>=0.9){
                G[i]=0;
              } else if (d[1+(3*i)]>=0.9){
                G[i]=1;
              } else if (d[2+(3*i)]>=0.9){
                G[i]=1;
              } else {
                G[i]=-9;
              }            
            }
            return G;     
}')

for(j in 1:length(files)){
  file <- paste0(dosage_path,"chr", chr,files[j])
  if(file.exists(file)){
    dos <- read.table(paste0(dosage_path,"chr", chr,files[j]))
      result <- do.call(rbind,lapply(1:nrow(dos), function(j){
        D <- dos[j,]
        d <- as.vector(unlist(D[,4:length(D)]))
        personloop2(d)
      }) ) 
      trunc <- cbind(dos[,1:3],result) #Put together with SNP name and info
      saveRDS(trunc, paste0(output_path,"G.chr", chr,files[j]))
  } else 
    print(paste0("No Such File:",file))
}