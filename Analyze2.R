############################################
# ANALYZE THE RDS. FILES IN DOMINANT FORMAT 
############################################
#Format_2_wCovariates.R
#################################
# Analyze Data in Plink Format
################################
#Version Control
#Made Format_2 for individual chromosomes
################################
#mainDir <- "C:/Users/Lilith Moss/Documents/MATH/RESEARCH/BMA_Project/#Real_Data/GEN2-Code/"
mainDir <- "/auto/pmd-01/chemenya/CHS/"
#source("/home/pmd-01/chemenya/CHS/Real.functions.R") #Source the functions
source(paste0(mainDir,"Real.functions_CapPP.R"))
library(data.table)
library(qqman)
library(parallel)
library(BMA)

chr <- K #Set which chromosome will be processed
script <- J #Set which script this is 
i=chr

#Assign Prior Probabilities
I_ge=0.5;I_g=0.5;I_gxe=0.5 #Set prior probabilities

epi_path <- "/auto/pmd-01/chemenya/CHS/CHS_GxEScan/Data/" #Path from which to read cov/phen files
fam_path <- "/home/pmd-01/chemenya/CHS/txtDosage/" #Path from which to read whole chromosome data, fam, map
output_path <- "/home/pmd-01/chemenya/CHS/PM25/Imputed/" #Path to write results to    
G_path <- "/home/pmd-01/chemenya/CHS/Split_Imputed_Results/" #Path where rds files are
dimension_path <- "/home/pmd-01/chemenya/CHS/txtDosage/dimensions/"

#Read in Epi files
phen <- read.table(paste0(epi_path,"chs.pheno"),header=T) #Read in phenotype #3000 samples
fam <- read.table(paste0(fam_path,"chs3000.fam"),header=F) #Read in fam file #3000 Samples
names(fam) <- c("FID","IID","PID","MID","SEX","PHEN")
cov <- read.table(paste0(epi_path,"chs.cov"),header=T) #Read in Covariates #3000 Samples

#Match all data up to fam order
cov.matched <- cov[match(fam$IID,cov$IID),]
phen.matched <- phen[match(fam$IID,phen$IID),]

#Create useful Covariate file with only relevant variables
Cov <- data.frame( cbind(cov.matched$male,cov.matched$afr,cov.matched$natam,
                         cov.matched$asian,cov.matched$ses) )
names(Cov) <- c("male","afr","natam","asian","ses")
#fh = family history
#hw = hispanic whites

#Assign Vectors
#E <- cov$hw #Environmental Factor = hispanic whites
E <- ifelse(cov.matched$pm25 <= median(cov.matched$pm25),0,1) #Environmental Factor = dichotomized pm25
Y <- phen.matched$asthma #Phenotype = asthma  

#Categorize continuous variables into dichotomies
Cov$afr <- ifelse(Cov$afr<median(Cov$afr),0,1)
Cov$natam <- ifelse(Cov$natam<median(Cov$natam),0,1)
Cov$asian <- ifelse(Cov$asian<median(Cov$asian),0,1)

#Cov$ses1 <- If all are zero
Cov$ses2 <- ifelse(Cov$ses==2,1,0)
Cov$ses3 <- ifelse(Cov$ses==3,1,0)
Cov$ses4 <- ifelse(Cov$ses==4,1,0)
Cov$ses5 <- ifelse(Cov$ses==5,1,0)

cov1 <- Cov$male
cov2 <- Cov$afr
cov3 <- Cov$natam
cov4 <- Cov$asian
ses2 <- Cov$ses2
ses3 <- Cov$ses3
ses4 <- Cov$ses4
ses5 <- Cov$ses5

#Read in the number of SNPs to be read in the dosage file
dimension <- as.numeric(read.table(paste0(dimension_path,"chr",chr,".row"))[1])

#Read in how many files there are for this chromosome
num.files <- ceiling(dimension/5000)

#Which files to read
scripts <- ceiling(num.files/20)
read.these <- c((20*(script-1)):(20+20*(script-1)-1))
files <- sprintf(".%0004d", read.these)

#Read in file, and analyze and write results - loop
for(j in 1:length(files)){ #Do batches of about 20 files
  file <- paste0(G_path,"G.chr", chr,files[j])
  if(file.exists(file)){
    #Read in imputed RDS file and structure for analysis 
    gen <- readRDS(paste0(G_path,"G.chr", chr,files[j]))
    names(gen) <- c("Marker","A1","A2",as.character(fam$IID)) #Name the columns of gen file
    gen.mat <- as.matrix( gen[,4:ncol(gen)] ) #Take only the genotype information from subset.gen
    lab.mat <- as.matrix( gen[,1:3] ) #SNP Intro Matrix
    gen.final <- t(gen.mat) #IID by SNP matrix, ordered by IID, with IID removed
    
    #Analyze each SNP in the gen file
    system.time( SNP.res <- apply(gen.final,2,function(G){  #Run loop by columns of gen matrix
    #system.time( SNP.res <- apply(gen.final[,1:2],2,function(G){  #Run loop by columns of gen matrix
      G[G==-9] <- NA
      dat <- as.data.frame( cbind(Y,G,E,cov1,cov2,cov3,cov4,ses2,ses3,ses4,ses5) )
      Dat <- dat[complete.cases(dat),]
      #One-Step methods
      #########################
      BMA_DF2.MVW <- matrix(run.BMA.Mult(Dat,T,F,F),ncol=7) #BMA_DF2.MVW will output pval, PPM1(CC),PPM2(CO),PPSNP,PPSNPxE,GlimEst.SNP,GlimEst.SNPxE
      BMA2 <- matrix(run.BMA2(Dat,T,F,F),ncol=4)
      CC <- matrix(run.CC(Dat),ncol=4)
      CO <- matrix(run.CO(Dat),ncol=4)      
      GENOTYPE.Y <- matrix(run.GENOTYPE.Y(Dat),ncol=4)
      GENOTYPE.E <- matrix(run.GENOTYPE.E(Dat),ncol=4)
      DF2 <- matrix(run.DF2(Dat),ncol=1)
      CO_DF2 <- matrix(t(t(run.CO_DF2(GENOTYPE.Y[1],GENOTYPE.Y[2],CO[1],CO[2]))),ncol=1)
      line <- cbind(nrow(Dat),nrow(Dat[Dat$Y==1,]),BMA_DF2.MVW,BMA2,CC,CO,GENOTYPE.Y,GENOTYPE.E,DF2,CO_DF2)
    }) )
   #After all SNPs in file are analyzed, put it all together
    SNP.res <- as.data.frame( cbind(lab.mat,t(SNP.res)) )
    #SNP.res <- as.data.frame( cbind(lab.mat[1:2,],t(SNP.res)) )
    names(SNP.res) <- c("SNP","A1","A2","Samples","Cases",   #Intro
                        "BMA_P","BMA_PPM1","BMA_PPM2","BMA_GlimEst.SNP1","BMA_GlimEst.SNP2",  #BMA_DF2 (1)
                        "BMA_GlimEst.SNPxE1","BMA_GlimEst.SNPxE2",                            #BMA_DF2 (2)
                        "BMA2_Int.est", "BMA2_Int.sd", "BMA2_Z.score","BMA2_P",               #BMA2
                        "CC_Int.est", "CC_Int.sd", "CC_Z.score","CC_P",                       #CC
                        "CO_Int.est", "CO_Int.sd", "CO_Z.score","CO_P",                       #CO
                        "GENY_Int.est", "GENY_Int.sd", "GENY_Z.score","GENY_P",               #GENOTYPE.Y
                        "GENE_Int.est", "GENE_Int.sd", "GENE_Z.score","GENE_P",               #GENOTYPE.E
                        "DF2_P",                                                              #DF2
                        "CODF2_P")
    
    #Gather Results per method  
    BMA_DF2.MVW.res <- SNP.res[,c(1:5,6:12)]
    BMA2.res <- SNP.res[,c(1:5,13:16)]
    CC.res <- SNP.res[,c(1:5,17:20)]
    CO.res <- SNP.res[,c(1:5,21:24)]
    GENOTYPE.Y.res <- SNP.res[,c(1:5,25:28)]
    GENOTYPE.E.res <- SNP.res[,c(1:5,29:32)]
    DF2.res <- SNP.res[,c(1:5,33)]
    CO_DF2.res <- SNP.res[,c(1:5,34)]
    
    #Write results per method
    saveRDS(BMA_DF2.MVW.res,file=paste0(output_path,"BMA_DF2.chr",i,files[j],".rds"))
    saveRDS(BMA2.res,file=paste0(output_path,"BMA2.chr",i,files[j],".rds"))
    saveRDS(CC.res,file=paste0(output_path,"CC.chr",i,files[j],".rds"))
    saveRDS(CO.res,file=paste0(output_path,"CO.chr",i,files[j],".rds"))
    saveRDS(GENOTYPE.Y.res,file=paste0(output_path,"GENOTYPE.Y.chr",i,files[j],".rds"))
    saveRDS(GENOTYPE.E.res,file=paste0(output_path,"GENOTYPE.E.chr",i,files[j],".rds"))
    saveRDS(DF2.res,file=paste0(output_path,"DF2.chr",i,files[j],".rds"))
    saveRDS(CO_DF2.res,file=paste0(output_path,"CO_DF2.chr",i,files[j],".rds"))
   } else 
    print(paste0("No Such File:",file))
}






