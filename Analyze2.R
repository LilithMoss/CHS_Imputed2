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

CHR <- 22 #Set which chromosome will be processed
i=CHR
I_ge=0.5;I_g=0.5;I_gxe=0.5 #Set prior probabilities
#path <- "/auto/pmd-01/chemenya/CHS/CHS_GxEScan/Data/" #Path from which to read map/fam files
#path2 <- "/auto/pmd-01/chemenya/CHS/CHS_GxEScan/"     #Path from which to read subsetted gen files by chromosomes
#path3 <- "/auto/pmd-01/chemenya/CHS/PM25/"        #Path to write results to    

path <- "/auto/pmd-01/chemenya/CHS/CHS_GxEScan/Data/" #Path from which to read cov/phen files
path2 <- "/home/pmd-01/chemenya/CHS/txtDosage/" #Path from which to read whole chromosome data, fam, map
path3 <- "/home/pmd-01/chemenya/CHS/PM25/Imputed/" #Path to write results to    
path4 <- "/home/pmd-01/chemenya/CHS/Parallel/rds_files/"
dimension_path <- "/home/pmd-01/chemenya/CHS/txtDosage/dimensions/"

phen <- read.table(paste0(path,"chs.pheno"),header=T) #Read in phenotype #3000 samples
fam <- read.table(paste0(path2,"chs3000.fam"),header=F) #Read in fam file #3000 Samples
names(fam) <- c("FID","IID","PID","MID","SEX","PHEN")

#map <- read.table(paste0(path,"chs.map"),header=F) #Read in map file #630600 SNPs
#names(map) <- c("Chr","SNP","Dist","BP")
cov <- read.table(paste0(path,"chs.cov"),header=T) #Read in Covariates #3000 Samples

#Match all data up to fam order
cov.matched <- cov[match(fam$IID,cov$IID),]
phen.matched <- phen[match(fam$IID,phen$IID),]

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

# #Will read subsetted files
# chr <- do.call(rbind, lapply(1:26,function(i){
#   assign(paste("CHR",i,sep=""),range( which(map$Chr==i,)) )
# }) )

#Read in RDS files
#Read in the number of SNPs to be read in the dosage file
dimension <- as.numeric(read.table(paste0(dimension_path,"chr",CHR,".row"))[1])
#Set the number of Lines to be read from dosage file at a time
Lines <- 2000
#Calculate the number of rounds needed to read in the final file
loops <- ceiling(dimension/Lines)

#Read map
map <- read.table(paste0(path2,"chr",CHR,".map"))

##########################
# SERIAL
##########################
#for(grps in 1:loops){
for(grps in 1:4){
  system.time(gen <- readRDS(paste0(path4,"chr",CHR,"_",grps,".rds"))  )
  names(gen) <- c("Marker","A1","A2",as.character(fam$IID)) #Name the columns of gen file
  gen.mat <- as.matrix( gen[,4:ncol(gen)] ) #Take only the genotype information from subset.gen
  lab.mat <- as.matrix( gen[,1:3] ) #SNP Intro Matrix
  gen.final <- t(gen.mat) #IID by SNP matrix, ordered by IID, with IID removed
  strt <- 1+2000*(grps-1)
  stp <- 2000*(grps-1)+ncol(gen.final)
  gen.map <- map[strt:stp,]


    #system.time( SNP.res <- apply(gen.final,2,function(G){  #Run loop by columns of gen matrix
    system.time( SNP.res <- apply(gen.final[,1:2],2,function(G){  #Run loop by columns of gen matrix
      #G <- ifelse(C==0,0,ifelse(C==1|C==2,1,NA))
      Frequency <- mean(G[G>=0])/2 #C is the dosage and we use it to get the Frequency
      G[G==-9] <- NA
      dat <- as.data.frame( cbind(Y,G,E,cov1,cov2,cov3,cov4,ses2,ses3,ses4,ses5) )
      Dat <- dat[complete.cases(dat),]
      #One-Step methods
      #########################
      #BMA_DF2.MVW will output pval, PPM1(CC),PPM2(CO),PPSNP,PPSNPxE,GlimEst.SNP,GlimEst.SNPxE
      BMA_DF2.MVW <- matrix(run.BMA.Mult(Dat,T,F,F),ncol=7)
      BMA2 <- matrix(run.BMA2(Dat,T,F,F),ncol=4)
      CC <- matrix(run.CC(Dat),ncol=4)
      CO <- matrix(run.CO(Dat),ncol=4)      
      GENOTYPE.Y <- matrix(run.GENOTYPE.Y(Dat),ncol=4)
      GENOTYPE.E <- matrix(run.GENOTYPE.E(Dat),ncol=4)
      DF2 <- matrix(run.DF2(Dat),ncol=1)
      CO_DF2 <- matrix(t(t(run.CO_DF2(GENOTYPE.Y[1],GENOTYPE.Y[2],CO[1],CO[2]))),ncol=1)
      line <- cbind(nrow(Dat),nrow(Dat[Dat$Y==1,]),Frequency,BMA_DF2.MVW,BMA2,CC,CO,GENOTYPE.Y,GENOTYPE.E,DF2,CO_DF2)
    }) ) 
    SNP.res2 <- cbind(lab.mat, t(SNP.res) )
    SNP.res2 <- cbind(gen.map,lab.mat,SNP.res)
    SNP.res2 <- cbind(gen.map[1:2,],lab.mat[1:2,],t(SNP.res))
    
    names(SNP.res2) <- c("CHR","SNP","DIST","BP","GenSNP","A1","A2","Samples","Cases","MAF",   #Intro
                         "BMA_P","BMA_PPM1","BMA_PPM2","BMA_GlimEst.SNP1","BMA_GlimEst.SNP2",  #BMA_DF2 (1)
                         "BMA_GlimEst.SNPxE1","BMA_GlimEst.SNPxE2",                            #BMA_DF2 (2)
                         "BMA2_Int.est", "BMA2_Int.sd", "BMA2_Z.score","BMA2_P",               #BMA2
                         "CC_Int.est", "CC_Int.sd", "CC_Z.score","CC_P",                       #CC
                         "CO_Int.est", "CO_Int.sd", "CO_Z.score","CO_P",                       #CO
                         "GENY_Int.est", "GENY_Int.sd", "GENY_Z.score","GENY_P",               #GENOTYPE.Y
                         "GENE_Int.est", "GENE_Int.sd", "GENE_Z.score","GENE_P",               #GENOTYPE.E
                         "DF2_P",                                                              #DF2
                         "CODF2_P")                                                            #CO_DF2

} #End of all grps for this chromosome






########################
# APPLY - Output Matrix
########################
# m.labels <- map[map$Chr==CHR,] #Subset the SNP names from the relevant chromosome
# gen.mat <- as.matrix( subset.gen[,4:ncol(subset.gen),with=FALSE] ) #Take only the genotype information from subset.gen
# lab.mat <- cbind( m.labels,as.matrix( subset.gen[,1:3,with=FALSE] ) ) #SNP Intro Matrix
# nam <- as.numeric(colnames(gen.mat)) #Take IIDs in order of the gen file
# tgen.mat <- cbind(nam,t(gen.mat))    #Transpose the gen matrix and add IIDs as row names
# tgen.mat.sort <- tgen.mat[order(tgen.mat[,1]),] #Order by IID
# gen.final <- tgen.mat.sort[,-1] #IID by SNP matrix, ordered by IID, with IID removed
# 
# #Garbage Collection to free up memory
# rm(gen)
# rm(map)
# rm(gen.mat)
# rm(tgen.mat)
# rm(tgen.mat.sort)

#APPLY FUNCTION
#gen.test <- gen.final[,1000:1100]
#lab.mat <- lab.mat[1000:1100,]

#system.time( SNP.res <- apply(gen.test,2,function(C){
#C=gen.final[,150] #Remove when not testing
system.time( SNP.res <- apply(gen.final,2,function(C){  #Run loop by columns of gen matrix
  G <- ifelse(C==0,0,ifelse(C==1|C==2,1,NA))
  dat <- as.data.frame( cbind(Y,G,E,cov1,cov2,cov3,cov4,ses2,ses3,ses4,ses5) )
  Dat <- dat[complete.cases(dat),]
  Frequency <- mean(C[C>=0])/2 #C is the dosage and we use it to get the Frequency
  #One-Step methods
  #########################
  #BMA_DF2.MVW will output pval, PPM1(CC),PPM2(CO),PPSNP,PPSNPxE,GlimEst.SNP,GlimEst.SNPxE
  BMA_DF2.MVW <- matrix(run.BMA.Mult(Dat,T,F,F),ncol=7)
  BMA2 <- matrix(run.BMA2(Dat,T,F,F),ncol=4)
  CC <- matrix(run.CC(Dat),ncol=4)
  CO <- matrix(run.CO(Dat),ncol=4)      
  GENOTYPE.Y <- matrix(run.GENOTYPE.Y(Dat),ncol=4)
  GENOTYPE.E <- matrix(run.GENOTYPE.E(Dat),ncol=4)
  DF2 <- matrix(run.DF2(Dat),ncol=1)
  CO_DF2 <- matrix(t(t(run.CO_DF2(GENOTYPE.Y[1],GENOTYPE.Y[2],CO[1],CO[2]))),ncol=1)
  line <- cbind(nrow(Dat),nrow(Dat[Dat$Y==1,]),Frequency,BMA_DF2.MVW,BMA2,CC,CO,GENOTYPE.Y,GENOTYPE.E,DF2,CO_DF2)
}) )
#lab.mat2 <- lab.mat[1000:1100,]
SNP.res2 <- cbind(lab.mat, t(SNP.res) )
names(SNP.res2) <- c("CHR","SNP","DIST","BP","GenSNP","A1","A2","Samples","Cases","MAF",   #Intro
                     "BMA_P","BMA_PPM1","BMA_PPM2","BMA_GlimEst.SNP1","BMA_GlimEst.SNP2",  #BMA_DF2 (1)
                     "BMA_GlimEst.SNPxE1","BMA_GlimEst.SNPxE2",                            #BMA_DF2 (2)
                     "BMA2_Int.est", "BMA2_Int.sd", "BMA2_Z.score","BMA2_P",               #BMA2
                     "CC_Int.est", "CC_Int.sd", "CC_Z.score","CC_P",                       #CC
                     "CO_Int.est", "CO_Int.sd", "CO_Z.score","CO_P",                       #CO
                     "GENY_Int.est", "GENY_Int.sd", "GENY_Z.score","GENY_P",               #GENOTYPE.Y
                     "GENE_Int.est", "GENE_Int.sd", "GENE_Z.score","GENE_P",               #GENOTYPE.E
                     "DF2_P",                                                              #DF2
                     "CODF2_P")                                                            #CO_DF2

BMA_DF2.MVW.res <- SNP.res2[,c(1:10,11:17)]
BMA2.res <- SNP.res2[,c(1:10,18:21)]
CC.res <- SNP.res2[,c(1:10,22:25)]
CO.res <- SNP.res2[,c(1:10,26:29)]
GENOTYPE.Y.res <- SNP.res2[,c(1:10,30:33)]
GENOTYPE.E.res <- SNP.res2[,c(1:10,34:37)]
DF2.res <- SNP.res2[,c(1:10,38)]
CO_DF2.res <- SNP.res2[,c(1:10,39)]

write.table(BMA_DF2.MVW.res,file=paste0(path3,"BMA_DF2.res.CHR",i,".txt"),row.names=F)
write.table(BMA2.res,file=paste0(path3,"BMA2.res.CHR",i,".txt"),row.names=F)
write.table(CC.res,file=paste0(path3,"CC.res.CHR",i,".txt"),row.names=F)
write.table(CO.res,file=paste0(path3,"CO.res.CHR",i,".txt"),row.names=F)
write.table(GENOTYPE.Y.res,file=paste0(path3,"GENOTYPE.Y.res.CHR",i,".txt"),row.names=F)
write.table(GENOTYPE.E.res,file=paste0(path3,"GENOTYPE.E.res.CHR",i,".txt"),row.names=F)
write.table(DF2.res,file=paste0(path3,"DF2.res.CHR",i,".txt"),row.names=F)
write.table(CO_DF2.res,file=paste0(path3,"CO_DF2.res.CHR",i,".txt"),row.names=F)


