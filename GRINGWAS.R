
#download phenotype data
setwd("~/Downloads")
library(readxl)
caro<- read_excel("carotenoid_quant_final_mean.xlsx")

#### Load the genotype file
setwd("~/Desktop/GWAS")
# read geno
#always read as delim and use strings as factors
read.delim("geno_R_no_duplicates_MAF.hmp.txt", stringsAsFactors=F)-> geno
geno<- rbind( sub('\\..+', '', names(geno)), geno)


#if genotype file is missing PIs it will give an error that 'x' must be numeric--> 
# filter out phenotype data with missing genotype data
genoonly<-data.frame(colnames(geno[,12:407]))
colnames(genoonly)[1]= "PI"
colnames(caro)[1]= "PI"

#pheno file only has phenotypes that have genotypes #312 
pheno <- merge(caro, genoonly, by = "PI")


setwd("~/Desktop/GRIN/ColorLab") # or where ever you what the files to go

#sometimes you have to coerce packages to work, i usually dont need to use this
#if (!require("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
#BiocManager::install("multtest")

#load dependent packages
library(multtest)
library(gplots)
library(LDheatmap)
library(genetics)
library(ape)
library(EMMREML)
library(compiler) #this library is already installed in R
library("scatterplot3d")
library("ggplot2")
library(devtools)

## GETTING GAPIT SOURCE CODE
source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("https://zzlab.net/GAPIT/gapit_functions.txt")


# Prepping our variables for GWAS, needs to be numeric
labpheno=pheno[,c(1,3:5,18)]
labpheno$Color=as.factor(labpheno$Color)
levels(labpheno$Color)
summary(labpheno$Color)
#1brown #2red #3white #4yellow
labpheno$Color=as.numeric(labpheno$Color)


#run this for only white and yellow grain
labpheno <- subset(labpheno,         
                Color == "4"|
                  Color == "3")

#############Quantitative GWAS for Carotenoid content######

#phenotype or Y
caropheno<- pheno[,c(1,6,7,9)]
caropheno$Total= pheno$Lutein + pheno$Zeaxanthin + pheno$BCarotene 

setwd("~/Desktop/GRIN/Carotenoid")
myGAPIT_caro <- GAPIT(
  Y=caropheno,
  G=geno,
  PCA.total = 7,
  SNP.MAF= 0.05, 
  model="BLINK")

############# Quantitative Color GWAS + Color

setwd("~/Desktop/GRIN/ColorLAB")
myGAPIT_caro <- GAPIT(
  Y=labpheno,
  G=geno,
  PCA.total = 7,
  SNP.MAF= 0.05, 
  model="BLINK")

###########################################################

#PCAs from GAPIT data
test<- myGAPIT_LAB$PCA
colnames(test)[1]= "PI"
pheno <- merge(together, test, by = "PI")
  ggplot(pheno,
         aes(
           PC1, 
           PC2,
           color= "Color"
         ))+
  geom_point(alpha=0.8)

  test2<- myGAPIT_LAB$PCA
  colnames(test2)[1]= "PI"
  pheno2 <- merge(labpheno, test2, by = "PI")
  pheno3 <- merge(pheno3, caropheno, by = "PI")
  pheno$Color <- sub(pattern = "1", replacement = "brown",pheno$Color)
  pheno$Color <- sub(pattern = "2", replacement = "red",pheno$Color)
  pheno$Color <- sub(pattern = "3", replacement = "white",pheno$Color)
  pheno$Color <- sub(pattern = "4", replacement = "yellow",pheno$Color)
  
  peri_color_lookup <- c(
    red = 'red',
    white = 'gray',
    brown = 'brown',
    yellow = 'orange')
  
  
  peri_clr_vec <- peri_color_lookup[match(pheno$Color, names(peri_color_lookup))]
  
  
  pheno=pheno[,c(1,5,10:12)]
  
  pheno=na.omit(pheno)
  
  
  
  pheno2$Color <-as.factor(pheno2$Color)
  ggplot(pheno,
         aes(PC1,PC3,
           color=Color,))+
    geom_point(alpha=0.8, size=3) + scale_colour_manual(values = peri_clr_vec, breaks=c("brown","red","white", "yellow")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                            panel.background = element_blank(),axis.line=element_line(colour="black")) + theme(axis.text.x = element_text(size = 20, colour="black")) + theme(axis.text.y = element_text(size = 20, colour="black")) + theme(text=element_text(size=20,colour="black")) 
  
 #################################### MANHATTAN PLOTS 
# load needed packages
install.packages("qqman")
library(qqman)

install.packages("glmnet", dep=TRUE) 
library(glmnet)

##########    ZEAXANTHIN.  ##################
zeaBLINK <-read.csv('~/Desktop/GRIN/Carotenoid/GAPIT.Association.GWAS_Results.BLINK.Zeaxanthin.csv',head=TRUE,sep=',',stringsAsFactors=F)
qq(zeaBLINK$P.value, main = "Q-Q plot of GWAS p-values")
zeaBLINK$FDR=p.adjust(zeaBLINK$P.value, method="fdr")

zeaBLINK_snpsOfInterest <- zeaBLINK[zeaBLINK$FDR <  0.1, "SNP"]
print(zeaBLINK_snpsOfInterest)

zeaBLINK <- zeaBLINK[,c(1:4,9)]
colnames(zeaBLINK)[2] <- "CHR"
colnames(zeaBLINK)[3] <- "BP"
colnames(zeaBLINK)[4] <- "P"
colnames(zeaBLINK)[5] <- "FDR"


manhattan(zeaBLINK, highlight= zeaBLINK_snpsOfInterest, main= "Zeaxanthin", suggestiveline=F, genomewideline=F)
abline(h=(FDR5),col="blue", lwd=1)
abline(h=(FDR10), col= "blue", lwd=1)


print(zeaBLINK_snpsOfInterest)


FDR5=-log10(min(subset(zeaBLINK,FDR>0.049)$P))
FDR10=-log10(min(subset(zeaBLINK,FDR>0.1)$P))

###### CHR 1 LUT1
manhattan(subset(zeaBLINK, CHR == 1), highlight= zeaBLINK_snpsOfInterest, main= "Zeaxanthin BLINK", suggestiveline = FALSE, genomewideline = FALSE) 
abline(v=(64901748), col= "hotpink", lwd=1)
abline(h=(FDR5), col= "blue", lwd=1)
abline(h=(FDR10), col= "blue", lwd=1)

####### CHR 3 DXR1
manhattan(subset(zeaBLINK, CHR == 3), highlight= zeaBLINK_snpsOfInterest, main= "Zeaxanthin BLINK", suggestiveline = FALSE, genomewideline = FALSE) 
abline(v=(9218437), col= "hotpink", lwd=1)
abline(v=(4446539), col= "hotpink", lwd=1)
abline(h=(FDR5), col= "blue", lwd=1)
abline(h=(FDR10), col= "blue", lwd=1)

###### CH6 ZEP
manhattan(subset(zeaBLINK, CHR == 6), highlight= zeaBLINK_snpsOfInterest, main= "Zeaxanthin BLINK", suggestiveline = FALSE, genomewideline = FALSE) 
abline(v=(46715530), col= "hotpink", lwd=1)
abline(h=(FDR5), col= "blue", lwd=1)
abline(h=(FDR10), col= "blue", lwd=1)

###### CH7 CCD8
manhattan(subset(zeaBLINK, CHR == 7), highlight= zeaBLINK_snpsOfInterest, main= "Zeaxanthin BLINK", suggestiveline = FALSE, genomewideline = FALSE) 
abline(v=(60503028), col= "hotpink", lwd=1)
abline(h=(FDR5), col= "blue", lwd=1)
abline(h=(FDR10), col= "blue", lwd=1)

################################### LUTEIN ######################

lutBLINK <-read.csv('~/Desktop/GRIN/Carotenoid/GAPIT.Association.GWAS_Results.BLINK.Lutein.csv',head=TRUE,sep=',',stringsAsFactors=F)
qq(lutBLINK$P.value, main = "Q-Q plot of GWAS p-values")
lutBLINK$FDR=p.adjust(lutBLINK$P.value, method="fdr")

lutBLINK_snpsOfInterest <- lutBLINK[lutBLINK$FDR <  0.1, "SNP"]
print(lutBLINK_snpsOfInterest)

lutBLINK <- lutBLINK[,c(1:4,9)]

colnames(lutBLINK)[2] <- "CHR"
colnames(lutBLINK)[3] <- "BP"
colnames(lutBLINK)[4] <- "P"
colnames(lutBLINK)[5] <- "FDR"

FDR5=-log10(min(subset(lutBLINK,FDR>0.049)$P))
FDR10=-log10(min(subset(lutBLINK,FDR>0.1)$P))


manhattan(lutBLINK, ylim = c(0, 15), genomewideline=FALSE,suggestiveline=FALSE,highlight=lutBLINK_snpsOfInterest, main="Lutein")
abline(h=(FDR5), col= "blue", lwd=1)
abline(h=(FDR10), col= "blue", lwd=1)





manhattan(subset(lutBLINK, CHR == 10), highlight= lutBLINK_snpsOfInterest, main= "Lutein Carotenoid BLINK") 
abline(v=(60960084), col= "hotpink", lwd=1)

####### CHR 2 DXS
manhattan(subset(lutBLINK, CHR == 2), highlight= lutBLINK_snpsOfInterest, main= "Lutein BLINK", suggestiveline = FALSE, genomewideline = FALSE) 
abline(v=(6275655), col= "hotpink", lwd=1)
abline(h=(FDR5), col= "blue", lwd=1)
abline(h=(FDR10), col= "blue", lwd=1)

####### CHR 3 LYCE and CYP86A1
manhattan(subset(lutBLINK, CHR == 3), highlight= lutBLINK_snpsOfInterest, main= "Lutein BLINK", suggestiveline = FALSE, genomewideline = FALSE) 
abline(v=(67816770), col= "hotpink", lwd=1)
abline(v=(52269412), col= "hotpink", lwd=1)
abline(h=(FDR5), col= "blue", lwd=1)
abline(h=(FDR10), col= "blue", lwd=1)

###### CH 4 NCED
manhattan(subset(lutBLINK, CHR == 4), highlight= lutBLINK_snpsOfInterest, main= "Lutein BLINK", suggestiveline = FALSE, genomewideline = FALSE) 
abline(v=(61268234), col= "hotpink", lwd=1)
abline(h=(FDR5), col= "blue", lwd=1)
abline(h=(FDR10), col= "blue", lwd=1)

###### CH 10 PSY
manhattan(subset(lutBLINK, CHR == 10), highlight= lutBLINK_snpsOfInterest, main= "Lutein BLINK", suggestiveline = FALSE, genomewideline = FALSE) 
abline(v=(60960084), col= "hotpink", lwd=1)
abline(h=(FDR5), col= "blue", lwd=1)
abline(h=(FDR10), col= "blue", lwd=1)



#################################### BETACAROTENE #############################
bcBLINK <-read.csv('~/Desktop/GRIN/Carotenoid/GAPIT.Association.GWAS_Results.BLINK.BCarotene.csv',head=TRUE,sep=',',stringsAsFactors=F)
qq(bcBLINK$P.value, main = "Q-Q plot of GWAS p-values")
bcBLINK$FDR=p.adjust(bcBLINK$P.value, method="fdr")


colnames(bcBLINK)[2] <- "CHR"
colnames(bcBLINK)[3] <- "BP"
colnames(bcBLINK)[4] <- "P"
colnames(bcBLINK)[9] <- "FDR"

bcBLINK <- bcBLINK[,c(1:4,9)]
bcBLINK_snpsOfInterest <- bcBLINK[bcBLINK$FDR <  0.1, "SNP"]
print(bcBLINK_snpsOfInterest)

FDR5=-log10(min(subset(bcBLINK,FDR>0.049)$P))
FDR10=-log10(min(subset(bcBLINK,FDR>0.1)$P))

manhattan(bcBLINK,highlight= bcBLINK_snpsOfInterest, main= "Beta-carotene", suggestiveline = FALSE, genomewideline = FALSE) 
abline(h=(FDR5), col= "blue", lwd=1)
abline(h=(FDR10), col= "blue", lwd=1)

################# CHR 2 ZN finger protein
manhattan(subset(bcBLINK, CHR == 2), highlight= bcBLINK_snpsOfInterest, main= "Beta-carotene BLINK", suggestiveline = FALSE, genomewideline = FALSE) 
abline(v=(60899412), col= "hotpink", lwd=1)
abline(h=(FDR5), col= "blue", lwd=1)
abline(h=(FDR10), col= "blue", lwd=1)


######## CH 10 DXS
manhattan(subset(bcBLINK, CHR == 10), highlight= bcBLINK_snpsOfInterest, main= "Beta-carotene BLINK", suggestiveline = FALSE, genomewideline = FALSE, ylim=c(0,7)) 
abline(v=(2611880), col= "hotpink", lwd=1)
abline(h=(FDR5), col= "blue", lwd=1)
abline(h=(FDR10), col= "blue", lwd=1)

######################################## TOTAL CAROTENOIDS ###########
totBLINK <-read.csv('~/Desktop/GRIN/Carotenoid/GAPIT.Association.GWAS_Results.BLINK.Total.csv',head=TRUE,sep=',',stringsAsFactors=F)
qq(totBLINK$P.value, main = "Q-Q plot of GWAS p-values")
totBLINK$FDR=p.adjust(totBLINK$P.value, method="fdr")


colnames(totBLINK)[2] <- "CHR"
colnames(totBLINK)[3] <- "BP"
colnames(totBLINK)[4] <- "P"
colnames(totBLINK)[9] <- "FDR"

totBLINK <- totBLINK[,c(1:4,9)]
totBLINK_snpsOfInterest <- totBLINK[totBLINK$FDR <  0.1, "SNP"]
print(totBLINK_snpsOfInterest)
FDR5=-log10(min(subset(totBLINK,FDR>0.049)$P))
FDR10=-log10(min(subset(totBLINK,FDR>0.1)$P))


manhattan(totBLINK, highlight= totBLINK_snpsOfInterest, main= "Total Carotenoids", suggestiveline = FALSE, genomewideline = FALSE)
abline(h=(FDR5), col= "blue", lwd=1)
abline(h=(FDR10), col= "blue", lwd=1)

######## CHR 2 CCD
manhattan(subset(totBLINK, CHR == 2), highlight= totBLINK_snpsOfInterest, main= "Total Carotenoid BLINK", suggestiveline = FALSE, genomewideline = FALSE) 
abline(v=(53189951), col= "hotpink", lwd=1)
abline(h=(FDR5), col= "blue", lwd=1)
abline(h=(FDR10), col= "blue", lwd=1)

######### CHR 3 DXR
manhattan(subset(totBLINK, CHR == 3), highlight= totBLINK_snpsOfInterest, main= "Total Carotenoid BLINK", suggestiveline = FALSE, genomewideline = FALSE) 
abline(v=(9218437), col= "hotpink", lwd=1)
abline(h=(FDR5), col= "blue", lwd=1)
abline(h=(FDR10), col= "blue", lwd=1)

######## CHR 6 ZEP
manhattan(subset(totBLINK, CHR == 6), highlight= totBLINK_snpsOfInterest, main= "Total Carotenoid BLINK", suggestiveline = FALSE, genomewideline = FALSE) 
abline(v=(46715530), col= "hotpink", lwd=1)
abline(h=(FDR5), col= "blue", lwd=1)
abline(h=(FDR10), col= "blue", lwd=1)


####### CHR 8 Transparent Testa
manhattan(subset(totBLINK, CHR == 8), highlight= totBLINK_snpsOfInterest, main= "Total Carotenoid BLINK", suggestiveline = FALSE, genomewideline = FALSE) 
abline(v=(60536621), col= "hotpink", lwd=1)
abline(h=(FDR5), col= "blue", lwd=1)
abline(h=(FDR10), col= "blue", lwd=1)


################################################# COLOR ###########

colorBLINK <-read.csv('~/Desktop/GRIN/ColorLab/GAPIT.Association.GWAS_Results.BLINK.Color.csv',head=TRUE,sep=',',stringsAsFactors=F)
qq(colorBLINK$P.value, main = "Q-Q plot of GWAS p-values")
colorBLINK$FDR=p.adjust(colorBLINK$P.value, method="fdr")
colnames(colorBLINK)[2] <- "CHR"
colnames(colorBLINK)[3] <- "BP"
colnames(colorBLINK)[4] <- "P"
colnames(colorBLINK)[9] <- "FDR"

colorBLINK <- colorBLINK[,c(1:4,9)]
colorBLINK_snpsOfInterest <- colorBLINK[colorBLINK$FDR <  0.1, "SNP"]
print(colorBLINK_snpsOfInterest)
FDR5=-log10(min(subset(colorBLINK,FDR>0.049)$P))
FDR10=-log10(min(subset(colorBLINK,FDR>0.1)$P))

manhattan(colorBLINK, highlight= colorBLINK_snpsOfInterest, main= "Color", suggestiveline = FALSE, genomewideline = FALSE) 
abline(h=(FDR5), col= "blue", lwd=1)
abline(h=(FDR10), col= "blue", lwd=1)

### CHR 4 LYCb, lut2/lut5, glutathione S-transferase
manhattan(subset(colorBLINK, CHR == 4), highlight= colorBLINK_snpsOfInterest, main= "COLOR BLINK", suggestiveline = FALSE, genomewideline = FALSE) 
abline(v=(4637479), col= "hotpink", lwd=1)
abline(v=(6008130), col= "hotpink", lwd=1)
abline(v=(57980100), col= "hotpink", lwd=1)
abline(h=(FDR5), col= "blue", lwd=1)
abline(h=(FDR10), col= "blue", lwd=1)


###### CHR 6 PAP/fibrillin
manhattan(subset(colorBLINK, CHR == 6), highlight= colorBLINK_snpsOfInterest, main= "COLOR BLINK", suggestiveline = FALSE, genomewideline = FALSE) 
abline(v=(59325649), col= "hotpink", lwd=1)
abline(h=(FDR5), col= "blue", lwd=1)
abline(h=(FDR10), col= "blue", lwd=1)


######## CHR 9 ABA receptor

manhattan(subset(colorBLINK, CHR == 9), highlight= colorBLINK_snpsOfInterest, main= "COLOR BLINK", suggestiveline = FALSE, genomewideline = FALSE) 
abline(v=(52641821), col= "hotpink", lwd=1)
abline(v=(6867136), col= "hotpink", lwd=1)
abline(h=(FDR5), col= "blue", lwd=1)
abline(h=(FDR10), col= "blue", lwd=1)

######################################## MEAN L ########

lBLINK <-read.csv('~/Desktop/GRIN/ColorLab/GAPIT.Association.GWAS_Results.BLINK.MeanL.csv',head=TRUE,sep=',',stringsAsFactors=F)
qq(lBLINK$P.value, main = "Q-Q plot of GWAS p-values")
lBLINK$FDR=p.adjust(lBLINK$P.value, method="fdr")
colnames(lBLINK)[2] <- "CHR"
colnames(lBLINK)[3] <- "BP"
colnames(lBLINK)[4] <- "P"
colnames(lBLINK)[9] <- "FDR"

lBLINK <- lBLINK[,c(1:4,9)]
lBLINK_snpsOfInterest <- lBLINK[lBLINK$FDR <  0.1, "SNP"]
print(lBLINK_snpsOfInterest)

FDR5=-log10(min(subset(lBLINK,FDR>0.049)$P))
FDR10=-log10(min(subset(lBLINK,FDR>0.1)$P))

manhattan(lBLINK, highlight= lBLINK_snpsOfInterest, main= " L ", suggestiveline = FALSE, genomewideline = FALSE) 
abline(h=(FDR5), col= "blue", lwd=1)
abline(h=(FDR10), col= "blue", lwd=1)

####### CHR 2 PSY
manhattan(subset(lBLINK, CHR == 2), highlight= lBLINK_snpsOfInterest, main= "L BLINK", suggestiveline = FALSE, genomewideline = FALSE) 
abline(v=(67031564), col= "hotpink", lwd=1)
abline(h=(FDR5), col= "blue", lwd=1)
abline(h=(FDR10), col= "blue", lwd=1)

####### CHR 5 Flavonol 3-O-glucosyltransferase
manhattan(subset(lBLINK, CHR == 5), highlight= lBLINK_snpsOfInterest, main= "L BLINK", suggestiveline = FALSE, genomewideline = FALSE) 
abline(v=(3055093), col= "hotpink", lwd=1)
abline(h=(FDR5), col= "blue", lwd=1)
abline(h=(FDR10), col= "blue", lwd=1)

###### CHR 6 ZDS

manhattan(subset(lBLINK, CHR == 6), highlight= lBLINK_snpsOfInterest, main= "L BLINK", suggestiveline = FALSE, genomewideline = FALSE) 
abline(v=(53228503), col= "hotpink", lwd=1)
abline(h=(FDR5), col= "blue", lwd=1)
abline(h=(FDR10), col= "blue", lwd=1)

##### CHR 9 ABA receptor
manhattan(subset(lBLINK, CHR == 9), highlight= lBLINK_snpsOfInterest, main= "L BLINK", suggestiveline = FALSE, genomewideline = FALSE) 
abline(v=(52641821), col= "hotpink", lwd=1)
abline(h=(FDR5), col= "blue", lwd=1)
abline(h=(FDR10), col= "blue", lwd=1)

########################################## MEAN A #####

aBLINK <-read.csv('~/Desktop/GRIN/ColorLab/GAPIT.Association.GWAS_Results.BLINK.MeanA.csv',head=TRUE,sep=',',stringsAsFactors=F)
qq(aBLINK$P.value, main = "Q-Q plot of GWAS p-values")
aBLINK$FDR=p.adjust(aBLINK$P.value, method="fdr")
colnames(aBLINK)[2] <- "CHR"
colnames(aBLINK)[3] <- "BP"
colnames(aBLINK)[4] <- "P"
colnames(aBLINK)[9] <- "FDR"

aBLINK <- aBLINK[,c(1:4,9)]
aBLINK_snpsOfInterest <- aBLINK[aBLINK$FDR <  0.1, "SNP"]
print(aBLINK_snpsOfInterest)
FDR5=-log10(min(subset(aBLINK,FDR>0.049)$P))
FDR10=-log10(min(subset(aBLINK,FDR>0.1)$P))


manhattan(aBLINK, highlight= aBLINK_snpsOfInterest, main= "A", ) 
abline(h=(FDR5), col= "blue", lwd=1)
abline(h=(FDR10), col= "blue", lwd=1)




#### CHR 2 mo cofactor sulfurase

manhattan(subset(aBLINK, CHR == 2), highlight= aBLINK_snpsOfInterest, main= "a BLINK", suggestiveline = FALSE, genomewideline = FALSE) 
abline(v=(67340632), col= "hotpink", lwd=1)
abline(h=(FDR5), col= "blue", lwd=1)
abline(h=(FDR10), col= "blue", lwd=1)

######################################## MEAN B ######

bBLINK <-read.csv('~/Desktop/GRIN/ColorLab/GAPIT.Association.GWAS_Results.BLINK.MeanB.csv',head=TRUE,sep=',',stringsAsFactors=F)
qq(bBLINK$P.value, main = "Q-Q plot of GWAS p-values")
bBLINK$FDR=p.adjust(bBLINK$P.value, method="fdr")
colnames(bBLINK)[2] <- "CHR"
colnames(bBLINK)[3] <- "BP"
colnames(bBLINK)[4] <- "P"
colnames(bBLINK)[9] <- "FDR"

bBLINK <- bBLINK[,c(1:4,9)]
bBLINK_snpsOfInterest <- bBLINK[bBLINK$FDR <  0.1, "SNP"]
print(bBLINK_snpsOfInterest)
FDR5=-log10(min(subset(bBLINK,FDR>0.049)$P))
FDR10=-log10(min(subset(bBLINK,FDR>0.1)$P))


manhattan(bBLINK, highlight= bBLINK_snpsOfInterest, main= "B", suggestiveline = FALSE, genomewideline = FALSE) 
abline(h=(FDR5), col= "blue", lwd=1)
abline(h=(FDR10), col= "blue", lwd=1)

############################Just Yellow and White############

ywcolorBLINK <-read.csv('~/Desktop/GRIN/YWonlyColorLab/GAPIT.Association.GWAS_Results.BLINK.Color.csv',head=TRUE,sep=',',stringsAsFactors=F)
qq(ywcolorBLINK$P.value, main = "Q-Q plot of GWAS p-values")
ywcolorBLINK$FDR=p.adjust(ywcolorBLINK$P.value, method="fdr")
colnames(ywcolorBLINK)[2] <- "CHR"
colnames(ywcolorBLINK)[3] <- "BP"
colnames(ywcolorBLINK)[4] <- "P"
colnames(ywcolorBLINK)[9] <- "FDR"

ywcolorBLINK <- ywcolorBLINK[,c(1:4,9)]
ywcolorBLINK_snpsOfInterest <- ywcolorBLINK[ywcolorBLINK$FDR <  0.1, "SNP"]
print(ywcolorBLINK_snpsOfInterest)
FDR5=-log10(min(subset(ywcolorBLINK,FDR>0.049)$P))
FDR10=-log10(min(subset(ywcolorBLINK,FDR>0.1)$P))

manhattan(ywcolorBLINK, highlight= ywcolorBLINK_snpsOfInterest, main= "Color: Yellow and White", suggestiveline = FALSE, genomewideline = FALSE) 
abline(h=(FDR5), col= "blue", lwd=1)
abline(h=(FDR10), col= "blue", lwd=1)


####### CHR 4 CYP97A3/LUT5

manhattan(subset(ywcolorBLINK, CHR == 4), highlight= ywcolorBLINK_snpsOfInterest, main= "YW Color BLINK", suggestiveline = FALSE, genomewideline = FALSE) 
abline(v=(67575995), col= "hotpink", lwd=1)
abline(h=(FDR5), col= "blue", lwd=1)
abline(h=(FDR10), col= "blue", lwd=1)

###### CHR 10 DXS
manhattan(subset(ywcolorBLINK, CHR == 10), highlight= ywcolorBLINK_snpsOfInterest, main= "YW Color BLINK", suggestiveline = FALSE, genomewideline = FALSE) 
abline(v=(2611880), col= "hotpink", lwd=1)
abline(h=(FDR5), col= "blue", lwd=1)
abline(h=(FDR10), col= "blue", lwd=1)


############################# YW MEAN L ###

ywlBLINK <-read.csv('~/Desktop/GRIN/YWonlyColorLab/GAPIT.Association.GWAS_Results.BLINK.MeanL.csv',head=TRUE,sep=',',stringsAsFactors=F)
qq(ywlBLINK$P.value, main = "Q-Q plot of GWAS p-values")
ywlBLINK$FDR=p.adjust(ywlBLINK$P.value, method="fdr")
colnames(ywlBLINK)[2] <- "CHR"
colnames(ywlBLINK)[3] <- "BP"
colnames(ywlBLINK)[4] <- "P"
colnames(ywlBLINK)[9] <- "FDR"

ywlBLINK <- ywlBLINK[,c(1:4,9)]
ywlBLINK_snpsOfInterest <- ywlBLINK[lBLINK$FDR <  0.1, "SNP"]
print(ywlBLINK_snpsOfInterest)

FDR5=-log10(min(subset(ywlBLINK,FDR>0.049)$P))
FDR10=-log10(min(subset(ywlBLINK,FDR>0.1)$P))

manhattan(ywlBLINK, highlight= ywlBLINK_snpsOfInterest, main= "L", suggestiveline = FALSE, genomewideline = FALSE)
abline(h=(FDR5), col= "blue", lwd=1)
abline(h=(FDR10), col= "blue", lwd=1)

manhattan(subset(ywlBLINK, CHR == 2), highlight= ywlBLINK_snpsOfInterest, main= "YW L BLINK") 
abline(v=(67031564), col= "hotpink", lwd=1)

######################## YW MEAN A ########

ywaBLINK <-read.csv('~/Desktop/GRIN/YWonlyColorLab/GAPIT.Association.GWAS_Results.BLINK.MeanA.csv',head=TRUE,sep=',',stringsAsFactors=F)
qq(ywaBLINK$P.value, main = "Q-Q plot of GWAS p-values")
ywaBLINK$FDR=p.adjust(ywaBLINK$P.value, method="fdr")
colnames(ywaBLINK)[2] <- "CHR"
colnames(ywaBLINK)[3] <- "BP"
colnames(ywaBLINK)[4] <- "P"
colnames(ywaBLINK)[9] <- "FDR"

ywaBLINK <- ywaBLINK[,c(1:4,9)]
ywaBLINK_snpsOfInterest <- ywaBLINK[ywaBLINK$FDR <  0.1, "SNP"]
print(ywaBLINK_snpsOfInterest)
FDR5=-log10(min(subset(ywaBLINK,FDR>0.049)$P))
FDR10=-log10(min(subset(ywaBLINK,FDR>0.1)$P))

manhattan(ywaBLINK, highlight= ywaBLINK_snpsOfInterest, main= "A", suggestiveline = FALSE, genomewideline = FALSE)
abline(h=(FDR5), col= "blue", lwd=1)
abline(h=(FDR10), col= "blue", lwd=1)
###############################YW MEAN B ######

ywbBLINK <-read.csv('~/Desktop/GRIN/YWonlyColorLab/GAPIT.Association.GWAS_Results.BLINK.MeanB.csv',head=TRUE,sep=',',stringsAsFactors=F)
qq(ywbBLINK$P.value, main = "Q-Q plot of GWAS p-values")
ywbBLINK$FDR=p.adjust(ywbBLINK$P.value, method="fdr")
colnames(ywbBLINK)[2] <- "CHR"
colnames(ywbBLINK)[3] <- "BP"
colnames(ywbBLINK)[4] <- "P"
colnames(ywbBLINK)[9] <- "FDR"

ywbBLINK <- ywbBLINK[,c(1:4,9)]
ywbBLINK_snpsOfInterest <- ywbBLINK[ywbBLINK$FDR <  0.1, "SNP"]
print(ywbBLINK_snpsOfInterest)
FDR5=-log10(min(subset(ywbBLINK,FDR>0.049)$P))
FDR10=-log10(min(subset(ywbBLINK,FDR>0.1)$P))

manhattan(ywbBLINK, highlight= ywbBLINK_snpsOfInterest, main= "B", suggestiveline = FALSE, genomewideline = FALSE) 
abline(h=(FDR5), col= "blue", lwd=1)
abline(h=(FDR10), col= "blue", lwd=1)
