setwd("/Volumes/HD3/NGS/G_RNASeq/Cufflinks")
getwd()
library(cummeRbund)
gf = "/Volumes/HD3/UCSC/mm10/Annotation/Archives/archive-2015-07-17-14-33-26/Genes/genes.gtf"
rm(cuff)
cuff = readCufflinks(gtfFile = gf, genome = 'mm10')
cuff
disp<-dispersionPlot(genes(cuff)) # evaluate the quality of the model fitting
disp
dispersionPlot(cuff) #directly will allow you to visualize the full model fit. 
genes.scv<-fpkmSCVPlot(genes(cuff)) #The squared coefficient of variation is a normalized measure of cross-replicate variability that can be useful for evaluating the quality your RNA-seq data. Differences in CV 2 can result in lower numbers of differentially expressed genes due to a higher degree of variability between replicate fpkm estimates.

isoforms.scv<-fpkmSCVPlot(isoforms(cuff)) 

dens<-csDensity(genes(cuff)) # To assess the distributions of FPKM scores across samples
dens

b<-csBoxplot(genes(cuff))
b
s<-csScatterMatrix(genes(cuff)) # matrix of pairwise scatterplots 
s
v<-csVolcanoMatrix(genes(cuff))
v
runInfo(cuff) #Run-level information such as run parameters, and sample information can be accessed from a CuffSet object by using the runInfo and replicates methods

#Retrive significant gene IDs (XLOC) with a pre-specified alpha

diffGeneIDs <- getSig(cuff,level="genes",alpha=0.05) 
#Use returned identifiers to create a CuffGeneSet object with all relevant info for given genes
diffGenes<-getGenes(cuff,diffGeneIDs)

#gene_short_name values (and corresponding XLOC_* values) can be retrieved from the CuffGeneSet by using:
featureNames(diffGenes)
head(diffGenes)
mySigMat<-sigMatrix(cuff,level='genes',alpha=0.05)
mySigMat

mySigGeneIds<-getSig(cuff,alpha=0.05,level='genes')
length(mySigGeneIds)
Y_MINvsO_MIN = getSig(cuff, x = "OLD_PLUS", y = "YOUNG_PLUS", alpha = 0.05, level = "isoforms")
diffY_MINvsO_MIN = getGenes(cuff, Y_MINvsO_MIN)
featureNames(diffY_MINvsO_MIN)
#NTvsIL = getSig(cuff, x= "NT", y = "IL4", alpha = 0.05, level = "genes")
myDistHeat<-csDistHeat(genes(cuff))
myDistHeat
#myGeneIds = c('Il6', 'Cdkn2a', 'Il4', 'Il1b', 'Mmp2', 'Mmp3', 'Trp53', 'Stat3', 'Stat6', 'Glb1', 'Cdkn2b', 'Il8')
myGeneIds = c("Eef1a1","Cox6c","Eif3l","Rpl12","Gtf2h1","Tmem258","Ctsa","Anxa2","Pla2g7","Thbs3","Hdlbp","Atp6v1b2","Nfe2l1","Cirh1a","Dnajc7","Gnptg","Ddx6","Ssh2")
myGenes<-getGenes(cuff,myGeneIds)
head(fpkm(myGenes))
h<-csHeatmap(myGenes,cluster='both')
h
b<-expressionBarplot(myGenes)
b



feature.level <- "genes" # create variable for the script below
# ... or "isoforms", "TSS", "CDS"
idColumnName <- "gene_id" # create variable for the script below
# ... or "isoform_id", "TSS_group_id", "CDS_id" to match above
report.name <- paste0('DiffExp_', feature.level,'.txt') # create variable for the script below
# ... or whatever you want

#cuff <- readCufflinks()
sigIDs <- getSig(cuff,level=feature.level,alpha=0.05)
if (NROW(sigIDs) > 0) {
        sigFeatures <- getFeatures(cuff,sigIDs,level=feature.level)
        sigData <- diffData(sigFeatures)
        sigData <- subset(sigData, (significant == 'yes'))
        names <- featureNames(sigFeatures)
        sigOutput <- merge(names, sigData, by.x="tracking_id", 
                           by.y=idColumnName)
        
        # Patch the merged table to have the original name for the ID column.  
        # This is always the first column for the examples we've seen.
        colnames(sigOutput)[1] <- idColumnName
        write.table(sigOutput, report.name, sep='\t', row.names = F, 
                    col.names = T, quote = F)
}
nameID=names$gene_short_name # just to extract gene names
nameID #print it
write.csv(nameID, file='nameID.csv',  eol = "\r", col.names = TRUE)
