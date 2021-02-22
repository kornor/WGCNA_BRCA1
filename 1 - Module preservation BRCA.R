
## SCript for analysing some other expression data to look for expression modules as in the TCGA / Metabric set

#set for the BRCA folder
setwd("~/Bioinformatics Work/Meth & RNA/BRCA1 data");
# Load the package
library(WGCNA);
library(flashClust)
library(dplyr)
library(tidyr)
# String settings
options(stringsAsFactors = FALSE);

# load the metabric / tcga stuff

load(file = "MetaAnalysis_trimmed_input.RData")
load(file = "Modules_DS0.RData")
colorsA1 = names(table(modules1))

#datExpr1 = TCGA
#datExpre2 = Metabric

## Load in the PDX data
brca_e <- read.table("BRCA_exp1.txt", sep = "\t", header = TRUE, row.names = 1, check.names = FALSE) 

######## initial load had duplicate rownames:
## removed these
#brca1 <- brca_e[!duplicated(brca_e$Gene),]
# then would not run goodgenes because of "non-numeric value"
# initially I thought this was NAs, but turns out to have been one "+" sign


### preprocess the data using goodgenes (WCGNA)

gsg = goodSamplesGenes(brca_e,verbose = 5);
gsg$allOK

## if return is true, all good to go
### Otherwise

if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(brca_e)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(brca_e)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  brca1 = brca_e[gsg$goodSamples, gsg$goodGenes]
}



########## check that gene names match 
##brca1 now the active expression file
commongenes <- intersect(rownames(datExpr1), rownames(brca1))
brca1 <- brca1[commongenes,]

#transpose again?
datExpr1_t <- as.data.frame(t(datExpr1))
brca_t <- as.data.frame(t(brca1))

######### Calculation of module preservation

setLabels = c("TCGA","BRCA");
multiExpr = list(TCGA= list(data = datExpr1_t), BRCA = list(data = brca_t));
multiColor = list(TCGA = modules1);

### check for presevation

system.time( {
  mp = modulePreservation(multiExpr, multiColor,
                          referenceNetworks = 1,
                          nPermutations = 200,
                          randomSeed = 1,
                          quickCor = 0,
                          verbose = 3)
} );
# Save the results
save(mp, file = "modulePreservation_BRCA.RData");

ref = 1
test = 2
statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1]);

# Compare preservation to quality:
print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
             signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )

##graphing
# Module labels and module sizes are also contained in the results
modColors = rownames(mp$preservation$observed[[ref]][[test]])
moduleSizes = mp$preservation$Z[[ref]][[test]][, 1];

# leave grey and gold modules out
plotMods = !(modColors %in% c("grey", "gold"));
# Text labels for points
text = modColors[plotMods];
# Auxiliary convenience variable
plotData = cbind(mp$preservation$observed[[ref]][[test]][, 2], mp$preservation$Z[[ref]][[test]][, 2])
# Main titles for the plot
mains = c("Preservation Median rank", "Preservation Zsummary");
# Start the plot
sizeGrWindow(10, 5);
#pdf(fi="Plots/TCGA_PDX-modulePreservation-Zsummary-medianRank.pdf", wi=10, h=5)
par(mfrow = c(1,2))
par(mar = c(4.5,4.5,2.5,1))
for (p in 1:2)
{
  min = min(plotData[, p], na.rm = TRUE);
  max = max(plotData[, p], na.rm = TRUE);
  # Adjust ploting ranges appropriately
  if (p==2)
  {
    if (min > -max/10) min = -max/10
    ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
  } else
    ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
  plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
       main = mains[p],
       cex = 2.4,
       ylab = mains[p], xlab = "Module size", log = "x",
       ylim = ylim,
       xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
  labelPoints(moduleSizes[plotMods], plotData[plotMods, p], text, cex = 1, offs = 0.08);
  # For Zsummary, add threshold lines
  if (p==2)
  {
    abline(h=0)
    abline(h=2, col = "blue", lty = 2)
    abline(h=10, col = "darkgreen", lty = 2)
  }
}
# If plotting into a file, close it
dev.off();


#### now to calculate the module eigengene for each module as in TCGA:

## first need to trim so the genes in pdx_t are the same as in module colours (otherwise it will fail)

mod_gene <- as.data.frame(modules1)
row.names(mod_gene) <- row.names(datExpr1)

mod_brca <- mod_gene[commongenes,]


###now the calcululation:

PCs_brca <- moduleEigengenes(brca_t, 
                            mod_brca, 
                            impute = TRUE, 
                            nPC = 1, 
                            align = "along average", 
                            excludeGrey = FALSE, 
                            grey = if (is.numeric(mod_pdx)) 0 else "grey",
                            subHubs = TRUE,
                            softPower = 6,
                            scale = TRUE,
                            verbose = 0, indent = 0)

####works!!

ME_brca    = PCs_brca$eigengenes
distPC_brca = 1-abs(cor(ME_brca,use="p"))
distPC_brca = ifelse(is.na(distPC_brca), 0, distPC_brca)
pcTree_brca = hclust(as.dist(distPC_brca),method="a") 
MDS_brca  = cmdscale(as.dist(distPC_brca),2)
colorsBRCA = names(table(mod_brca))



### add in the rownames

rownames(ME_brca) <- colnames(brca_e)

write.table(ME_brca, file = "ME_BRCA.txt", sep = "\t")
save(PCs_brca, file = 'BRCA_moduledata.Rdata')



##### MDS of the module expression in the PDX group:

par()
plot(MDS_brca, col= colorsBRCA,  main="MDS plot", cex=2, pch=19,
     xlab = "Principal component 1", ylab = "Principal component 2")


######### look at the samples clustering?

## Create a dendro of the samples to look for outliers

sampleTree = flashClust(dist(brca_t), method = "average");


# Plot the sample tree

sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, 
     main = "Sample clustering", 
     sub="", xlab="", 
     cex.lab = 1.5,
     cex.axis = 1.5, 
     cex.main = 2, labels = FALSE)

