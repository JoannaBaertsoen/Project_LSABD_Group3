#step 1: load the needed packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("edgeR")
library(edgeR)

#step 2: counts => matrix: gene_id as rowname so that the matrix only has numbers in it 
counts <- read.csv("GSE80655_counts (1).csv")
rownames(counts) <- counts$gene_id
counts$gene_id <- NULL
counts <- as.matrix(counts)

#step 3: aligning counts and metadata: so that they have the same collums 
meta <- read.csv("GSE80655_metadata (1).csv")
all(colnames(counts) == meta$description)
rownames(meta) <- meta$description
meta <- meta[colnames(counts), ]

#step 4: check the data for missing values 
sum(is.na(counts))
colSums(is.na(meta))
#tells us that there are 3 missing values => more specific in the collum of brain_pH, so we have to 
#identify these missing values with the following code:
meta[is.na(meta$brain_pH), c("description", "brain_region", "diagnosis", "donor_id")]

#step 5: the solution for the missing values: imputate with the median 
meta$brain_pH[is.na(meta$brain_pH)] <- median(meta$brain_pH, na.rm = TRUE)

#step 6: check for duplicates in the data
sum(duplicated(rownames(counts)))
sum(duplicated(colnames(counts)))
sum(duplicated(meta$GSM))
counts <- counts[!duplicated(rownames(counts)), ]
table(table(meta$donor_id))
#this tells us that there are no duplicates, but there are 7 people with 2 samples and 89 people with 3 samples

#step 7: cleaning the data and informing R about all the values 
meta$diagnosis <- relevel(factor(meta$diagnosis), ref = "Control")
meta$brain_region <- factor(meta$brain_region)
meta$sex <- factor(meta$sex)
meta$ethnicity <- factor(meta$ethnicity)
meta$donor_id <- factor(meta$donor_id)
str(meta)

#step 8: remove genes with a very low expression level
sum(rowSums(counts) == 0)
y <- DGEList(counts = counts, group = meta$diagnosis)
keep <- filterByExpr(y, group = meta$diagnosis)
table(keep)
y <- y[keep, , keep.lib.sizes = FALSE]
#after this 33914 genes were removed, we have 23991 genes left 

#step 9: normalisation, a lot of genes have a different range of reads (11,7 - 73,2 million)
#we can not compare them properly => we have to normalize the samples
library(edgeR)
y <- normLibSizes(y)
logcpm <- cpm(y, log = TRUE, prior.count = 2)
View(logcpm)
dim(logcpm)

#step 10: qualitycontrol: looking for outliers, checking librarysize
barplot(y$samples$lib.size/1e6, las = 2, ylab = "Library size (miljoen)")
boxplot(logcmp, las = 2, outline = FALSE)
plotMDS(y, col = as.numeric(meta$brain_region), pch = 16)
legend("topright", legend = levels(meta$brain_region), col = 1:3, pch = 16)





















#next part = splitting the data in the different psychiatric diseases and in the different affected brainregions
