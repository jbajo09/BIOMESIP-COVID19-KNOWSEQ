set.seed(60)

# LOAD PACKAGES
require(KnowSeq)
require(GEOquery)
require(caret)
require(M3C)
require(class)
require(ggplot2)

source('calculateGeneExpressionValues.R')

# LABELS FROM GEO 
GSE152075 <- getGEO("GSE152075", destdir = '')
labels <- GSE152075$GSE152075_series_matrix.txt.gz$`sars-cov-2 positivity:ch1`

# LOAD COUNTS MATRIX FROM GEO.
countsMatrix <- as.matrix(read.table('GSE152075_raw_counts_GEO.txt', header =TRUE))

# Get GC content  
Annotation_gene <- getGenesAnnotation(rownames(countsMatrix), filter = 'external_gene_name')

# FROM COUNTS TO GENE EXPRESSION
expressionMatrix <- calculateGeneExpressionValues(countsMatrix, Annotation_gene, genesNames = TRUE, Ensembl_ID = FALSE)

# REMOVE OUTLIERS
outliers <- RNAseqQA(expressionMatrix, toRemoval = TRUE, toPDF = FALSE, toPNG = FALSE) # 30 outliers
expressionMatrix_out <-  outliers$matrix
labels_out <- labels[-which(colnames(expressionMatrix)%in%outliers$outliers)]

# BATCH EFFECT TREATMENT WITH SVA
expressionMatrix_fix <- batchEffectRemoval(expressionMatrix_out,labels_out, method = 'sva')

# TRAINING-TEST SPLIT (80%-20%)
Index_train_test <- createDataPartition(labels_out, p = .80, list = FALSE, times = 1)
train_labels <- labels_out[Index_train_test]
test_labels <- labels_out[-Index_train_test]
train_matrix <- expressionMatrix_fix[,Index_train_test]
test_matrix <- expressionMatrix_fix[,-Index_train_test]

# DEGs extraction
DEGs <- DEGsExtraction(train_matrix, as.factor(train_labels), lfc=2, pvalue = 0.01, number = Inf, CV=TRUE,numFolds = 5)

# Feature selection by mRMR
mRMR_rank <- featureSelection(t(train_matrix),train_labels,DEGs$Common_DEGs, mode ='mrmr')

# Training k-NN CV
knn_CV <- knn_trn(t(train_matrix), as.factor(train_labels), names(mRMR_rank), 10) 
plot(knn_CV$accuracyInfo$meanAccuracy, type = 'l', col= 'black', ylab='Metric Performance', xlab='Genes', lwd=2.5, ylim = c(0.85,1), panel.first = grid(col='gray45'),cex.axis=1.4,cex.lab=1.6)
lines(knn_CV$sensitivityInfo$meanSensitivity, col='blue', lwd=2.5, lty=2)
lines(knn_CV$specificityInfo$meanSpecificity, col='#FF8B00', lwd=2.5, lty=4)
lines(knn_CV$F1Info$meanF1, col='red', lwd=2.5, lty=4)
legend(x=5.56 ,y =0.8915, c('Accuracy', 'Sensitivity','Specificity','F1-Score'), lty = c(1,2,4,5), col = c('black','blue','#FF8B00','red'), cex=1.6)

# Test 
knn_test <- knn_test(t(train_matrix), as.factor(train_labels), t(test_matrix), as.factor(test_labels), names(mRMR_rank), bestK =  knn_CV$bestK)

#cf Matrix
dataPlot(knn_test$cfMats[[4]]$table,test_labels,mode = 'confusionMatrix',toPNG = FALSE, toPDF = FALSE)

#T-SNE
tsne(expressionMatrix_fix[which(rownames(expressionMatrix_fix)%in%names(mRMR_rank)[1:4]),],labels=as.factor(labels_out),controlscale=TRUE, scale=3, colvec = c('blue','red'),seed = 3, axistextsize = 0)

#boxplot
dataPlot(expressionMatrix_fix[which(rownames(expressionMatrix_fix) %in% names(mRMR_rank[1:4])),],labels_out,mode = "genesBoxplot",toPNG = FALSE, colours = c("red", "blue"))


