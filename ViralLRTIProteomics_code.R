library(tidyverse)
library(dplyr)
library(limma)
library(ggplot2)
library(MSnSet.utils)
library(edgeR)
library(WebGestaltR)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(DOSE)
library(yardstick)
library(glmnet)
library(plotROC)
library(pROC)
library(ggpubr)
library(ggsignif)
library(WebGestaltR)

#reading in and manipulating necessary files for analysis
metadata <- read.csv("ViralLRTIProteomics_metadata.csv")
TAdf <- read.csv("TAdf.csv", check.names = F)
rownames(TAdf) <- TAdf[,1]
TAdf <- TAdf[,-1]
plasmadf <- read.csv("plasmadf.csv", check.names = F)
rownames(plasmadf) <- plasmadf[,1]
plasmadf <- plasmadf[,-1]
protein_list <- colnames(TAdf)

#################################################################################

###Tracheal aspirate differential expression###
###Generates all plots in Figure 2, plus table S1 and figure S1###

#combining metadata and protein data
TAdfwithmetadata = cbind(Age = NA, LRTIstatus = NA, LRTIsubtype = NA, ViralRPMSum = NA, TAdf)
for (x in 1:nrow(TAdfwithmetadata))
{
  rowIndex = match(as.character(rownames(TAdfwithmetadata)[x]), as.character(metadata$SubjectID))
  TAdfwithmetadata$Age[x] = metadata$Age[rowIndex]
  TAdfwithmetadata$LRTIstatus[x] = metadata$LRTIStatus[rowIndex]
  TAdfwithmetadata$LRTIsubtype[x] = metadata$LRTIsubtype[rowIndex]
  TAdfwithmetadata$ViralRPMSum[x] = metadata$Viral_RPM_sum[rowIndex]
}

#limma to get list of DE proteins (no age adjustment)
transposedTAdf <- t(TAdfwithmetadata) #transposes to have columns as samples, proteins as rows (necessary for limma)
group = factor(transposedTAdf[2,]) #generating design matrix
group = factor(group, levels = c("NoLRTI", "LRTI")) #setting NoLRTI as the baseline comparison in design matrix
design = model.matrix(~group)
transposedTAdf <- transposedTAdf[-c(1:4),] #gets rid of metadata to get df ready for input in limma
class(transposedTAdf) = "numeric" #sets entire matrix to be numeric values, rather than strings, for limma
model_TA_unadjusted = lmFit(transposedTAdf, design) #limma
model_TA_unadjusted = eBayes(model_TA_unadjusted)
summary(decideTests(model_TA_unadjusted)) #outputs table of DE proteins

#generating list of all proteins with logFC/P value/T statistic, labeled by whether significantly upregulated or downregulated
model_TA_unadjusted_results = topTable(model_TA_unadjusted, coef = 2, number = 1305, adjust.method = "BH")
model_TA_unadjusted_results <- model_TA_unadjusted_results %>% mutate(gene_type= case_when(logFC >0 & adj.P.Val <= 0.05 ~ "Upregulated",
                                                           logFC <=0 & adj.P.Val <= 0.05 ~ "Downregulated",
                                                           TRUE ~ "Not significant")) ##adds categories for coloring volcano plot appropriately
model_TA_unadjusted_results <- tibble::rownames_to_column(model_TA_unadjusted_results, "Protein")

#generating volcano plot (Fig 2A)
cols = c("Upregulated" = "red", "Downregulated" = "blue", "Not significant" = "grey")
proteins_to_label = model_TA_unadjusted_results %>% top_n(-10, adj.P.Val) #labeling top 10 DE proteins
vol_plot1 <- model_TA_unadjusted_results %>%
  ggplot(aes(x=logFC,
             y = -log10(adj.P.Val),
             col = gene_type)) +
  geom_point() +
  scale_color_manual(values=cols) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_label_repel(max.overlaps = Inf, data = proteins_to_label, aes(label = Protein), size = 3, fill = "white") +
  xlab("Log(fold change)") + 
  ylab("-Log10(adjusted P value)") +
  theme_bw() +
  theme(legend.position = "none")
vol_plot1

#generating heat map (Fig 2B)
top20_DEproteins <- model_TA_unadjusted_results %>% top_n(-20, adj.P.Val)
#sig_proteins <- significant_results$Protein
rows_to_keep <- which(row.names(transposedTAdf) %in% top20_DEproteins$Protein)
heatmap_data <- as.data.frame(transposedTAdf[rows_to_keep,]) #creating new protein dataframe with just the top 20 DE proteins

#adding tstat corresponding to protein temporarily for re-ordering purposes (will then be removed)
heatmap_data$tstat = NA
for(x in 1:nrow(heatmap_data)) {
  index = match(row.names(heatmap_data)[x], model_TA_unadjusted_results$Protein)
  heatmap_data$tstat[x] <- model_TA_unadjusted_results$t[index]
}
heatmap_data <- heatmap_data[order(heatmap_data$tstat),] #reordering accordingly
heatmap_data <- as.matrix(heatmap_data[,-56]) #gets rid of tstatistic and converts to matrix

#generating side color vector indicating whether LRTI or no LRTI (will go above the columns)
heatmap_LRTI_labels <- as.data.frame(colnames(heatmap_data)) #extracting subject IDs
heatmap_LRTI_labels$LRTIstatus <- NA
heatmap_LRTI_labels$colors <- "blue" #setting default to no LRTI, which is blue (red = LRTI)
for(x in 1:nrow(heatmap_LRTI_labels)) {
  index = match(heatmap_LRTI_labels$`colnames(heatmap_data)`[x], metadata$SubjectID)
  heatmap_LRTI_labels$LRTIstatus[x] <- metadata$LRTIStatus[index]
  if (metadata$LRTIStatus[index] == "LRTI") {
    heatmap_LRTI_labels$colors[x] = "red"
  }
}

#code to generate heatmap, using complementary colors orange and purple to represent down- and up-regulation, respectively
colMain <- colorRampPalette(c("#e08214", "#f7f7f7", "#542788"))(25)
heatmap(heatmap_data, Rowv = NA, col=colMain, ColSideColors = heatmap_LRTI_labels$colors)
legend(legend =  c("vLRTI", "No LRTI"), x="topright", inset=c(-0.1, -0.1), fill = c("red", "blue"), box.lty = 0, bg = "transparent")

#performing enrichment analysis using WebGestaltR package: https://bzhanglab.github.io/WebGestaltR/ (Figure 3c)
#using uniprotein key to use UniProt IDs, which are needed for enrichment. 
#Correlation of protein aptamer name to UniProt ID provided by SomaScan (attached in needed files)
uniProteinIds <- read.csv("Uniprotein key.csv")
model_TA_unadjusted_results <- cbind(model_TA_unadjusted_results, "UniProt" = NA)
for (x in 1:nrow(model_TA_unadjusted_results)) {
  index = match(model_TA_unadjusted_results$Protein[x], uniProteinIds$Protein.Name)
  model_TA_unadjusted_results$UniProt[x] <- uniProteinIds$Uni.Protein.ID[index]
}

#generating rank file with UniProt ID, logFC, ranked by T statistic
#need to export then reimport to keep the .rnk extension (needed for WebGestalt function to run properly)
model_TA_unadjusted_results <- model_TA_unadjusted_results[order(model_TA_unadjusted_results$t, decreasing = TRUE), ] #ordering df based on T statistic
GSEA_input <- cbind(model_TA_unadjusted_results$UniProt, model_TA_unadjusted_results$logFC)
write.table(GSEA_input, file = "GSEA input.rnk", quote = F, sep="\t", row.names = F, col.names = F)
rankFilePath <- "GSEA input.rnk"

#running enrichment - note that results will be slightly different each time, as WebGestalt
#uses random permutation for significance evaluation, though the normalized enrichment score should be the same
#note: Webgestalt package will save a folder in your working directory with more details re: enrichment analysis
enrichment_result <- WebGestaltR(
  enrichMethod = "GSEA",
  organism = "hsapiens",
  enrichDatabase = "pathway_Reactome",
  interestGeneFile = rankFilePath,
  interestGeneType = "uniprotswissprot",
  sigMethod = "top",
  topThr = 10,
  minNum = 5
)

#generating enrichment plot (Figure 2C)
enrichment_result <- head(enrichment_result, 10)
enrichment_plot <- ggplot(enrichment_result, aes(x=normalizedEnrichmentScore, y = reorder(description, normalizedEnrichmentScore), size=size, color = FDR)) + 
  geom_point() + scale_size(name = "Count", range = c(1,10)) + scale_color_gradient(low="red", high = "blue") + xlim(-3, 3) +
  theme_bw() + xlab("Normalized Enrichment Score") +ylab(" ") + theme(axis.text.y = element_text(size=12))# + theme(axis.text.x = element_text(size = 14))
enrichment_plot

##Model generation for tracheal aspirate D1 LRTI vs no LRTI classifier 
##using 5 fold cross-validation within cohort to generate performance characteristics
#confidence intervals generated with bootstrapping
#Generating figure 2D

#first, need to generate the folds. will randomly select folds to keep similar breakdown (n=3 no LRTI cases at least)
outcome_and_fold <- factor(TAdfwithmetadata$LRTIstatus, levels = c("NoLRTI", "LRTI"))
outcome_and_fold <- as.data.frame(outcome_and_fold)
outcome_and_fold <- cbind(rownames(TAdfwithmetadata), outcome_and_fold)
colnames(outcome_and_fold) <- c("SubjectID", "LRTIStatus")

#while loop to interate through to generate folds 
success <- FALSE
while (!success) {
  outcome_and_fold %>% 
    dplyr::mutate(fold = sample(rep(1:5, length.out = nrow(.)))) -> cv_folds
  
  cv_folds %>% dplyr::select(LRTIStatus, fold) %>% 
    table() -> cv_fold_table
  
  cv_fold_table %>%
    .["NoLRTI", ] %>%
    min() %>%
    {. > 2} ->
    success
  
  if (!success) {
    print("At least one fold has too few No Evidence samples. Regenerating folds...")
  }
}

#input files for model generation - adding fold to proteomic dataframe (temporarily) 
#generating empty predictions file
classifier_input <-cbind(cv_folds$fold, TAdf)
colnames(classifier_input)[1] <- "fold"
predictions <- as.data.frame(cbind("SubjectID" = NA, "LRTIStatus" = NA, fold = NA, pred = NA))

#for loop for each fold to generate model on ~80% of the subjects and predict on the held-out 20% of the subjects
#uses the function glmnet for both model coefficient selection and prediction (table S1)
for (x in 1:5) {
  test_dataset <- classifier_input %>% filter(fold == x)
  test_LRTIstatus <- cv_folds %>% filter(fold == x)
  train_dataset <- classifier_input %>% filter(fold != x)
  train_LRTIstatus <- cv_folds %>% filter(fold != x)
  model = glmnet::cv.glmnet(as.matrix(train_dataset[,-1]), train_LRTIstatus$LRTIStatus, family = "binomial")
  pred = predict(model, as.matrix(test_dataset[, -1]), type='response',
                 s='lambda.1se', gamma=c("gamma.1se"))[,1]
  predictions_fold_metadata = cbind(test_LRTIstatus, pred)
  predictions <- rbind(predictions, predictions_fold_metadata) #adding the new predictions for each fold
}

#if desired, can extract the coefficients from each of the folds using this code (used in table S1)
coefficients <- as.data.frame(coef(model, s='lambda.1se', gamma=c("gamma.1se"))[,1] %>%
                                # Filter to nonzero coefficients
                                .[. != 0])

##for entire dataset, generating a complete model to see coefficient weights (used in table S1)
complete_model <- cv.glmnet(as.matrix(classifier_input[,-1]), outcome_and_fold$LRTIStatus, family = "binomial")
all_coeffs <- as.data.frame(coef(complete_model, s='lambda.1se', gamma=c("gamma.1se"))[,1] %>%
                                # Filter to nonzero coefficients
                                .[. != 0])

#using predictions to generate a ROC object using R package ROC
predictions <- predictions[-1, ] #gets rid of the first row - all NAs
predictions$LRTIStatus <- factor(predictions$LRTIStatus, levels =c("NoLRTI", "LRTI"))
ROC_curve <- roc(predictions$LRTIStatus ~ as.numeric(predictions$pred), ci=TRUE)
ROC_curve$ci
ROC_curve$auc

##generating CI data using bootstrapping
ciobj <- ci.se(ROC_curve, specificities = seq(0,1,l=25))
dat.ci <- data.frame(x=as.numeric(rownames(ciobj)),
                     lower = ciobj[,1],
                     upper = ciobj[,3])

#plotting ROC curve (Figure 2D)
ROCplot_formatted <- ggroc(ROC_curve, col = "red", alpha = 0.7) +
  geom_ribbon(data = dat.ci, aes(x=x, ymin = lower, ymax = upper), fill = "red", alpha = 0.2) +
  theme_bw() +
  labs(x="False positive rate (1 - specificity)", y = "True positive rate (sensitivity)") +
  annotate("label", x=0.5, y=0.75, label = paste("AUC = ", round(ROC_curve$auc, 2), "\n95% CI", round(ROC_curve$ci[1], 2), "-", round(ROC_curve$ci[2], 2)))
ROCplot_formatted

##repeating limma on tracheal aspirate vLRTI vs no LRTI comparison, but adding age as a continuous variable into the model
##generates Figure S1
transposedTAdf <- t(TAdfwithmetadata) #transposes to have columns as samples, proteins as rows
group = factor(transposedTAdf[2,]) #generating design matrix with LRTI status
group = factor(group, levels = c("NoLRTI", "LRTI")) #setting NoLRTI as the baseline comparison
age = transposedTAdf[1,] #extracting age
class(age) = "numeric" #changing age from string to numeric variable
design = model.matrix(~group+age) #adding age to design matrix
transposedTAdf <- transposedTAdf[-c(1:4),] #gets rid of metadata that you needed to generate design matrix, leaving just protein values
class(transposedTAdf) = "numeric" #forces all values to numeric, previously strings

#Generating the age-adjusted model
model_TA_adjusted = lmFit(transposedTAdf, design)
model_TA_adjusted = eBayes(model_TA_adjusted)
summary(decideTests(model_TA_adjusted)) #outputs table of differentially expressed proteins

#extracting DE results from the model and annotating for volcano plot
model_TA_adjusted_results = topTable(model_TA_adjusted, coef = 2, number = 1305, adjust.method = "BH")
model_TA_adjusted_results <- model_TA_adjusted_results %>% mutate(gene_type= case_when(logFC >0 & adj.P.Val <= 0.05 ~ "Upregulated",
                                                           logFC <=0 & adj.P.Val <= 0.05 ~ "Downregulated",
                                                           TRUE ~ "Not significant")) ##adds categories for coloring volcano plot appropriately
model_TA_adjusted_results <- tibble::rownames_to_column(model_TA_adjusted_results, "Protein")

#generating volcano plot for age-adjusted proteins (Figure S1)
cols = c("Upregulated" = "red", "Downregulated" = "blue", "Not significant" = "grey")
proteins_to_label = model_TA_adjusted_results %>% top_n(-10, adj.P.Val)
vol_plot2 <- model_TA_adjusted_results %>%
  ggplot(aes(x=logFC,
             y = -log10(adj.P.Val),
             col = gene_type)) +
  geom_point() +
  scale_color_manual(values=cols) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_label_repel(max.overlaps = Inf, data = proteins_to_label, aes(label = Protein), size = 3, fill = "white") +
  xlab("Log(fold change)") + 
  ylab("-Log10(adjusted P value)") +
  theme_bw() +
  theme(legend.position = "none")
vol_plot2

#################################################################################

###Plasma differential expression###
###Generates all plots in Figure 2, plus table S1 and figure S1###

#combining metadata and protein data
plasmadfwithmetadata = cbind(Age = NA, LRTIstatus = NA, LRTIsubtype = NA, ViralRPMSum = NA, plasmadf)
for (x in 1:nrow(plasmadfwithmetadata))
{
  rowIndex = match(as.character(rownames(plasmadfwithmetadata)[x]), as.character(metadata$SubjectID))
  plasmadfwithmetadata$Age[x] = metadata$Age[rowIndex]
  plasmadfwithmetadata$LRTIstatus[x] = metadata$LRTIStatus[rowIndex]
  plasmadfwithmetadata$LRTIsubtype[x] = metadata$LRTIsubtype[rowIndex]
  plasmadfwithmetadata$ViralRPMSum[x] = metadata$Viral_RPM_sum[rowIndex]
}

#limma to get list of DE proteins in plasma (no age adjustment)
plasma_transposeddf <- t(plasmadfwithmetadata) #transposes to have columns as samples, proteins as rows (necessary for limma)
group = factor(plasma_transposeddf[2,]) #generating design matrix
group = factor(group, levels = c("NoLRTI", "LRTI")) #setting NoLRTI as the baseline comparison in design matrix
design = model.matrix(~group)
plasma_transposeddf <- plasma_transposeddf[-c(1:4),] #gets rid of metadata to get df ready for input in limma
class(plasma_transposeddf) = "numeric" #sets entire matrix to be numeric values, rather than strings, for limma
model_plasma_unadjusted = lmFit(plasma_transposeddf, design) #limma
model_plasma_unadjusted = eBayes(model_plasma_unadjusted)
summary(decideTests(model_plasma_unadjusted)) #outputs table of DE proteins

#generating list of all proteins with logFC/P value/T statistic, labeled by whether significantly upregulated or downregulated
model_plasma_unadjusted_results = topTable(model_plasma_unadjusted, coef = 2, number = 1305, adjust.method = "BH")
model_plasma_unadjusted_results <- model_plasma_unadjusted_results %>% mutate(gene_type= case_when(logFC >0 & adj.P.Val <= 0.05 ~ "Upregulated",
                                                                                           logFC <=0 & adj.P.Val <= 0.05 ~ "Downregulated",
                                                                                           TRUE ~ "Not significant")) ##adds categories for coloring volcano plot appropriately
model_plasma_unadjusted_results <- tibble::rownames_to_column(model_plasma_unadjusted_results, "Protein")

#generating volcano plot for DE proteins in plasma, age unadjusted (Fig 3A)
cols = c("Upregulated" = "red", "Downregulated" = "blue", "Not significant" = "grey")
proteins_to_label = model_plasma_unadjusted_results %>% top_n(-10, adj.P.Val) #labeling top 10 DE proteins
vol_plot3 <- model_plasma_unadjusted_results %>%
  ggplot(aes(x=logFC,
             y = -log10(adj.P.Val),
             col = gene_type)) +
  geom_point() +
  scale_color_manual(values=cols) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_label_repel(max.overlaps = Inf, data = proteins_to_label, aes(label = Protein), size = 3, fill = "white") +
  xlab("Log(fold change)") + 
  ylab("-Log10(adjusted P value)") +
  theme_bw() +
  theme(legend.position = "none")
vol_plot3

##repeating limma on plasma vLRTI vs no LRTI comparison, but adding age as a continuous variable into the model
##generates Figure S2
plasma_transposeddf <- t(plasmadfwithmetadata) #transposes to have columns as samples, proteins as rows (necessary for limma)
group = factor(plasma_transposeddf[2,]) #generating design matrix
group = factor(group, levels = c("NoLRTI", "LRTI")) #setting NoLRTI as the baseline comparison in design matrix
age = plasma_transposeddf[1,] #extracting age
class(age) = "numeric" #changing age from string to numeric variable
design = model.matrix(~group+age) #adding age to design matrix
plasma_transposeddf <- plasma_transposeddf[-c(1:4),] #gets rid of metadata that you needed to generate design matrix, leaving just protein values
class(plasma_transposeddf) = "numeric" #forces all values to numeric, previously strings

#Generating the age-adjusted model
model_plasma_adjusted = lmFit(plasma_transposeddf, design)
model_plasma_adjusted = eBayes(model_plasma_adjusted)
summary(decideTests(model_plasma_adjusted)) #outputs table of differentially expressed proteins

#extracting DE results from the model and annotating for volcano plot
model_plasma_adjusted_results = topTable(model_plasma_adjusted, coef = 2, number = 1305, adjust.method = "BH")
model_plasma_adjusted_results <- model_plasma_adjusted_results %>% mutate(gene_type= case_when(logFC >0 & adj.P.Val <= 0.05 ~ "Upregulated",
                                                                                       logFC <=0 & adj.P.Val <= 0.05 ~ "Downregulated",
                                                                                       TRUE ~ "Not significant")) ##adds categories for coloring volcano plot appropriately
model_plasma_adjusted_results <- tibble::rownames_to_column(model_plasma_adjusted_results, "Protein")

#generating volcano plot for age-adjusted proteins (Figure S2)
cols = c("Upregulated" = "red", "Downregulated" = "blue", "Not significant" = "grey")
proteins_to_label = model_plasma_adjusted_results %>% top_n(-1, adj.P.Val)
vol_plot4 <- model_plasma_adjusted_results %>%
  ggplot(aes(x=logFC,
             y = -log10(adj.P.Val),
             col = gene_type)) +
  geom_point() +
  scale_color_manual(values=cols) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_label_repel(max.overlaps = Inf, data = proteins_to_label, aes(label = Protein), size = 3, fill = "white") +
  xlab("Log(fold change)") + 
  ylab("-Log10(adjusted P value)") +
  theme_bw() +
  theme(legend.position = "none")
vol_plot4

#################################################################################

###Comparative analyses between plasma and tracheal aspirate###

#Looking at ISG-15 first, since this is the only DE protein in plasma when using age adjustment 
##and one of the most significantly DE in trach aspirate

#Generates box plot for ISG-15 (figure 3B)
ISG15_TA <- as.data.frame(cbind(TAdfwithmetadata$LRTIstatus, sampletype = "TA", as.numeric(TAdfwithmetadata$ISG15)))
ISG15_plasma <- cbind(plasmadfwithmetadata$LRTIstatus, sampletype = "Plasma", as.numeric(plasmadfwithmetadata$ISG15))
ISG15_all <- rbind(ISG15_TA, ISG15_plasma) #generates dataframe of all ISG-15 values, labeled by whether they are TA or plasma and whether they are vLRTI or No LRTI
colnames(ISG15_all) <- c("LRTI_Status", "Sample_Type", "ISG15")
ISG15_all <-as.data.frame(ISG15_all)
ISG15_all$ISG15 <- as.numeric(ISG15_all$ISG15)
ISG15_all$LRTI_Status <- factor(ISG15_all$LRTI_Status, levels = c("NoLRTI", "LRTI"))
ISG15_boxplot <- ggplot(ISG15_all, aes(x=Sample_Type, y=ISG15, color = LRTI_Status)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge()) +
  scale_color_manual(values=c("blue", "red"), labels = c('No LRTI', 'vLRTI')) +
  stat_compare_means(aes(group = LRTI_Status), label = "p.format", label.y=17.5) +
  theme_bw() +
  labs(y="ISG15 (log2(RFU)", x="Sample Type") +
  ylim(c(8, 18.7)) +
  theme(legend.position= c(0.5, 0.94), legend.direction = "horizontal", legend.title = element_blank(), legend.background = element_rect(color = "darkgrey", fill = "white", linetype = "solid"))
ISG15_boxplot

#Generating ROC curves for tracheal aspirate ISG15 then plasma ISG15 (generates figure S3)
colnames(ISG15_TA) <- c("LRTI_Status", "Sample_Type", "ISG15")
ISG15_TA$LRTI_Status <- factor(ISG15_TA$LRTI_Status, levels = c("NoLRTI", "LRTI"))
ROC_ISG15_TA <- roc(ISG15_TA$LRTI_Status ~ as.numeric(ISG15_TA$ISG15), ci=TRUE)
ROC_ISG15_TA$ci
ROC_ISG15_TA$auc

ROCplot_ISG15_TA <- ggroc(ROC_ISG15_TA, col = "red", alpha = 1) +
  theme_bw() +
  labs(x="False positive rate (1 - specificity)", y = "True positive rate (sensitivity)") +
  annotate("label", x=0.5, y=0.75, label = paste("AUC = ", round(ROC_ISG15_TA$auc, 2), "\n95% CI", round(ROC_ISG15_TA$ci[1], 2), "-", round(ROC_ISG15_TA$ci[3], 2))) +
  ggtitle("Tracheal Aspirate ISG15") +
  theme(plot.title = element_text(hjust = 0.5))
ROCplot_ISG15_TA

ISG15_plasma <- ISG15_all %>% filter(Sample_Type == "Plasma")
ROC_ISG15_plasma <- roc(ISG15_plasma$LRTI_Status ~ as.numeric(ISG15_plasma$ISG15), ci=TRUE)
ROC_ISG15_plasma$ci
ROC_ISG15_plasma$auc

ROCplot_ISG15_plasma <- ggroc(ROC_ISG15_plasma, col = "red", alpha = 1) +
  theme_bw() +
  labs(x="False positive rate (1 - specificity)", y = "True positive rate (sensitivity)") +
  annotate("label", x=0.5, y=0.75, label = paste("AUC = ", round(ROC_ISG15_plasma$auc, 2), "\n95% CI", round(ROC_ISG15_plasma$ci[1], 2), "-", round(ROC_ISG15_plasma$ci[3], 2))) +
  ggtitle("Plasma ISG15") +
  theme(plot.title = element_text(hjust = 0.5))
ROCplot_ISG15_plasma

#Comparing paired protein values between trach aspirate and plasma
#First generating a table with all paired values for all proteins (takes some time to run)
log_log_table <- data.frame(matrix(nrow=0, ncol = 5))
for (x in 1:nrow(TAdfwithmetadata)) {
  index = NA
  index <- match(rownames(TAdfwithmetadata)[x], rownames(plasmadfwithmetadata))
  if (is.na(index) == FALSE) { #aka if we find a matching subject that has both TA and plasma data
    for (y in 5:ncol(TAdfwithmetadata)) {
      log_log_table <- rbind(log_log_table, c(rownames(TAdfwithmetadata)[x], colnames(TAdfwithmetadata)[y], as.numeric(TAdfwithmetadata[x, y]), as.numeric(plasmadfwithmetadata[index,y]), TAdfwithmetadata$LRTIstatus[x]))
    }
  }
}
colnames(log_log_table) <- c("subject", "protein", "TA_value", "plasma_value", "LRTI_status")

#Generating a log-log/correlation plot for ISG15 first (Figure 3C)
log_log_table_ISG15 <- log_log_table %>% filter(protein == "ISG15")
log_log_table_ISG15$TA_value <- as.numeric(log_log_table_ISG15$TA_value)
log_log_table_ISG15$plasma_value <- as.numeric(log_log_table_ISG15$plasma_value)
log_log_plot_ISG15 <- ggplot(log_log_table_ISG15, aes(x=TA_value, y=plasma_value)) +
  geom_point(size=1) +
  labs(x="Log2(ISG15 RFUs), plasma", y="Log2(ISG15 RFUs), TA") +
  geom_smooth(method="lm", se=TRUE) +
  coord_cartesian(ylim = c(8,16), xlim = c(8,16)) +
  theme_bw() +
  geom_label(x=9.5, y=11.9, label = paste("Cor =", round(cor.test(log_log_table_ISG15$TA_value, log_log_table_ISG15$plasma_value)$estimate, 2), "\np =", signif(cor.test(log_log_table_ISG15$TA_value, log_log_table_ISG15$plasma_value)$p.value, digits = 4)))
log_log_plot_ISG15

#Generating correlation coefficients for all proteins, in sum and stratified by LRTI status
#Data for density plot (Figure 3E)
protein_correlations_all <- data.frame(protein = protein_list, corr_all = "NA", pvalue_all = "NA", corr_LRTI = "NA", pvalue_LRTI = "NA", corr_noLRTI = "NA", pvalue_noLRTI = "NA")
for (x in 1:nrow(protein_correlations_all)) { ##adding correlations for total dataset of paired samples first
  log_log_protein_table <- log_log_table %>% filter(protein == protein_correlations_all$protein[x])
  log_log_protein_table$TA_value <- as.numeric(log_log_protein_table$TA_value)
  log_log_protein_table$plasma_value <- as.numeric(log_log_protein_table$plasma_value)
  corrtest <- cor.test(log_log_protein_table$plasma_value, log_log_protein_table$TA_value)
  protein_correlations_all$corr_all[x] <- round(corrtest$estimate, digits = 2)
  protein_correlations_all$pvalue_all[x] <- round(corrtest$p.value, digits = 5)
}
protein_correlations_all$corr_all = as.numeric(protein_correlations_all$corr_all)
protein_correlations_all$pvalue_all = as.numeric(protein_correlations_all$pvalue_all)

#adding correlations for LRTI paired samples
log_log_table_LRTI <- log_log_table %>% filter(LRTI_status == "LRTI")
for (x in 1:nrow(protein_correlations_all)) {
  log_log_protein_table <- log_log_table_LRTI %>% filter(protein == protein_correlations_all$protein[x])
  log_log_protein_table$TA_value <- as.numeric(log_log_protein_table$TA_value)
  log_log_protein_table$plasma_value <- as.numeric(log_log_protein_table$plasma_value)
  corrtest <- cor.test(log_log_protein_table$plasma_value, log_log_protein_table$TA_value)
  protein_correlations_all$corr_LRTI[x] <- round(corrtest$estimate, digits = 2)
  protein_correlations_all$pvalue_LRTI[x] <- round(corrtest$p.value, digits = 5)
}
protein_correlations_all$corr_LRTI = as.numeric(protein_correlations_all$corr_LRTI)
protein_correlations_all$pvalue_LRTI = as.numeric(protein_correlations_all$pvalue_LRTI)

#adding correlations for noLRTI paired samples
log_log_table_noLRTI <- log_log_table %>% filter(LRTI_status == "NoLRTI")
for (x in 1:nrow(protein_correlations_all)) {
  log_log_protein_table <- log_log_table_noLRTI %>% filter(protein == protein_correlations_all$protein[x])
  log_log_protein_table$TA_value <- as.numeric(log_log_protein_table$TA_value)
  log_log_protein_table$plasma_value <- as.numeric(log_log_protein_table$plasma_value)
  corrtest <- cor.test(log_log_protein_table$plasma_value, log_log_protein_table$TA_value)
  protein_correlations_all$corr_noLRTI[x] <- round(corrtest$estimate, digits = 2)
  protein_correlations_all$pvalue_noLRTI[x] <- round(corrtest$p.value, digits = 5)
}
protein_correlations_all$corr_noLRTI = as.numeric(protein_correlations_all$corr_noLRTI)
protein_correlations_all$pvalue_noLRTI = as.numeric(protein_correlations_all$pvalue_noLRTI)

#Generating density plot. setting dashed lines to indicate typical cut points for weak, moderate, and strong correlation.
colors <- c("vLRTI" = "red", "No LRTI" = "blue")
correlation_density_plot <- protein_correlations_all %>%
  ggplot() +
  geom_density(aes(x=corr_LRTI, fill = "vLRTI"), alpha = 0.9) +
  geom_density(aes(x=corr_noLRTI, fill = "No LRTI"), alpha = 0.7) +
  scale_fill_manual(values=colors) +
  labs(x="Pearson's correlation values across paired samples", y="Protein Density") +
  ylim(0, 2.2) +
  xlim (-0.91, 0.91) +
  theme_bw() +
  geom_vline(xintercept=-0.5, color = "grey", linetype = "dashed") +
  geom_vline(xintercept=-0.3, color = "grey", linetype = "dashed") +
  geom_vline(xintercept=0.3, color = "grey", linetype = "dashed") +
  geom_vline(xintercept=0.5, color = "grey", linetype = "dashed") +
  theme(legend.position= c(0.88, 0.88), legend.title = element_blank(), legend.background = element_rect(color = "darkgrey", fill = "white", linetype = "solid"))
correlation_density_plot

#Comparison of DE proteins in tracheal aspirate to plasma (age unadjusted)
#Generates Figure 3D
combined_TA_plasma <- merge(model_TA_unadjusted_results[,-c(9,10)], model_plasma_unadjusted_results, by="Protein")
combined_TA_plasma <- combined_TA_plasma %>% rename(gene_type.x = "protein type trach aspirate", gene_type.y = "protein type plasma")
combined_TA_plasma <- combined_TA_plasma %>% mutate(protein_type= case_when(`protein type trach aspirate`== "Upregulated" & `protein type plasma` == "Upregulated" ~ "Both Upregulated",
                                                                            `protein type trach aspirate`== "Downregulated" & `protein type plasma` == "Downregulated" ~ "Both Downregulated",
                                                                            `protein type trach aspirate`== "Upregulated" & `protein type plasma` == "Downregulated" ~ "Opposite",
                                                                            `protein type trach aspirate`== "Downregulated" & `protein type plasma` == "Upregulated" ~ "Opposite",
                                                                            TRUE ~ "Not significant"))

#scatter plot comparing TA and plasma differentially expressed proteins (plotting T statistic by T statistic)
#intercepts chosen incidate the corresponding T statistic that equates to adjusted p-value of 0.05
cols = c("Both Upregulated" = "red", "Both Downregulated" = "blue", "Opposite" = "blueviolet", "Not significant" = "grey")
proteins_to_label <- combined_TA_plasma %>% filter(protein_type != "Not significant")
DEproteins_TAvsplasma_plot <- combined_TA_plasma %>%
  ggplot(aes(x=t.x, y=t.y, col = protein_type)) + 
  geom_vline(xintercept=-2.75, color = "grey") +
  geom_vline(xintercept=2.75, color = "grey") +
  geom_hline(yintercept=-3.21, color = "grey") +
  geom_hline(yintercept=3.21, color = "grey")+
  geom_point(size = 0.5) +
  scale_color_manual(values=cols) +
  theme_bw() +
  geom_label_repel(max.overlaps = Inf, data = proteins_to_label, aes(label = Protein), size = 4, fill = "white") +
  labs(x="Tracheal Aspirate T-statistic", y="Plasma T-statistic") +
  theme(legend.position = "none") +
  xlim(-10, 10) +
  ylim(-10, 10)
DEproteins_TAvsplasma_plot


#################################################################################

###Subanalysis within vLRTI to identify DE proteins between viral infection and bacterial-viral co-infection###
###Within tracheal aspirate samples###

TA_LRTI <- TAdfwithmetadata %>% filter(LRTIstatus == "LRTI")
transposed_TA_LRTI <- t(TA_LRTI) #transposes to have columns as samples, proteins as rows
#colnames(TA_LRTI_d1) <- TA_LRTI_d1[1,]
group = factor(transposed_TA_LRTI[3,]) #generating design matrix
group = factor(group, levels = c("Viral", "Coinfection")) #setting viral as the baseline comparison
design = model.matrix(~group)
transposed_TA_LRTI <- transposed_TA_LRTI[-c(1:4),] #gets rid of metadata that you needed to generate design matrix
class(transposed_TA_LRTI) = "numeric" #previously was string

#generating the model
model_TA_coinfection = lmFit(transposed_TA_LRTI, design)
model_TA_coinfection = eBayes(model_TA_coinfection)
summary(decideTests(model_TA_coinfection))

#generating a table of results from the model
model_TA_coinfection_results = topTable(model_TA_coinfection, coef = 2, number = 1305, adjust.method = "BH")
model_TA_coinfection_results <- tibble::rownames_to_column(model_TA_coinfection_results, "Protein")

#adding the uniProt ID for performing enrichment analysis
model_TA_coinfection_results <- cbind(model_TA_coinfection_results, "UniProt" = NA)
for (x in 1:nrow(model_TA_coinfection_results)) {
  index = match(model_TA_coinfection_results$Protein[x], uniProteinIds$Protein.Name)
  model_TA_coinfection_results$UniProt[x] <- uniProteinIds$Uni.Protein.ID[index]
}

#generating rank file with UniProt ID, logFC, ranked by T statistic
#need to export then reimport to keep the .rnk extension (needed for WebGestalt function to run properly)
model_TA_coinfection_results <- model_TA_coinfection_results[order(model_TA_coinfection_results$t, decreasing = TRUE), ] #ordering df based on T statistic
GSEA_input_coinfection <- cbind(model_TA_coinfection_results$UniProt, model_TA_coinfection_results$logFC)
write.table(GSEA_input_coinfection, file = "GSEA input_coinfection.rnk", quote = F, sep="\t", row.names = F, col.names = F)
rankFilePath_coinfection <- "GSEA input_coinfection.rnk"

#running enrichment - note that results will be slightly different each time, as WebGestalt
#uses random permutation for significance evaluation, though the normalized enrichment score should be the same
#note: Webgestalt package will save a folder in your working directory with more details re: enrichment analysis
enrichment_result_coinfection <- WebGestaltR(
  enrichMethod = "GSEA",
  organism = "hsapiens",
  enrichDatabase = "pathway_Reactome",
  interestGeneFile = rankFilePath_coinfection,
  interestGeneType = "uniprotswissprot",
  sigMethod = "top",
  topThr = 10,
  minNum = 5
)

#plotting coinfection enrichment results (Figure 4B)
enrichment_plot_coinfection <- ggplot(enrichment_result_coinfection, aes(x=normalizedEnrichmentScore, y = reorder(description, normalizedEnrichmentScore), size=size, color = FDR)) + 
  geom_point() + scale_size(name = "Count", range = c(1,10)) + scale_color_gradient(low="purple", high = "orange") + xlim(-3, 3) +
  theme_bw() + xlab("Normalized Enrichment Score") +ylab(" ") + theme(axis.text.y = element_text(size=12)) +
  geom_vline(xintercept = 0)
enrichment_plot_coinfection

##generating box plots for TSG6 and CRP (two of the most DE proteins, upregulated in coinfection, based on limma)
#TSG-6 first (Figure 4A)
TSG6 <- cbind(TA_LRTI$LRTIsubtype, TA_LRTI$`TSG-6`)
colnames(TSG6) <- c("Infection_Status", "TSG6")
TSG6 <-as.data.frame(TSG6)
TSG6$TSG6<- as.numeric(TSG6$TSG6)
TSG6$Infection_Status <- factor(TSG6$Infection_Status, levels = c("Viral", "Coinfection"))
my_comparisons <- list(c("Viral", "Coinfection"))
TSG6_boxplot <- ggplot(TSG6, aes(x=Infection_Status, y=TSG6, color = Infection_Status)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge()) +
  scale_color_manual(values=c("purple", "orange")) +
  theme_bw() +
  ylim(8, 13.4) +
  geom_signif(comparisons = my_comparisons, map_signif_level = function(x) paste("Adj. p =", round(model_TA_coinfection_results[1,6], 2)), color = "black") +
  labs(y="TSG-6 (log2(RFU)", x="Infection Type", title="TSG-6") +
  theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank())
TSG6_boxplot

#then CRP (Figure 4A)
CRP <- cbind(TA_LRTI$LRTIsubtype, TA_LRTI$CRP)
colnames(CRP) <- c("Infection_Status", "CRP")
CRP <-as.data.frame(CRP)
CRP$CRP<- as.numeric(CRP$CRP)
CRP$Infection_Status <- factor(CRP$Infection_Status, levels = c("Viral", "Coinfection"))
CRP_boxplot <- ggplot(CRP, aes(x=Infection_Status, y=CRP, color = Infection_Status)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge()) +
  scale_color_manual(values=c("purple", "orange")) +
  theme_bw() +
  ylim(9.8, 18) +
  geom_signif(comparisons = my_comparisons, map_signif_level = function(x) paste("Adj. p =", round(model_TA_coinfection_results[2,6], 2)), color = "black") +
  labs(y="CRP (log2(RFU)", x="Infection Type", title="CRP") +
  theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank())
CRP_boxplot

#################################################################################

###Subanalysis within vLRTI (subset that were positive for virus on mNGS) to investigate
###correlation of specific proteins with viral load

#Correlation calculations for log(viral RPM sum) and log2(protein concentrations)
TAdf_viralmNGSpos <- TAdfwithmetadata %>% filter(ViralRPMSum != 0 & LRTIstatus == "LRTI") #filtering for subjects with LRTI with viral mNGS+
protein_RPM_correlations <- data.frame(protein = protein_list, corr = "NA", pvalue = "NA")
for (x in 1:nrow(protein_RPM_correlations)) { ##cycling through all proteins, adding correlations with viral load
  corrtest <- cor.test(TAdf_viralmNGSpos[,4+x], log(TAdf_viralmNGSpos$ViralRPMSum))
  protein_RPM_correlations$corr[x] <- round(corrtest$estimate, digits = 2)
  protein_RPM_correlations$pvalue[x] <- round(corrtest$p.value, digits = 5)
}
protein_RPM_correlations$corr = as.numeric(protein_RPM_correlations$corr)
protein_RPM_correlations$pvalue = as.numeric(protein_RPM_correlations$pvalue)

#now selecting a few of the biologically relevant more and least correlated proteins
#ISG15 (positively correlated)
ISG15_viral_scatter <- ggplot(TAdf_viralmNGSpos, aes(x=log(ViralRPMSum), y=ISG15)) +
  geom_point(size=1) +
  geom_smooth(method="lm", se=TRUE, color = "indianred3", fill = "indianred3") +
  theme_bw() +
  geom_label(x=3.1, y=12.5, label = paste("Cor =", round(cor.test(log(TAdf_viralmNGSpos$ViralRPMSum), TAdf_viralmNGSpos$ISG15)$estimate, 2),"\np = ", round(cor.test(log(TAdf_viralmNGSpos$ViralRPMSum), TAdf_viralmNGSpos$ISG15)$p.value, 3))) +
  labs(x = "Log(total viral load)", y = "Log2(ISG15)", title = "ISG15") +
  theme(plot.title = element_text(hjust = 0.5)) 
ISG15_viral_scatter

#IFN-lambda 1 (positively correlated)
IFNlambda1_viral_scatter <- ggplot(TAdf_viralmNGSpos, aes(x=log(ViralRPMSum), y=`IFN-lambda 1`)) +
  geom_point(size=1) +
  geom_smooth(method="lm", se=TRUE, color = "indianred3", fill = "indianred3") +
  theme_bw() +
  geom_label(x=2.9, y=8, label = paste("Cor =", round(cor.test(log(TAdf_viralmNGSpos$ViralRPMSum), TAdf_viralmNGSpos$`IFN-lambda 1`)$estimate, 2),"\np = ", round(cor.test(log(TAdf_viralmNGSpos$ViralRPMSum), TAdf_viralmNGSpos$`IFN-lambda 1`)$p.value, 6))) +
  labs(x = "Log(Total Viral Load)", y = "Log2(IFN-lambda 1 RFUs)", title = "IFN-lambda 1") +
  theme(plot.title = element_text(hjust = 0.5)) 
IFNlambda1_viral_scatter

#MCP2 (positively correlated)
MCP2_viral_scatter <- ggplot(TAdf_viralmNGSpos, aes(x=log(ViralRPMSum), y=`MCP-2`)) +
  geom_point(size=1) +
  geom_smooth(method="lm", se=TRUE, color = "indianred3", fill = "indianred3") +
  theme_bw() +
  geom_label(x=3, y=7.2, label = paste("Cor =", round(cor.test(log(TAdf_viralmNGSpos$ViralRPMSum), TAdf_viralmNGSpos$`MCP-2`)$estimate, 2),"\np = ", round(cor.test(log(TAdf_viralmNGSpos$ViralRPMSum), TAdf_viralmNGSpos$`MCP-2`)$p.value, 6))) +
  labs(x = "Log(Total Viral Load)", y = "Log2(MCP-2 RFUs)", title = "MCP-2") +
  theme(plot.title = element_text(hjust = 0.5)) 
MCP2_viral_scatter

#GI-24
GI24_viral_scatter <- ggplot(TAdf_viralmNGSpos, aes(x=log(ViralRPMSum), y=GI24)) +
  geom_point(size=1) +
  geom_smooth(method="lm", se=TRUE, color = "cornflowerblue", fill = "cornflowerblue") +
  theme_bw() +
  geom_label(x=10, y=14, label = paste("Cor =", round(cor.test(log(TAdf_viralmNGSpos$ViralRPMSum), TAdf_viralmNGSpos$GI24)$estimate, 2),"\np = ", round(cor.test(log(TAdf_viralmNGSpos$ViralRPMSum), TAdf_viralmNGSpos$GI24)$p.value, 6))) +
  labs(x = "Log(Total Viral Load)", y = "Log2(GI_24 RFUs)", title = "GI-24") +
  theme(plot.title = element_text(hjust = 0.5)) 
GI24_viral_scatter

#Activin AB
activin_viral_scatter <- ggplot(TAdf_viralmNGSpos, aes(x=log(ViralRPMSum), y=`Activin AB`)) +
  geom_point(size=1) +
  geom_smooth(method="lm", se=TRUE, color = "cornflowerblue", fill = "cornflowerblue") +
  theme_bw() +
  geom_label(x=10.2, y=10, label = paste("Cor =", round(cor.test(log(TAdf_viralmNGSpos$ViralRPMSum), TAdf_viralmNGSpos$`Activin AB`)$estimate, 2),"\np = ", round(cor.test(log(TAdf_viralmNGSpos$ViralRPMSum), TAdf_viralmNGSpos$`Activin AB`)$p.value, 6))) +
  labs(x = "Log(Total Viral Load)", y = "Log2(Activin AB RFUs)", title = "Activin AB") +
  theme(plot.title = element_text(hjust = 0.5)) 
activin_viral_scatter

#CD177
CD177_viral_scatter <- ggplot(TAdf_viralmNGSpos, aes(x=log(ViralRPMSum), y=CD177)) +
  geom_point(size=1) +
  geom_smooth(method="lm", se=TRUE, color = "cornflowerblue", fill = "cornflowerblue") +
  theme_bw() +
  geom_label(x=10.7, y=15.7, label = paste("Cor =", round(cor.test(log(TAdf_viralmNGSpos$ViralRPMSum), TAdf_viralmNGSpos$CD177)$estimate, 2),"\np = ", round(cor.test(log(TAdf_viralmNGSpos$ViralRPMSum), TAdf_viralmNGSpos$CD177)$p.value, 6))) +
  labs(x = "Log(Total Viral Load)", y = "Log2(CD177 RFUs)", title = "CD177") +
  theme(plot.title = element_text(hjust = 0.5)) 
CD177_viral_scatter

