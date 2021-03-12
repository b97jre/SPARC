#----------------------------------------------# basics

setwd("~/Desktop/SPARC_Rcode_variation/")
library(genefilter)

own_residuals <- function(x, y, log=F, p=1, own_weights=F, return_labels=F, labels=c())
{
  filter = is.na(x) | is.na(y) | is.infinite(x) | is.infinite(y) | (x==0); x=x[!filter]; y=y[!filter]; if (own_weights!=F){own_weights=own_weights[!filter]}
  if (log==T){x = log2(x); y=log2(y)}
  if (p==1){if (own_weights==F){fit = lm(y ~ x)} else {fit = lm(y ~ x, weights=own_weights)}}
  if (p==2){if (own_weights==F){fit = lm(y ~ I(x^0) + I(x^1) + I(x^2))} else {fit = lm(y ~ I(x^0) + I(x^1) + I(x^2), weights=own_weights)}}
  if (p==3){if (own_weights==F){fit = lm(y ~ I(x^0) + I(x^1) + I(x^2) + I(x^3))} else {fit = lm(y ~ I(x^0) + I(x^1) + I(x^2) + I(x^3), weights=own_weights)}}
  if (p==4){if (own_weights==F){fit = lm(y ~ I(x^0) + I(x^1) + I(x^2) + I(x^3) + I(x^4))} else {fit = lm(y ~ I(x^0) + I(x^1) + I(x^2) + I(x^3) + I(x^4), weights=own_weights)}}
  if (p==5){if (own_weights==F){fit = lm(y ~ I(x^0) + I(x^1) + I(x^2) + I(x^3) + I(x^4) + I(x^5))} else {fit = lm(y ~ I(x^0) + I(x^1) + I(x^2) + I(x^3) + I(x^4) + I(x^5), weights=own_weights)}}
  if (p==6){if (own_weights==F){fit = lm(y ~ I(x^0) + I(x^1) + I(x^2) + I(x^3) + I(x^4) + I(x^5) + I(x^6))} else {fit = lm(y ~ I(x^0) + I(x^1) + I(x^2) + I(x^3) + I(x^4) + I(x^5) + I(x^6), weights=own_weights)}}
  if (p==7){if (own_weights==F){fit = lm(y ~ I(x^0) + I(x^1) + I(x^2) + I(x^3) + I(x^4) + I(x^5) + I(x^6) + I(x^7))} else {fit = lm(y ~ I(x^0) + I(x^1) + I(x^2) + I(x^3) + I(x^4) + I(x^5) + I(x^6) + I(x^7), weights=own_weights)}}
  if (!return_labels)
  {
    if (log==T){return(list(2**x, residuals(fit)))}
    else{return(list(x, residuals(fit)))}
  }
  else
  {
    names(x) = labels[!filter]
    z = residuals(fit)
    names(z) = labels[!filter]
    if (log==T){return(list(2**x, z))}
    else{return(list(x, z))}    
  }
}

least_total_squares <- function(x, y)
{
  filter = (!is.na(x)) & (!is.na(y)); x = x[filter]; y = y[filter]
  vector <- prcomp(cbind(x, y))$rotation
  beta <- (-vector[(-ncol(vector)),ncol(vector)] / vector[ncol(vector),ncol(vector)])
  inter <- mean(y) - beta * mean(x)
  return(c(inter, beta))
}

#----------------------------------------------# loading

load(file='UU_RNA_all.RData') # RNA_all
load(file='UU_cell_info.RData') # cell_info
load(file='UU_gene_info.RData') # gene_info

load(file='UU_protein_all.RData') # protein_all
load(file='UU_cell_info_protein.RData') # cell_info_protein
load(file='UU_gene_info_protein.RData') # gene_info_protein

#----------------------------------------------# filtering

shared_cells = intersect(colnames(RNA_all), colnames(protein_all))
shared_genes = intersect(rownames(RNA_all), rownames(protein_all))

# RNA_all = t(t(RNA_all)/colSums(RNA_all))*(10**6) # old, simple RPM

RNA_all <- data.frame(read.table('RPKMs.table.QC.egn.tsv'))

length(rownames(RNA_all)); length(unique(rownames(RNA_all)))
RNA_all = RNA_all[rowSums(RNA_all)>0,]

#----------------------------------------------# time selection

filtered_cells_0h = shared_cells[(cell_info[shared_cells,'pseudotime_scorpius']>=0.8) & (cell_info_protein[shared_cells,'pseudotime_scorpius']<1.1) & (cell_info_protein[shared_cells,'time_given']=='0h')]
filtered_cells_24h = shared_cells[(cell_info[shared_cells,'pseudotime_scorpius']>=0.45) & (cell_info_protein[shared_cells,'pseudotime_scorpius']<0.8) & (cell_info_protein[shared_cells,'time_given']=='24h')]
filtered_cells_48h = shared_cells[(cell_info[shared_cells,'pseudotime_scorpius']>=0.0) & (cell_info_protein[shared_cells,'pseudotime_scorpius']<0.55) & (cell_info_protein[shared_cells,'time_given']=='48h')]
filtered_cells_72h = rownames(cell_info_protein)[(cell_info_protein[,'time_given']=='72h')]
filtered_cells_all = c(filtered_cells_0h, filtered_cells_24h, filtered_cells_48h)

RNA_0h_filtered = RNA_all[,filtered_cells_0h]
RNA_24h_filtered = RNA_all[,filtered_cells_24h]
RNA_48h_filtered = RNA_all[,filtered_cells_48h]

protein_0h_filtered = protein_all[,filtered_cells_0h]
protein_24h_filtered = protein_all[,filtered_cells_24h]
protein_48h_filtered = protein_all[,filtered_cells_48h]
protein_72h_filtered = protein_all[,filtered_cells_72h]

#----------------------------------------------# further filtering (low protein abundance, cell cycle, etc.)

data_RNA = RNA_0h_filtered[,intersect(colnames(RNA_0h_filtered), rownames(cell_info))]

# FACS data
protein_FACS = read.table(file='scProtein_0_72h_ESC_QC.txt')
protein_FACS = t(protein_FACS); protein_FACS = protein_FACS[,protein_FACS['type',]=='FACS']; protein_FACS = protein_FACS[5:(dim(protein_FACS)[1]),]; class(protein_FACS) <- "numeric" 
protein_FACS_0h = protein_FACS[,substr(colnames(protein_FACS),1,7)=='Comb_0h']

data_protein = protein_0h_filtered[,intersect(colnames(protein_0h_filtered), rownames(cell_info))] # looks better

data_protein = data_protein[rowMeans(data_protein)!=0,]

low_abundance_filter = intersect(rownames(data_protein), rownames(protein_FACS_0h)[rank(rowMeans(protein_FACS_0h))>6])
data_protein = data_protein[low_abundance_filter,]

data_RNA = data_RNA[rownames(data_protein),]

#----------------------------------------------# gene expression noise

data_median_RNA = rowMeans(data_RNA)
data_variation_RNA = (rowSds(as.matrix(data_RNA))/rowMeans(data_RNA))**2
data_residuals_RNA = own_residuals(data_median_RNA, data_variation_RNA, log=T, p=3)[[2]]

data_median_protein = rowMeans(data_protein)
data_variation_protein = (rowSds(data_protein)/rowMeans(data_protein))**2
data_residuals_protein = own_residuals(data_median_protein, data_variation_protein, log=T, p=3)[[2]]

shared_genes = intersect(names(data_residuals_RNA), names(data_residuals_protein))

load(file='PaxDB_ES.RData')

shared_genes_PaxDB = intersect(shared_genes, rownames(PaxDB_ES))

noise_ratio = (data_residuals_protein[shared_genes_PaxDB]) - (data_residuals_RNA[shared_genes_PaxDB]) - median(data_residuals_protein[shared_genes_PaxDB]) + median(data_residuals_RNA[shared_genes_PaxDB])

expression_ratio = log2(PaxDB_ES[shared_genes_PaxDB,4] / rowMeans(data_RNA[shared_genes_PaxDB,]))

function_to_optimize <- function(data, par){-cor.test((par[1]*data[,1])+(par[2]*data[,2]), data[,3], method='spearman')$estimate}
data_for_optimization = cbind(data_residuals_RNA[shared_genes_PaxDB], expression_ratio[shared_genes_PaxDB], data_residuals_protein[shared_genes_PaxDB])
para = optim(c(1,1), function_to_optimize, data=data_for_optimization, method="SANN")$par; para # used to be SANN

#----------------------------------------------# set up figure

par(mfrow=c(2,2), cex.axis=1.5, cex.main=1.5, cex.lab=1.5)

library(scales); rbPal <- colorRampPalette(c('cornflowerblue','grey','coral'))

dec_col <- rbPal(100)[as.numeric(cut(rank(data_residuals_RNA), breaks = 100))]
plot(log(data_median_RNA)[names(data_residuals_RNA)], log(data_variation_RNA)[names(data_residuals_RNA)], col=dec_col, pch=16, frame.plot=F, xlab='Mean RNA expression (RPKM)', ylab='RNA expression noise (CV^2)', main='Gene expression and gene expression noise', cex=1.5)
points(log(data_median_RNA)[names(data_residuals_RNA)], log(data_variation_RNA)[names(data_residuals_RNA)], cex=1.5, lwd=0.5)
#abline(lts[1], lts[2], lty='dotted')
text(log(data_median_RNA)[names(data_residuals_RNA)], log(data_variation_RNA)[names(data_residuals_RNA)], names(data_residuals_RNA), pos=1)
legend("bottomleft", inset=0.05, legend=c("highly variable", "stably expressed"), fill=c("coral", "cornflowerblue"), cex=1)

dec_col <- rbPal(100)[as.numeric(cut(rank(data_residuals_protein), breaks = 100))]
plot(log(data_median_protein)[names(data_residuals_protein)], log(data_variation_protein)[names(data_residuals_protein)], col=dec_col, pch=16, frame.plot=F, xlab='Mean protein expression (2^dCq)', ylab='protein expression noise (CV^2)', main='Gene expression and gene expression noise', cex=1.5)
points(log(data_median_protein)[names(data_residuals_protein)], log(data_variation_protein)[names(data_residuals_protein)], cex=1.5, lwd=0.5)
#abline(lts[1], lts[2], lty='dotted')
text(log(data_median_protein)[names(data_residuals_protein)], log(data_variation_protein)[names(data_residuals_protein)], names(data_residuals_protein), pos=1)
legend("bottomleft", inset=0.05, legend=c("highly variable", "stably expressed"), fill=c("coral", "cornflowerblue"), cex=1)

shared_genes = intersect(names(data_residuals_RNA), names(data_residuals_protein))
rho = round(cor.test(data_residuals_RNA[shared_genes], data_residuals_protein[shared_genes], method='spearman')$estimate,3)
rho2 = round(cor.test(data_residuals_RNA[shared_genes_PaxDB], data_residuals_protein[shared_genes_PaxDB], method='spearman')$estimate,3)
plot(data_residuals_RNA[shared_genes], data_residuals_protein[shared_genes], frame.plot=F, pch=16, xlab='normalized RNA expression variation (A.U.)', ylab='normalized protein expression variation (A.U.)', main=paste('Gene expression variation\nrho_1: ',rho,"; rho_2:",rho2,sep=''),  cex=1.5, col="grey")
points(data_residuals_RNA[shared_genes], data_residuals_protein[shared_genes], cex=1.5, lwd=0.5)
points(data_residuals_RNA[shared_genes_PaxDB], data_residuals_protein[shared_genes_PaxDB], cex=1.5, pch=16)
lts = least_total_squares(data_residuals_RNA[shared_genes], data_residuals_protein[shared_genes]); inter = lts[1]; slope = lts[2]
abline(inter, slope, lty='dotted')
text(data_residuals_RNA[shared_genes], data_residuals_protein[shared_genes], shared_genes, pos=1)

shared_genes_PaxDB = intersect(shared_genes, rownames(PaxDB_ES))
expression_ratio = log2(PaxDB_ES[shared_genes_PaxDB,4] / rowMeans(data_RNA[shared_genes_PaxDB,]))

shared_genes = intersect(intersect(names(data_residuals_RNA), names(data_residuals_protein)), names(expression_ratio))
rho = round(cor.test(data_residuals_RNA[shared_genes]+para[2]/para[1]*expression_ratio[shared_genes], data_residuals_protein[shared_genes], method='spearman')$estimate,3)

plot(data_residuals_RNA[shared_genes]+para[2]/para[1]*expression_ratio[shared_genes], data_residuals_protein[shared_genes], pch=16, cex=1.5, frame.plot=F, col="grey", xlab='Combined RNA normalized variation and translation rate (A.U.)', ylab='Normalized protein expression variation (A.U.)', main=paste('Translation of gene expression variation\nrho: ',rho,sep=''))#, ylim=c(-1.5,1.5))#, xlim=c(-0.5,4.5))
points(data_residuals_RNA[shared_genes]+para[2]/para[1]*expression_ratio[shared_genes], data_residuals_protein[shared_genes], cex=1.5, lwd=0.5, pch=16)
lts = least_total_squares(data_residuals_RNA[shared_genes]+para[2]/para[1]*expression_ratio[shared_genes], data_residuals_protein[shared_genes]); inter = lts[1]; slope = lts[2]
abline(inter, slope, lty='dotted')
text(data_residuals_RNA[shared_genes]+para[2]/para[1]*expression_ratio[shared_genes],  data_residuals_protein[shared_genes], shared_genes, pos=1)
