library(DESeq2)

cnts = rbind(ORF_table, CUT_reads, NUT_reads, SUT_reads, XUT_reads)
cnts = round(cnts)
cnts = cnts[,-7]
cond <- factor(c('WT', 'WT', 'CCR4', 'DHH1', 'LSM1', 'RPB4', 'XRN1', 'XRN1'))
# object construction
dds <- DESeqDataSetFromMatrix(cnts, DataFrame(cond), ~ cond)
#sizeFactors(dds) = norm_factors
sizeFactors(dds) = norm_factors[-7]/exp(mean(log(norm_factors[-7])))
# standard analysis
dds <- DESeq(dds)

cond_mut = unique(cond[-c(1:2)])
                        
fc_mat = NULL
res_list = list()
i = 1
for(y in c('CCR4', 'DHH1', 'LSM1', 'RPB4', 'XRN1')){
  res <- results(dds, contrast = c('cond', y, 'WT'))
  res_list[[i]] = res
  i = i+1
  # moderated log2 fold changes
  # resLFC <- lfcShrink(dds, contrast = c('cond', 'WT', 'CCR4'))
  
  fc = res$stat
  fc_mat = cbind(fc_mat, fc)
  
  #write.table(x = xu_annot$name[order(fc, decreasing = T)], file = paste0(y, '_up.txt'), quote = F,
  #col.names = F, row.names = F)
  #write.table(x = xu_annot$name[order(fc, decreasing = F)], file = paste0(y, '_down.txt'), quote = F,
  #col.names = F, row.names = F)
}
colnames(fc_mat) = c('CCR4', 'DHH1', 'LSM1', 'RPB4', 'XRN1')
names(res_list) = colnames(fc_mat)
fc_mat = data.frame(fc_mat)

#hyper_overlap = function(vec_1, vec_2, cutoff){
#  l1 = length(intersect(which(vec_1 < cutoff), which(vec_2 < cutoff)))
#  l2 = length(intersect(which(vec_1 > cutoff), which(vec_2 < cutoff)))
#  l3 = length(intersect(which(vec_1 < cutoff), which(vec_2 > cutoff)))
#  l4 = length(intersect(which(vec_1 > cutoff), which(vec_2 > cutoff)))
#  fisher_test = fisher.test(matrix(c(l1, l2, l3, l4), ncol = 2))
#  return(fisher_test)
#}

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(cor(ORF_fc[,-5], method = 'spearman', use = 'complete.obs'), method="color", col=col(200),  
         type='full', order="hclust", 
         addCoef.col = "black", 
         tl.col="black", tl.srt=45, diag=T )

ORF_fc = fc_mat[1:4973,]
rownames(ORF_fc) = ORF_annot$commonName

CUT_fc = fc_mat[4974:5898,]
NUT_fc = fc_mat[5899:7424,]
SUT_fc = fc_mat[7425:8271,]
XUT_fc = fc_mat[8272:9929,]

i = 1
a1 = density(ORF_fc[,i], na.rm = T)
a2 = density(CUT_fc[,i], na.rm = T)
a3 = density(NUT_fc[,i], na.rm = T)
a4 = density(SUT_fc[,i], na.rm = T)
a5 = density(XUT_fc[,i], na.rm = T)
plot(a1, typ = 'l', xlim = c(-5,5), ylim = c(0,.5), lwd = 2,
     main = expression(paste('xrn1', Delta)), xlab = 'Standardized NET-seq FC',
     cex.axis = 2, cex.lab = 2, cex.main = 3)
lines(a2, col = 2, lwd = 2)
lines(a3, col = 3, lwd = 2)
lines(a4, col = 4, lwd = 2)
lines(a5, col = 5, lwd = 2)
legend(x = 'topright', legend = c('Genes', 'CUTs', 'NUTs', 'SUTs', 'XUTs'),
       col = 1:5, lwd = 2, cex = 2)
abline(v = 0, lty = 2, lwd = 2)


a1 = density(ORF_fc$CCR4, na.rm = T)
a2 = density(ORF_fc$DHH1, na.rm = T)
a3 = density(ORF_fc$LSM1, na.rm = T)
a4 = density(ORF_fc$RPB4, na.rm = T)
a5 = density(ORF_fc$XRN1, na.rm = T)

pdf('paper/DF2/Figures/main/other_fcs.pdf', height = 7, width = 7)
par(mar = c(5,5,5,3))
plot(a1, xlim = c(-5,5), xlab = 'Standardized Log2 NET-seq FC', main = 'Gene fold change distributions',
     col = 4, lwd = 2, ylim = c(0,.5), cex.lab = 2, cex.main = 2, cex.axis = 2)
lines(a2,col = 2, lwd = 2)
lines(a3,col = 3, lwd = 2)
#lines(a4,col = 4, lwd = 2)
#lines(a5,col = 5, lwd = 2)
abline(v = 0, col = 'black', lty = 2, lwd = 3)
#legend(x = 'topright', legend = c(expression(paste("ccr4", Delta)), expression(paste("dhh1", Delta)), expression(paste("lsm1", Delta))),
#       col = 1:3, lwd = 2, cex = 2)
legend(x = 'topright', legend = c(expression(paste("dhh1", Delta)), expression(paste("lsm1", Delta)), expression(paste("ccr4", Delta))),
       col = c(2,3,4), lwd = 2, cex = 2)
dev.off()

par(mfrow = c(2,3))
par(mar = c(5,5,5,3))
xrn1_rating_list = lapply(2:10, function(x) which(ORF_annot$name %in% xrn1_synth_ratings$ORF.name[xrn1_synth_ratings$responsiveness == x]))
boxplot(ORF_fc[xrn1_rating_list[[1]],1], ORF_fc[xrn1_rating_list[[2]],1], ORF_fc[xrn1_rating_list[[3]],1], 
        ORF_fc[xrn1_rating_list[[4]],1], ORF_fc[xrn1_rating_list[[5]],1], ORF_fc[xrn1_rating_list[[6]],1], 
        ORF_fc[xrn1_rating_list[[7]],1], ORF_fc[xrn1_rating_list[[8]],1], ORF_fc[xrn1_rating_list[[9]],1],
        ylim = c(-5,5), ylab = 'NET-seq FC', xlab = 'Xrn1 responsiveness', names = 2:10, main = expression(paste('ccr4', Delta)),
        cex.main = 3, cex.lab = 2, cex.axis = 2)
abline(h = 0, lwd = 2, lty = 2, col = 2, cex = 2)

boxplot(ORF_fc[xrn1_rating_list[[1]],2], ORF_fc[xrn1_rating_list[[2]],2], ORF_fc[xrn1_rating_list[[3]],2], 
        ORF_fc[xrn1_rating_list[[4]],2], ORF_fc[xrn1_rating_list[[5]],2], ORF_fc[xrn1_rating_list[[6]],2], 
        ORF_fc[xrn1_rating_list[[7]],2], ORF_fc[xrn1_rating_list[[8]],2], ORF_fc[xrn1_rating_list[[9]],2],
        ylim = c(-5,5), ylab = 'NET-seq FC', xlab = 'Xrn1 responsiveness', names = 2:10, main = expression(paste('dhh1', Delta)),
        cex.main = 3, cex.lab = 2, cex.axis = 2)
abline(h = 0, lwd = 2, lty = 2, col = 2)

boxplot(ORF_fc[xrn1_rating_list[[1]],3], ORF_fc[xrn1_rating_list[[2]],3], ORF_fc[xrn1_rating_list[[3]],3], 
        ORF_fc[xrn1_rating_list[[4]],3], ORF_fc[xrn1_rating_list[[5]],3], ORF_fc[xrn1_rating_list[[6]],3], 
        ORF_fc[xrn1_rating_list[[7]],3], ORF_fc[xrn1_rating_list[[8]],3], ORF_fc[xrn1_rating_list[[9]],3],
        ylim = c(-5,5), ylab = 'NET-seq FC', xlab = 'Xrn1 responsiveness', names = 2:10, main = expression(paste('lsm1', Delta)),
        cex.main = 3, cex.lab = 2, cex.axis = 2)
abline(h = 0, lwd = 2, lty = 2, col = 2)

boxplot(ORF_fc[xrn1_rating_list[[1]],4], ORF_fc[xrn1_rating_list[[2]],4], ORF_fc[xrn1_rating_list[[3]],4], 
        ORF_fc[xrn1_rating_list[[4]],4], ORF_fc[xrn1_rating_list[[5]],4], ORF_fc[xrn1_rating_list[[6]],4], 
        ORF_fc[xrn1_rating_list[[7]],4], ORF_fc[xrn1_rating_list[[8]],4], ORF_fc[xrn1_rating_list[[9]],4],
        ylim = c(-5,5), ylab = 'NET-seq FC', xlab = 'Xrn1 responsiveness', names = 2:10, main = expression(paste('rpb4', Delta)),
        cex.main = 3, cex.lab = 2, cex.axis = 2)
abline(h = 0, lwd = 2, lty = 2, col = 2)

boxplot(ORF_fc[xrn1_rating_list[[1]],5], ORF_fc[xrn1_rating_list[[2]],5], ORF_fc[xrn1_rating_list[[3]],5], 
        ORF_fc[xrn1_rating_list[[4]],5], ORF_fc[xrn1_rating_list[[5]],5], ORF_fc[xrn1_rating_list[[6]],5], 
        ORF_fc[xrn1_rating_list[[7]],5], ORF_fc[xrn1_rating_list[[8]],5], ORF_fc[xrn1_rating_list[[9]],5],
        ylim = c(-5,5), ylab = 'NET-seq FC', xlab = 'Xrn1 responsiveness', names = 2:10, main = expression(paste('xrn1', Delta)),
        cex.main = 3, cex.lab = 2, cex.axis = 2)
abline(h = 0, lwd = 2, lty = 2, col = 2)


boxplot(ORF_fc[SAGA_ind,1], ORF_fc[TFIID_ind,1], main = expression(paste('ccr4', Delta)), names = c('SAGA-dominated', 'TFIID-dominated'), xlab = 'Complex association', ylab = 'Log2 Pol II FC', ylim = c(-5,5))
boxplot(ORF_fc[SAGA_ind,2], ORF_fc[TFIID_ind,2], main = expression(paste('dhh1', Delta)), names = c('SAGA-dominated', 'TFIID-dominated'), xlab = 'Complex association', ylab = 'Log2 Pol II FC', ylim = c(-5,5))
boxplot(ORF_fc[SAGA_ind,3], ORF_fc[TFIID_ind,3], main = expression(paste('lsm1', Delta)), names = c('SAGA-dominated', 'TFIID-dominated'), xlab = 'Complex association', ylab = 'Log2 Pol II FC', ylim = c(-5,5))
boxplot(ORF_fc[SAGA_ind,4], ORF_fc[TFIID_ind,4], main = expression(paste('rpb4', Delta)), names = c('SAGA-dominated', 'TFIID-dominated'), xlab = 'Complex association', ylab = 'Log2 Pol II FC', ylim = c(-5,5))
boxplot(ORF_fc[SAGA_ind,5], ORF_fc[TFIID_ind,5], main = expression(paste('xrn1', Delta)), names = c('SAGA-dominated', 'TFIID-dominated'), xlab = 'Complex association', ylab = 'Log2 Pol II FC', ylim = c(-5,5))

q_list = lapply(1:5, function(x) cut(ORF_fc[,x], quantile(ORF_fc[,x], c(0:10)/10, na.rm = T), include.lowest = T, labels = F))

mat_list = list()
n = 1
for(i in 1:5){
  new_mat = matrix(0, nrow = 10, ncol = 10)
  for(j in c(1:4,6)){
    if(j<=i){
      
    } else {
      for(q1 in 1:10){
        for(q2 in 1:10){
          new_mat[q1,q2] = jaccard(which(q_list[[i]]==q1), which(q_list[[j]]==q2))
    }
      }
      mat_list[[n]] = new_mat
      n = n+1
      }
  }
}

sapply(mat_list, function(x) x[1,1])
sapply(mat_list, function(x) x[10,10])

y1 = our_sense_reads[,-c(1:4,6)]
y2 = their_sense_reads[,-c(1:4,6)]

y = merge(y1, y2)
rownames(y) = y$name
y = y[,-c(1,12)]

their_nf = estimateSizeFactorsForMatrix(y[,10:15])*13.06185/2.9090171
merged_nf = c(norm_factors, their_nf)
cnts = y_TSS_merged
cnts = round(cnts)
#cnts = cnts[,-7]
cond <- factor(c('WT', 'WT', 'CCR4', 'DHH1', 'LSM1', 'RPB4', 'SFP1', 'XRN1', 'XRN1',
                 'WT', 'RCO1', 'DST1', 'EAF3', 'SET1', 'SET2'))
# object construction
dds <- DESeqDataSetFromMatrix(cnts, DataFrame(cond), ~ cond)
sizeFactors(dds) = merged_nf / exp(mean(log(merged_nf)))
# standard analysis
dds <- DESeq(dds)

cond_mut = unique(cond[-c(1:2)])

fc_mat = NULL
res_list = list()
i = 1
for(y in c('CCR4', 'DHH1', 'LSM1', 'RPB4', 'SFP1', 'XRN1', 'RCO1', 'DST1', 'EAF3', 'SET1', 'SET2')){
  res <- results(dds, contrast = c('cond', y, 'WT'))
  res_list[[i]] = res
  i = i+1
  # moderated log2 fold changes
  # resLFC <- lfcShrink(dds, contrast = c('cond', 'WT', 'CCR4'))
  
  fc = res$stat
  fc_mat = cbind(fc_mat, fc)
  
  #write.table(x = xu_annot$name[order(fc, decreasing = T)], file = paste0(y, '_up.txt'), quote = F,
  #col.names = F, row.names = F)
  #write.table(x = xu_annot$name[order(fc, decreasing = F)], file = paste0(y, '_down.txt'), quote = F,
  #col.names = F, row.names = F)
}
colnames(fc_mat) = c('CCR4', 'DHH1', 'LSM1', 'RPB4', 'SFP1', 'XRN1', 'RCO1', 'DST1', 'EAF3', 'SET1', 'SET2')
names(res_list) = colnames(fc_mat)
fc_mat = data.frame(fc_mat)

fc_mat = fc_mat[,c(6,2,3,1,4,7:11)]
M = cor(fc_mat)
corrplot(M)

y1_TSS = cbind(colSums(WT_TSS[101:601,]),
colSums(WT_pilot_TSS[101:601,]),
colSums(CCR4_TSS[101:601,]),
colSums(DHH1_TSS[101:601,]),
colSums(LSM1_TSS[101:601,]),
colSums(RPB4_TSS[101:601,]),
colSums(SFP1_TSS[101:601,]),
colSums(XRN1_TSS[101:601,]),
colSums(XRN1_pilot_TSS[101:601,]))
y1_TSS = cbind.data.frame(ORF_annot$name, y1_TSS)
colnames(y1_TSS) = c('name', 'WT', 'WT_p', 'ccr4', 'dhh1', 'lsm1', 'rpb4', 'sfp1', 'xrn1', 'xrn1_p')

y2_TSS = cbind(colSums(WT_NC_TSS[101:601,]),
colSums(RCO1_TSS[101:601,]),
colSums(DST1_TSS[101:601,]),
colSums(EAF3_TSS[101:601,]),
colSums(SET1_TSS[101:601,]),
colSums(SET2_TSS[101:601,]))
y2_TSS = cbind.data.frame(ORF_annot_c$name, y2_TSS)
colnames(y2_TSS) = c('name', 'WT_NC', 'rco1', 'dst1', 'eaf3', 'set1', 'set2')

y_TSS_merged = merge(y1_TSS, y2_TSS)[,-1]

