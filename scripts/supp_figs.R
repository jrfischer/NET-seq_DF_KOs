library(corrplot)

#### Figure S1 ####
a = cbind.data.frame(ORF_annot$name, ORF_fc)
b = cbind.data.frame(ORF_annot_c$name, ORF_fc_c)
colnames(a)[1] = 'name'
colnames(b)[1] = 'name'
d = merge(a,b)

M2 = cor(d[,-1])
rownames(M2)<- c(":italic('ccr4') * Delta", ":italic('dhh1') * Delta", ":italic('lsm1') * Delta", ":italic('rpb4') * Delta", ":italic('xrn1') * Delta",
                 ":italic('rco1') * Delta", ":italic('dst1') * Delta", ":italic('eaf3') * Delta", ":italic('set1') * Delta", ":italic('set2') * Delta")
colnames(M2)<- rownames(M2)
corrplot(M2)

##### Figure S2 ####
a = cbind.data.frame(ORF_ids, DHH1_rnaseq$wt_dhh1_log2FoldChange, LSM1_rnaseq$wt_lsm1_log2FoldChange)
colnames(a) = c('name', 'Dhh1_RNA-seq_FC', 'Lsm1_RNA-seq_FC')

b = cbind.data.frame(as.character(ORF_annot$name), ORF_fc[,2:3])
colnames(b) = c('name', 'Dhh1_NET-seq_FC', 'Lsm1_NET-seq_FC')

d = merge(a,b, by = 'name')
d = d[!is.na(d[,2]),]
d = d[! d[,2] == Inf,]

pdf(file = 'netseq_vs_rnaseq_vs_gro.pdf', width = 10, height = 8)
par(mar = c(5,5,4,1))
par(mfrow = c(2,2))
plot(d[,2], d[,3], xlab  = expression(paste(italic('dhh1'), Delta)), ylab = expression(paste(italic('lsm1'), Delta)), main = 'NET-seq Log2 FC',
     col = rgb(0,0,0,.15), pch = 16, cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5, xlim = c(-4,4), ylim = c(-3,3))
text(x = -3, y = 2.8, paste('r =', round(cor(d[,2], d[,3], method = c('pearson')),2)), cex = 2)

plot(d[,4], d[,5], xlab  = expression(paste(italic('dhh1'), Delta)), ylab = expression(paste(italic('lsm1'), Delta)), main = 'RNA-seq Log2 FC',
     col = rgb(0,0,0,.15), pch = 16, cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5, xlim = c(-6,6), ylim = c(-5,5))
text(x = -4.75, y = 4.75, paste('r =', round(cor(d[,4], d[,5], use = 'complete.obs', method = c('pearson')),2)), cex = 2)

plot(NET_GRO_FC$RPB4_NET, NET_GRO_FC$XRN1_NET, xlab = expression(paste(italic('rpb4'), Delta)), ylab = expression(paste(italic('xrn1'), Delta)), main = 'NET-seq Log2 FC',
     col = rgb(0,0,0,.15), pch = 16, cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5, xlim = c(-6,6), ylim = c(-5,5))
text(x = -4.75, y = 4.75, paste('r =', round(cor(NET_GRO_FC$RPB4_NET, NET_GRO_FC$XRN1_NET, use = 'complete.obs', method = c('pearson')),2)), cex = 2)

plot(NET_GRO_FC$RPB4_GRO, NET_GRO_FC$XRN1_GRO, xlab = expression(paste(italic('rpb4'), Delta)), ylab = expression(paste(italic('xrn1'), Delta)), main = 'GRO Log2 FC',
     col = rgb(0,0,0,.15), pch = 16, cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5, xlim = c(-5,3), ylim = c(-5,2))
text(x = -4, y = 1.75, paste('r =', round(cor(NET_GRO_FC$RPB4_GRO, NET_GRO_FC$XRN1_GRO, use = 'complete.obs', method = c('pearson')),2)), cex = 2)
dev.off()

pdf('netseq_vs_rnaseq_vs_gro2.pdf', width = 10, height = 8)
par(mar = c(5,5,4,1))
par(mfrow = c(2,2))
plot(d[,2], d[,4], xlab  = 'RNA-seq Log2 FC', ylab = 'NET-seq Log2 FC', main = expression(paste(italic('dhh1'), Delta)),
     col = rgb(0,0,0,.15), pch = 16, cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5, xlim = c(-3,3), ylim = c(-5,5))
text(x = -2.25, y = 4.7, paste('r =', round(cor(d[,2], d[,4], use = 'complete.obs', method = c('pearson')),2)), cex = 2)

plot(d[,3], d[,5], xlab  = 'RNA-seq Log2 FC', ylab = 'NET-seq Log2 FC', main = expression(paste(italic('lsm1'), Delta)),
     col = rgb(0,0,0,.15), pch = 16, cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5, xlim = c(-2,3), ylim = c(-4,4))
text(x = -1.4, y = 3.7, paste('r =', round(cor(d[,3], d[,5], use = 'complete.obs', method = c('pearson')),2)), cex = 2)

plot(NET_GRO_FC$RPB4_NET, NET_GRO_FC$RPB4_GRO, xlab = 'NET-seq Log2 FC', ylab = 'GRO Log2 FC', main = expression(paste(italic('rpb4'), Delta)),
     col = rgb(0,0,0,.15), pch = 16, cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5, xlim = c(-5,5), ylim = c(-4,3))
text(x = -3.85, y = 2.7, paste('r =', round(cor(NET_GRO_FC$RPB4_NET, NET_GRO_FC$RPB4_GRO, use = 'complete.obs', method = c('pearson')),2)), cex = 2)

plot(NET_GRO_FC$XRN1_NET, NET_GRO_FC$XRN1_GRO, xlab = 'NET-seq Log2 FC', ylab = 'GRO Log2 FC', main = expression(paste(italic('xrn1'), Delta)),
     col = rgb(0,0,0,.15), pch = 16, cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5, xlim = c(-5,5), ylim = c(-5,2))
text(x = -3.75, y = 1.7, paste('r =', round(cor(NET_GRO_FC$XRN1_NET, NET_GRO_FC$XRN1_GRO, use = 'complete.obs', method = c('pearson')),2)), cex = 2)
par(mfrow = c(1,1))
dev.off()

##### Figure S3 ####
par(mfrow = c(1,2))
m = merge(merged_mat_2, ORF_annot)
gl2 = log2(abs(m$start-m$end))
y = log2(m$`RPB4_GRO$VALUE`)-log2(m$`WT2_GRO$VALUE`)-log2(m$RPB4+1)+log2(m$WT+1)
plot(gl2, y, xlab = 'Log2 gene length', ylab = 'Log2 GRO/NET-seq FC', main = expression(paste(italic('rpb4'), Delta)),
     cex.main = 2, cex.axis = 2, cex.lab = 2, ylim = c(-6.25, 4.25))
n = intersect(which(!is.na(y)), which(is.finite(y)))
trend = lowess(gl2[n], y[n], f = .1)
lines(trend, col = 2)
lm(y[n] ~ gl2[n])

m = merge(merged_mat_2, ORF_annot)
gl2 = log2(abs(m$start-m$end))
y = log2(m$`XRN1_GRO$VALUE`)-log2(m$`WT_GRO$VALUE`)-log2(m$XRN1+1)+log2(m$WT+1)
plot(gl2, y, xlab = 'Log2 gene length', ylab = 'Log2 GRO/NET-seq FC', main = expression(paste(italic('xrn1'), Delta)),
     cex.main = 2, cex.axis = 2, cex.lab = 2, ylim = c(-6.25, 4.25))
n = intersect(which(!is.na(y)), which(is.finite(y)))
trend = lowess(gl2[n], y[n], f = .1)
lines(trend, col = 2)
lm(y[n] ~ gl2[n])

#### Figure S8 ####
a1 = density(CCR4_ratio - WT_ratio, na.rm = T)
a2 = density(DHH1_ratio - WT_ratio, na.rm = T)
a3 = density(LSM1_ratio - WT_ratio, na.rm = T)
a4 = density(RPB4_ratio - WT_ratio, na.rm = T)
a5 = density(XRN1_ratio - WT_ratio, na.rm = T)

plot(a1, xlim = c(-5,5), xlab = 'Log2 Pol II FC', main = "5'/3' ratio fold change distributions",
     col = 4, lwd = 2, ylim = c(0,.5), cex.main = 2, cex.axis = 2, cex.lab = 2)
lines(a2,col = 2, lwd = 2)
lines(a3,col = 3, lwd = 2)
lines(a4,col = 6, lwd = 2)
lines(a5,col = 1, lwd = 2)
abline(v = 0, col = 1, lty = 2)
legend(x = 'topright', legend = c(expression(paste(italic("xrn1"), Delta)), expression(paste(italic("dhh1"), Delta)), expression(paste(italic("lsm1"), Delta)), expression(paste(italic("ccr4"), Delta)), expression(paste(italic("rpb4"), Delta))),
       col = c(1:4,6), lwd = 2, cex = 2)

##### Figure S14 #####
pdf(file = '~/Documents/R_Projects/NETSEQ/paper/DF2/Figures/review/ChIP_exo.pdf', height = 12.75)
par(mfrow = c(3,1))
par(mar = c(5,5,5,3))
plot(ChIP_prop$Bin, ChIP_prop$XRN1_SAGA, typ = 'l', col = 4, lwd = 2,
     xlab = 'Distance from TSS', ylab = 'Fraction of bound genes',
     main = 'Xrn1 ChIP-exo binding', cex.main = 2, cex.lab = 2, cex.axis = 2, ylim = c(0, .2))
lines(ChIP_prop$Bin, ChIP_prop$XRN1_TFIID, col = 2, lwd = 2)
legend(x = 'topright', legend = c('SAGA-dominated', 'TFIID-dominated'),
       lwd = 2, col = c(4,2), cex = 1.5)

plot(ChIP_prop$Bin, ChIP_prop$LSM1_SAGA, typ = 'l', col = 4, lwd = 2,
     xlab = 'Distance from TSS', ylab = 'Fraction of bound genes',
     main = 'Lsm1 ChIP-exo binding', cex.main = 2, cex.lab = 2, cex.axis = 2, ylim = c(0, .2))
lines(ChIP_prop$Bin, ChIP_prop$LSM1_TFIID, col = 2, lwd = 2)
legend(x = 'topright', legend = c('SAGA-dominated', 'TFIID-dominated'),
       lwd = 2, col = c(4,2), cex = 1.5)

plot(ChIP_prop$Bin, ChIP_prop$DCP2_SAGA, typ = 'l', col = 4, lwd = 2,
     xlab = 'Distance from TSS', ylab = 'Fraction of bound genes',
     main = 'Dcp2 ChIP-exo binding', cex.main = 2, cex.lab = 2, cex.axis = 2, ylim = c(0, .2))
lines(ChIP_prop$Bin, ChIP_prop$DCP2_TFIID, col = 2, lwd = 2)
legend(x = 'topright', legend = c('SAGA-dominated', 'TFIID-dominated'),
       lwd = 2, col = c(4,2), cex = 1.5)
par(mfrow = c(1,1))
dev.off()

##### Figure S16 #####
pdf(file = 'ChIP_exo_ncRNA.pdf', height = 8)
par(mfrow = c(3,1))
par(mar = c(5,5,5,3))
plot(seq(-270, 270, 60), colMeans(DCP2_ChIP_mat  >= 10 ), typ = 'l', xlab = 'Position relative to TSS', ylab = 'Fraction of bound loci', main = 'Dcp2',
     cex.main = 2, cex.lab = 1.5, ylim = c(0,.13), cex.axis = 1.5)
lines(seq(-270, 270, 60), colMeans(DCP2_ChIP_mat_ncRNA >= 10), typ = 'l', col = 2)
lines(seq(-270, 270, 60), colMeans(control_ChIP_mat >= 10), lty = 2)
lines(seq(-270, 270, 60), colMeans(control_ChIP_mat_ncRNA >= 10), lty = 2, col = 2)
legend(x = 'topleft', legend = c('Coding genes', 'ncRNAs'), fill = c(1,2))
legend(x = 'topright', legend = c('Dcp2', 'control'), lty = c(1,2))
abline(v = 0, lty = 3, col = 1)

plot(seq(-270, 270, 60), colMeans(LSM1_ChIP_mat  >= 10 ), typ = 'l', xlab = 'Position relative to TSS', ylab = 'Fraction of bound loci', main = 'Lsm1',
     cex.main = 2, cex.lab = 1.5, ylim = c(0,.13), cex.axis = 1.5)
lines(seq(-270, 270, 60), colMeans(LSM1_ChIP_mat_ncRNA >= 10), typ = 'l', col = 2)
lines(seq(-270, 270, 60), colMeans(control_ChIP_mat >= 10), lty = 2)
lines(seq(-270, 270, 60), colMeans(control_ChIP_mat_ncRNA >= 10), lty = 2, col = 2)
legend(x = 'topleft', legend = c('Coding genes', 'ncRNAs'), fill = c(1,2))
legend(x = 'topright', legend = c('Lsm1', 'control'), lty = c(1,2))
abline(v = 0, lty = 3, col = 1)

plot(seq(-270, 270, 60), colMeans(XRN1_ChIP_mat  >= 10 ), typ = 'l', xlab = 'Position relative to TSS', ylab = 'Fraction of bound loci', main = 'Xrn1',
     cex.main = 2, cex.lab = 1.5, ylim = c(0,.13), cex.axis = 1.5)
lines(seq(-270, 270, 60), colMeans(XRN1_ChIP_mat_ncRNA >= 10), typ = 'l', col = 2)
lines(seq(-270, 270, 60), colMeans(control_ChIP_mat >= 10), lty = 2)
lines(seq(-270, 270, 60), colMeans(control_ChIP_mat_ncRNA >= 10), lty = 2, col = 2)
legend(x = 'topleft', legend = c('Coding genes', 'ncRNAs'), fill = c(1,2))
legend(x = 'topright', legend = c('Xrn1', 'control'), lty = c(1,2))
abline(v = 0, lty = 3, col = 1)
dev.off()

#############
pdf(file = 'TSS_metagenes_RiBi.pdf', width = 14)
end_metagene_plotter(gene_ind = which(ORF_annot$name %in% rb_genes$Gene.Systematic.Name), gene_title = 'RiBi genes', which_end = 'TSS')
dev.off()
pdf(file = 'PAS_metagenes_RiBi.pdf', width = 14)
end_metagene_plotter(gene_ind = which(ORF_annot$name %in% rb_genes$Gene.Systematic.Name), gene_title = 'RiBi genes', which_end = 'PAS')
dev.off()

pdf(file = 'TSS_metagenes_RP.pdf', width = 14)
end_metagene_plotter(gene_ind = which(ORF_annot$name %in% rp_genes$Gene.Systematic.Name), gene_title = 'RP genes', which_end = 'TSS')
dev.off()
pdf(file = 'PAS_metagenes_RP.pdf', width = 14)
end_metagene_plotter(gene_ind = which(ORF_annot$name %in% rp_genes$Gene.Systematic.Name), gene_title = 'RP genes', which_end = 'PAS')
dev.off()

pdf(file = 'TSS_metagenes_xrn1_synthegradon.pdf', width = 14)
end_metagene_plotter(gene_ind = which(ORF_annot$name %in% xrn1_synth_ratings$ORF.name[xrn1_synth_ratings$responsiveness >= 8]), gene_title = 'Xrn1 synthegradon', which_end = 'TSS')
dev.off()
pdf(file = 'PAS_metagenes_xrn1_synthegradon.pdf', width = 14)
end_metagene_plotter(gene_ind = which(ORF_annot$name %in% xrn1_synth_ratings$ORF.name[xrn1_synth_ratings$responsiveness >= 8]), gene_title = 'Xrn1 sythegradon', which_end = 'PAS')
dev.off()

pdf(file = 'TSS_metagenes_xrn1_antisynthegradon.pdf', width = 14)
end_metagene_plotter(gene_ind = which(ORF_annot$name %in% xrn1_synth_ratings$ORF.name[xrn1_synth_ratings$responsiveness <= 4]), gene_title = 'Xrn1 antisynthegradon', which_end = 'TSS')
dev.off()
pdf(file = 'PAS_metagenes_xrn1_antisynthegradon.pdf', width = 14)
end_metagene_plotter(gene_ind = which(ORF_annot$name %in% xrn1_synth_ratings$ORF.name[xrn1_synth_ratings$responsiveness <= 4]), gene_title = 'Xrn1 antisythegradon', which_end = 'PAS')
dev.off()

pdf(file = 'TSS_metagenes_SAGA.pdf', width = 14)
end_metagene_plotter(gene_ind = which(ORF_annot$name %in% SAGA_dominated$Name), gene_title = 'SAGA-dominated genes', which_end = 'TSS')
dev.off()
pdf(file = 'PAS_metagenes_SAGA.pdf', width = 14)
end_metagene_plotter(gene_ind = which(ORF_annot$name %in% SAGA_dominated$Name), gene_title = 'SAGA-dominated genes', which_end = 'PAS')
dev.off()

pdf(file = 'TSS_metagenes_TFIID.pdf', width = 14)
end_metagene_plotter(gene_ind = which(ORF_annot$name %in% TFIID_dominated$Name), gene_title = 'TFIID-dominated genes', which_end = 'TSS')
dev.off()
pdf(file = 'PAS_metagenes_TFIID.pdf', width = 14)
end_metagene_plotter(gene_ind = which(ORF_annot$name %in% TFIID_dominated$Name), gene_title = 'TFIID-dominated genes', which_end = 'PAS')
dev.off()

##########################################################################################
a = WT_TSS_prof_NUT + WT_pilot_TSS_prof_NUT
a[a > 250] = 0
plot(lowess(rowMeans(a)[401:1001], f = .05)$y/(nf[1]+nf[2]), typ = 'l')
a = XRN1_TSS_prof_NUT + XRN1_TSS_prof_NUT
a[a > 250] = 0
plot(lowess(rowMeans(a)[401:1001], f = .05)$y/(nf[1]+nf[2]), typ = 'l')
plot(lowess(rowMeans(a)[401:1001], f = .05)$y/(nf[1]+nf[2]), typ = 'l')
plot(lowess(rowMeans(a)[401:1001], f = .05)$y/(nf[1]+nf[2]), typ = 'l')
plot(lowess(rowMeans(a)[401:1001], f = .05)$y/(nf[1]+nf[2]), typ = 'l')
plot(lowess(rowMeans(a)[401:1001], f = .05)$y/(nf[1]+nf[2]), typ = 'l')

rs = sample(ncol(LSM1_TSS), 300)
plot(rowMeans(LSM1_TSS[,rs]), typ = 'l', lwd = 2)
lines(lowess(rowMeans(LSM1_TSS[,rs]), f = .01), col = 2, lwd = 2)
lines(lowess(rowMeans(LSM1_TSS[,rs]), f = .05), col = 3, lwd = 2)
lines(lowess(rowMeans(LSM1_TSS[,rs]), f = .15), col = 4, lwd = 2)
legend(x = 'topleft', col = 1:4, lwd = 2, legend = c(0.0, 0.01, 0.05, 0.15),
       title = 'Amount of smoothing')
##########################################################################################
# Fig S17

conv_col = function(pl = WT_pilot_plus, mi = WT_pilot_minus, nf = norm_factors[2], cap = 200){
  
  a = lapply(ind_list, function(x) lapply(x, function(y) interval_reads(start = pairs$Int1[y], end = pairs$Int2[y], strand = '+', plus_list = pl, minus_list = mi,
                                                                        chrom = pairs$chr[y], sums = F, plot_flag = F)))
  l1 = lapply(a, function(y) sapply(y, function(x) x[,1]))
  l2 = lapply(a, function(y) sapply(y, function(x) x[,2]))
  
  if(!is.null(cap)){
    l1 = lapply(l1, function(x){x[x > cap*nf] = 0; return(x)})
    l2 = lapply(l2, function(x){x[x > cap*nf] = 0; return(x)})
  }
  
  p1 = sapply(l1, rowMeans)
  p2 = sapply(l2, rowMeans)
  
  #par(mfrow = c(2,3))
  
  #r1 = lowess(p1[,1], f = .1)$y/nf
  #r2 = lowess(p2[,1], f = .1)$y/nf
  #plot(r1, typ = 'l', ylim = c(-max(r2),max(r1)), xlab = 'Dist from midpoint', ylab = 'Normalized metaprofile', main = 'd < -200')
  #lines(-r2, col = 'red')
  #abline(v = c(500, 600, 400), col = c('blue', 'black', 'red'), lty = c(1,2,2))
  
  #r1 = lowess(p1[,2], f = .1)$y/nf
  #r2 = lowess(p2[,2], f = .1)$y/nf
  #plot(r1, typ = 'l', ylim = c(-max(r2),max(r1)),
  #     xlab = 'Dist from midpoint', ylab = 'Normalized metaprofile', main = '-200 < d < -100')
  #lines(-r2, col = 'red')
  #abline(v = c(500, 575, 425), col = c('blue', 'black', 'red'), lty = c(1,2,2))
  
  #r1 = lowess(p1[,3], f = .1)$y/nf
  #r2 = lowess(p2[,3], f = .1)$y/nf
  #plot(r1, typ = 'l', ylim = c(-max(r2),max(r1)),
  #     xlab = 'Dist from midpoint', ylab = 'Normalized metaprofile', main = '-100 < d < 0')
  #lines(-r2, col = 'red')
  #abline(v = c(500, 525, 475), col = c('blue', 'black', 'red'), lty = c(1,2,2))
  
  #r1 = lowess(p1[,4], f = .1)$y/nf
  #r2 = lowess(p2[,4], f = .1)$y/nf
  #plot(r1, typ = 'l', ylim = c(-max(r2),max(r1)),
  #     xlab = 'Dist from midpoint', ylab = 'Normalized metaprofile', main = '0 < d < 100')
  #lines(-r2, col = 'red')
  #abline(v = c(500, 475, 525), col = c('blue', 'black', 'red'), lty = c(1,2,2))
  
  #r1 = lowess(p1[,5], f = .1)$y/nf
  #r2 = lowess(p2[,5], f = .1)$y/nf
  #plot(r1, typ = 'l', ylim = c(-max(r2),max(r1)),
  #     xlab = 'Dist from midpoint', ylab = 'Normalized metaprofile', main = '100 < d < 200')
  #lines(-r2, col = 'red')
  #abline(v = c(500, 425, 575), col = c('blue', 'black', 'red'), lty = c(1,2,2))
  
  #r1 = lowess(p1[,6], f = .1)$y/nf
  #r2 = lowess(p2[,6], f = .1)$y/nf
  #plot(r1, typ = 'l', ylim = c(-max(r2),max(r1)),
  #     xlab = 'Dist from midpoint', ylab = 'Normalized metaprofile', main = '200 < d < 300')
  #lines(-r2, col = 'red')
  #abline(v = c(500, 375, 625), col = c('blue', 'black', 'red'), lty = c(1,2,2))
  
  to_return = list(l1, l2, p1, p2)
  names(to_return) = c('plus_mat_list', 'minus_mat_list', 'plus_avg', 'minus_avg')
  return(to_return)
}

conv_sub = convergent_pairs[which(convergent_pairs$gene1 %in% ORF_annot$name & convergent_pairs$gene2 %in% ORF_annot$name),]
a = cbind(ORF_annot[match(conv_sub$gene1, ORF_annot$name),], ORF_annot[match(conv_sub$gene2, ORF_annot$name),])
a = a[,-c(2,7,8)]
a$midpoint = as.integer((a$end + a$start.1)/2)
a$dist = a$end - a$start.1

pairs = a[,c(1,10,11)]
pairs$Int1 = pairs$midpoint - 500
pairs$Int2 = pairs$midpoint + 500
ind_list = 1:nrow(pairs)

conv_col_WT = conv_col(pl = WT_plus, mi = WT_minus, nf = norm_factors[1])
conv_col_WT_pilot = conv_col(pl = WT_pilot_plus, mi = WT_pilot_minus, nf = norm_factors[2])
conv_col_ccr4 = conv_col(pl = CCR4_plus, mi = CCR4_minus, nf = norm_factors[3])
conv_col_dhh1 = conv_col(pl = DHH1_plus, mi = DHH1_minus, nf = norm_factors[4])
conv_col_lsm1 = conv_col(pl = LSM1_plus, mi = LSM1_minus, nf = norm_factors[5])
conv_col_rpb4 = conv_col(pl = RPB4_plus, mi = RPB4_minus, nf = norm_factors[6])
conv_col_sfp1 = conv_col(pl = SFP1_plus, mi = SFP1_minus, nf = norm_factors[7])
conv_col_xrn1 = conv_col(pl = XRN1_plus, mi = XRN1_minus, nf = norm_factors[8])
conv_col_xrn1_pilot = conv_col(pl = XRN1_pilot_plus, mi = XRN1_pilot_minus, nf = norm_factors[9])

wt_conv = (do.call(cbind, conv_col_WT$plus_mat_list)+do.call(cbind, conv_col_WT$minus_mat_list)+
             do.call(cbind, conv_col_WT_pilot$plus_mat_list)+do.call(cbind, conv_col_WT_pilot$minus_mat_list))/(norm_factors[1]+norm_factors[2])
ccr4_conv = do.call(cbind, conv_col_ccr4$plus_mat_list)/norm_factors[3]+do.call(cbind, conv_col_ccr4$minus_mat_list)/norm_factors[3]
dhh1_conv = do.call(cbind, conv_col_dhh1$plus_mat_list)/norm_factors[4]+do.call(cbind, conv_col_dhh1$minus_mat_list)/norm_factors[4]
lsm1_conv = do.call(cbind, conv_col_lsm1$plus_mat_list)/norm_factors[5]+do.call(cbind, conv_col_lsm1$minus_mat_list)/norm_factors[5]
rpb4_conv = do.call(cbind, conv_col_rpb4$plus_mat_list)/norm_factors[6]+do.call(cbind, conv_col_rpb4$minus_mat_list)/norm_factors[6]
sfp1_conv = do.call(cbind, conv_col_sfp1$plus_mat_list)/norm_factors[7]+do.call(cbind, conv_col_sfp1$minus_mat_list)/norm_factors[7]
xrn1_conv = (do.call(cbind, conv_col_xrn1$plus_mat_list)+do.call(cbind, conv_col_xrn1$minus_mat_list)+
               do.call(cbind, conv_col_xrn1_pilot$plus_mat_list)+do.call(cbind, conv_col_xrn1_pilot$minus_mat_list))/(norm_factors[8]+norm_factors[9])

a = wt_conv
inds = list(which(pairs$dist > 200), which(pairs$dist > 100 & pairs$dist < 200),
            which(pairs$dist > 0 & pairs$dist < 100), which(pairs$dist < 0 & pairs$dist > -100),
            which(pairs$dist > -200 & pairs$dist < -100), which(pairs$dist < -200 & pairs$dist > -400))


for(i in 1:5){
  b = list(ccr4_conv, dhh1_conv, lsm1_conv, rpb4_conv, xrn1_conv)[[i]]
  nm = c('ccr4', 'dhh1', 'lsm1', 'rpb4', 'xrn1')[i]
  
  mut = c(expression(paste(italic('ccr4'), Delta)),
          expression(paste(italic('dhh1'), Delta)),
          expression(paste(italic('lsm1'), Delta)),
          expression(paste(italic('rpb4'), Delta)),
          expression(paste(italic('xrn1'), Delta)))[i]
  
  pdf(paste0('conv_by_dist_', nm, '.pdf'), width = 10)
  
  par(mfrow = c(3,2))
  par(mar = c(5,5,5,3))
  x = lowess(rowMeans(a[,inds[[6]]]), f = .05)$y
  y = lowess(rowMeans(b[,inds[[6]]]), f = .05)$y
  plot(-500:500, x, typ = 'l', main = 'Gene distance between -200 and -400', ylab = 'Normalized mean reads',
       xlab = 'Distance from midpoint', ylim = c(0.2, 0.9), cex.lab = 1.7, cex.axis = 2, cex.main = 2)
  lines(-500:500, y, col = 'red')
  legend(x = 'topleft', cex = 1, col = c(1,2), lwd = 1, legend = c('WT', mut))
  
  x = lowess(rowMeans(a[,inds[[5]]]), f = .05)$y
  y = lowess(rowMeans(b[,inds[[5]]]), f = .05)$y
  plot(-500:500, x, typ = 'l', main = 'Gene distance between -200 & -100', ylab = 'Normalized mean reads',
       xlab = 'Distance from midpoint', ylim = c(.15,.4), cex.lab = 1.7, cex.axis = 2, cex.main = 2)
  lines(-500:500, y, col = 'red')
  legend(x = 'topright', cex = 1, col = c(1,2), lwd = 1, legend = c('WT', mut))
  
  x = lowess(rowMeans(a[,inds[[4]]]), f = .05)$y
  y = lowess(rowMeans(b[,inds[[4]]]), f = .05)$y
  plot(-500:500, x, typ = 'l', main = 'Gene distance between -100 and 0', ylab = 'Normalized mean reads',
       xlab = 'Distance from midpoint', ylim = c(.15, .5), cex.lab = 1.7, cex.axis = 2, cex.main = 2)
  lines(-500:500, y, col = 'red')
  legend(x = 'topright', cex = 1, col = c(1,2), lwd = 1, legend = c('WT', mut))
  
  x = lowess(rowMeans(a[,inds[[3]]]), f = .05)$y
  y = lowess(rowMeans(b[,inds[[3]]]), f = .05)$y
  plot(-500:500, x, typ = 'l', main = 'Gene distance between 0 and 100', ylab = 'Normalized mean reads',
       xlab = 'Distance from midpoint', ylim = c(.12, .3), cex.lab = 1.7, cex.axis = 2, cex.main = 2)
  lines(-500:500, y, col = 'red')
  legend(x = 'topright', cex = 1, col = c(1,2), lwd = 1, legend = c('WT', mut))
  
  x = lowess(rowMeans(a[,inds[[2]]]), f = .05)$y
  y = lowess(rowMeans(b[,inds[[2]]]), f = .05)$y
  plot(-500:500, x, typ = 'l', main = 'Gene distance between 100 and 200', ylab = 'Normalized mean reads',
       xlab = 'Distance from midpoint', ylim = c(.08, .2), cex.lab = 1.7, cex.axis = 2, cex.main = 2)
  lines(-500:500, y, col = 'red')
  legend(x = 'topright', cex = 1, col = c(1,2), lwd = 1, legend = c('WT', mut))
  
  x = lowess(rowMeans(a[,inds[[1]]]), f = .05)$y
  y = lowess(rowMeans(b[,inds[[1]]]), f = .05)$y
  plot(-500:500, x, typ = 'l', main = 'Gene distance between 200 and 400', ylab = 'Normalized mean reads',
       xlab = 'Distance from midpoint', ylim = c(.05, .2), cex.lab = 1.7, cex.axis = 2, cex.main = 2)
  lines(-500:500, y, col = 'red')
  legend(x = 'topright', cex = 1, col = c(1,2), lwd = 1, legend = c('WT', mut))
  dev.off()
}




