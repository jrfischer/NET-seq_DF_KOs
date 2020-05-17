big_annot = rbind(ORF_annot[,-c(5,6)], ncRNA_annot)
site_coverage_list = coverage_list(annot = big_annot)
full_neighbor_mat = data.frame(t(sapply(1:nrow(ORF_annot), function(x) gene_neighbors(start = ORF_annot$start[x], end = ORF_annot$end[x], chrom = ORF_annot$chr[x], strand = ORF_annot$strand[x], radius = 500, cov_list = site_coverage_list))))
rownames(full_neighbor_mat) = ORF_annot$name

CUT_overlaps = data.frame(overlap_full(my_annot = CUT_annot, ref_annot = ORF_annot))
NUT_overlaps = data.frame(overlap_full(my_annot = NUT_annot, ref_annot = ORF_annot))
SUT_overlaps = data.frame(overlap_full(my_annot = SUT_annot, ref_annot = ORF_annot))
XUT_overlaps = data.frame(overlap_full(my_annot = XUT_annot, ref_annot = ORF_annot))
ORF_overlaps = data.frame(overlap_full(my_annot = ORF_annot, ref_annot = ORF_annot))

ncRNA_overlaps = unique(c(CUT_overlaps$Overlap, NUT_overlaps$Overlap, SUT_overlaps$Overlap, XUT_overlaps$Overlap))
z = c(unique(CUT_overlaps$Overlap), unique(NUT_overlaps$Overlap), unique(SUT_overlaps$Overlap), unique(XUT_overlaps$Overlap), unique(ORF_overlaps$Overlap))

#table(1:4973 %in% z, full_neighbor_mat$AS_overlap)

nf = norm_factors/exp(mean(log(norm_factors)))
strains = c('WT', 'WT pilot', 'Ccr4', 'Dhh1', 'Lsm1', 'Rpb4', 'Sfp1', 'Xrn1', 'Xrn1 pilot')
#pdf('convergent_transcripts.pdf', width = 12)
par(mfrow = c(1,3))
for(i in 1:9){
  boxplot(log2(ORF_table[unique(CUT_overlaps$Overlap),i]/nf[i]+1), 
          log2(ORF_table[unique(NUT_overlaps$Overlap),i]/nf[i]+1), 
          log2(ORF_table[unique(SUT_overlaps$Overlap),i]/nf[i]+1),
          log2(ORF_table[unique(XUT_overlaps$Overlap),i]/nf[i]+1),
          log2(ORF_table[unique(ORF_overlaps$Overlap),i]/nf[i]+1),
          log2(ORF_table[setdiff(1:4973,z),i]/nf[i]+1),
          names = c('CUT', 'NUT', 'SUT', 'XUT', 'ORF', 'None'),
          main = strains[i], xlab = 'Convergent transcript type', ylab = 'Normalized log2 reads in genes')
}
#dev.off()

convergent_cor = function(overlap_mat = XUT_overlaps, s_reads = ORF_table, as_reads = XUT_reads, plot_flag = F){
  
  cor_vec = sapply(1:ncol(s_reads), function(x) cor(s_reads[overlap_mat$Overlap,x], as_reads[overlap_mat$Ind,x], method = 'spearman', use = 'complete.obs'))
  
  return(cor_vec)
  
}

convergent_cor_mat = cbind(convergent_cor(overlap_mat = CUT_overlaps, s_reads = our_sense_reads[,-c(1:6)], as_reads = CUT_reads, plot_flag = F),
                           convergent_cor(overlap_mat = NUT_overlaps, s_reads = our_sense_reads[,-c(1:6)], as_reads = NUT_reads, plot_flag = F),
                           convergent_cor(overlap_mat = SUT_overlaps, s_reads = our_sense_reads[,-c(1:6)], as_reads = SUT_reads, plot_flag = F),
                           convergent_cor(overlap_mat = XUT_overlaps, s_reads = our_sense_reads[,-c(1:6)], as_reads = XUT_reads, plot_flag = F),
                           convergent_cor(overlap_mat = ORF_overlaps, s_reads = our_sense_reads[,-c(1:6)], as_reads = our_sense_reads[,-c(1:6)], plot_flag = F))

rownames(convergent_cor_mat) = colnames(ORF_table)
colnames(convergent_cor_mat) = c('CUT', 'NUT', 'SUT', 'XUT', 'ORF')

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
ccm = convergent_cor_mat[-c(2,7,9),]
ccm[1,] = .5*(convergent_cor_mat[1,]+convergent_cor_mat[2,])
ccm[6,] = .5*(convergent_cor_mat[8,]+convergent_cor_mat[9,])
rownames(ccm) = c('WT', ":italic('ccr4') * Delta", ":italic('dhh1') * Delta", ":italic('lsm1') * Delta", ":italic('rpb4') * Delta", ":italic('xrn1') * Delta")
corrplot(ccm, method="color", col=col(200), addCoef.col = "black", tl.col="black", tl.srt=0, mar = c(1,1,2,1))

convergent_fc_cor_mat = cbind(convergent_cor(overlap_mat = CUT_overlaps, s_reads = ORF_fc, as_reads = CUT_fc, plot_flag = F),
                              convergent_cor(overlap_mat = NUT_overlaps, s_reads = ORF_fc, as_reads = NUT_fc, plot_flag = F),
                              convergent_cor(overlap_mat = SUT_overlaps, s_reads = ORF_fc, as_reads = SUT_fc, plot_flag = F),
                              convergent_cor(overlap_mat = XUT_overlaps, s_reads = ORF_fc, as_reads = XUT_fc, plot_flag = F),
                              convergent_cor(overlap_mat = ORF_overlaps, s_reads = ORF_fc, as_reads = ORF_fc, plot_flag = F))
rownames(convergent_fc_cor_mat) = colnames(ORF_fc)
colnames(convergent_fc_cor_mat) = c('CUT', 'NUT', 'SUT', 'XUT', 'ORF')

rownames(convergent_fc_cor_mat) = c(":italic('ccr4') * Delta", ":italic('dhh1') * Delta", ":italic('lsm1') * Delta", ":italic('rpb4') * Delta", ":italic('xrn1') * Delta")
corrplot(convergent_fc_cor_mat, method="color", col=col(200), addCoef.col = "black", tl.col="black", tl.srt=0)


no_overlap = setdiff(1:4973, c(unique(CUT_overlaps$Overlap), unique(NUT_overlaps$Overlap), unique(SUT_overlaps$Overlap), unique(XUT_overlaps$Overlap), unique(ORF_overlaps$Overlap)))

#par(mfrow = c(2,3))
par(mar = c(5,5,5,3))
for(i in 1:5){
  a = sapply(list(ORF_fc[unique(CUT_overlaps$Overlap),i], ORF_fc[unique(NUT_overlaps$Overlap),i],
                  ORF_fc[unique(SUT_overlaps$Overlap),i], ORF_fc[unique(XUT_overlaps$Overlap),i],
                  ORF_fc[unique(ORF_overlaps$Overlap),i], ORF_fc[unique(no_overlap),i]), function(y) median(y, na.rm = T))
  
  boxplot(ORF_fc[unique(CUT_overlaps$Overlap),i], ORF_fc[unique(NUT_overlaps$Overlap),i],
          ORF_fc[unique(SUT_overlaps$Overlap),i], ORF_fc[unique(XUT_overlaps$Overlap),i],
          ORF_fc[unique(ORF_overlaps$Overlap),i], ORF_fc[unique(no_overlap),i],
          xlab = 'Overlapping transcript type', names = c('CUT', 'NUT', 'SUT', 'XUT', 'ORF', 'None'),
          main = expression(paste(italic('xrn1'),Delta)), ylab = 'Gene log2 FC', ylim = c(-5,5), notch = T, cex.axis = 2, cex.lab = 2, cex.main = 2)
  abline(h = 0, lty = 2, col = 2)
  #text(x = 1:6, y = 5, round(a,2), col = 2, cex = 2)
}


kruskal.test(list(ORF_fc[unique(CUT_overlaps$Overlap),i], ORF_fc[unique(NUT_overlaps$Overlap),i],
                  ORF_fc[unique(SUT_overlaps$Overlap),i], ORF_fc[unique(XUT_overlaps$Overlap),i],
                  ORF_fc[unique(ORF_overlaps$Overlap),i], ORF_fc[unique(no_overlap),i]))

par(mfrow = c(3,3))
for(i in 1:9){
  a = cor(log2(ORF_table[,i]+1), log2(our_antisense_reads[,i+5]+1), method = 'spearman')
  plot(log2(ORF_table[,i]+1), log2(our_antisense_reads[,i+5]+1), col = (full_neighbor_mat$AS_overlap == 1)+1, xlab = 'Gene reads', ylab = 'Overlapping antisense reads', main = paste(colnames(ORF_table)[i], round(a,2)))
  legend(x = 'topleft', legend = c('No annotated AS transcript', 'Annotated AS transcript'), pch = 1, col = 1:2, cex = .65)
}

j = 1
plot(ORF_fc[CUT_overlaps$Overlap,j], CUT_fc[CUT_overlaps$Ind,j])
points(ORF_fc[NUT_overlaps$Overlap,j], NUT_fc[NUT_overlaps$Ind,j])
points(ORF_fc[SUT_overlaps$Overlap,j], SUT_fc[SUT_overlaps$Ind,j])
points(ORF_fc[XUT_overlaps$Overlap,j], XUT_fc[XUT_overlaps$Ind,j])
points(ORF_fc[ORF_overlaps$Overlap,j], ORF_fc[ORF_overlaps$Ind,j], col = 2)

plot(rowMeans(WT_pilot_TES[,-c(2447, 4133, 4898, 1431)]), typ = 'l')

plot(1, type="n", xlab="Distance from TSS", ylab="Average reads", xlim=c(0, 500), ylim=c(0, 4))
for(i in 1:5){
  r = which(ORF_annot$name %in% xrn1_synth_ratings$ORF.name[xrn1_synth_ratings$TR.change.ranking == i])
  lines(rowMeans(DHH1_TSS[,r]), col = i)
}

plot(1, type="n", xlab="Distance from TES", ylab="Average reads", xlim=c(-150, 150), ylim=c(0, 3))
for(i in 1:5){
  r = which(ORF_annot$name %in% xrn1_synth_ratings$ORF.name[xrn1_synth_ratings$TR.change.ranking == i])
  lines(-150:150, rowMeans(WT_pilot_TES[,r]), col = i)
}




a = lapply(1:5, function(x) xrn1_synth_ratings$ORF.name[xrn1_synth_ratings$TR.change.ranking == x])
a = lapply(a, function(x) x[x %in% ORF_annot$name])
b = lapply(a, function(x) match(x, ORF_annot$name))

overlap_props_synth = rbind(sapply(b, function(x) mean(x %in% CUT_overlaps$Overlap)),
                            sapply(b, function(x) mean(x %in% NUT_overlaps$Overlap)),
                            sapply(b, function(x) mean(x %in% SUT_overlaps$Overlap)),
                            sapply(b, function(x) mean(x %in% XUT_overlaps$Overlap)),
                            sapply(b, function(x) mean(x %in% ncRNA_overlaps)),
                            sapply(b, function(x) mean(x %in% ORF_overlaps$Overlap)),
                            sapply(b, function(x) mean(x %in% no_overlap)))

colnames(overlap_props_synth) = c('TR1', 'TR2', 'TR3', 'TR4', 'TR5')
rownames(overlap_props_synth) = c('CUT', 'NUT', 'SUT', 'XUT', 'ncRNA', 'ORF', 'None')

plot(overlap_props_synth[1,], typ = 'l', xlim = c(1,5), ylim = c(0,.6), lwd = 2,
     xlab = 'Xrn1 TR sensitivity', ylab = 'Proportion of genes with convergent AS transcript',
     main = 'Convergent transcript type and rates')
lines(overlap_props_synth[2,], col = 2, lwd = 2)
lines(overlap_props_synth[3,], col = 3, lwd = 2)
lines(overlap_props_synth[4,], col = 4, lwd = 2)
lines(overlap_props_synth[5,], col = 5, lwd = 2)
lines(overlap_props_synth[6,], col = 6, lwd = 2)
lines(overlap_props_synth[7,], col = 8, lwd = 2)
legend(x = 'topleft', lwd = 2, col = c(1:6,8), legend = c('CUT', 'NUT', 'SUT', 'XUT', 'ncRNA', 'ORF', 'None'),
       cex = .75)

par(mfrow = c(2,3))
i = 6
fc_list = lapply(1:5, function(x) list(ORF_fc[b[[x]] %in% ncRNA_overlaps,i], ORF_fc[b[[x]] %in% ORF_overlaps$Overlap,i], ORF_fc[b[[x]] %in% no_overlap,i]))
h = sapply(fc_list, function(x) sapply(x, function(e) median(e, na.rm = T)))
rownames(h) = c('ncRNA', 'ORF', 'None')
colnames(h) = 1:5

plot(h[1,], typ = 'l', col = 1, lwd = 2, ylim = c(min(h), max(h)), ylab = 'Median log2 FC', xlab = 'Xrn1 TR sensitivity',
     main = colnames(ORF_fc)[i])
lines(h[2,], col = 2, lwd = 2)
lines(h[3,], col = 3, lwd = 3)
legend(x = 'topleft', col = 1:3, lwd = 2, legend = c('ncRNA', 'ORF', 'None'), cex = .5)
