GRO_ref = read.delim2('~/Documents/R_Projects/NETSEQ/data/GSE29519_RAW/gro_reference.txt', header = F)
colnames(GRO_ref) = c('ID', 'Probe_Length', 'ORF', 'ORF_Name', 'Gene_Name', 'Molecular_Function', 'Biological_Process', 'SPOT_ID')

################################################################################################

WT2_GRO = read.delim('data/WT2_GRO.txt')
RPB4_GRO = read.delim('data/rpb4_GRO.txt')
SFP1_GRO = read.delim('data/sfp1_GRO.txt')
DST1_GRO = read.delim('data/dst1_GRO.txt')

WT_GRO = read.delim('data/WT_GRO.txt')
XRN1_GRO = read.delim('data/xrn1_GRO.txt')

cdta = read.delim('data/cdta.txt')
ccr4_norm_names = cdta$orf[cdta[,45] <= .25 & cdta[,45] >= -.25]
a = match(ccr4_norm_names, ORF_annot$name)
a = a[!is.na(a)]
n1 = estimateSizeFactorsForMatrix(ORF_table[a,c(2,3)])
estimateSizeFactorsForMatrix(ORF_table[,c(2,3)])
ccr4_norm_factor = n1[2]/n1[1]
plot(log2(ORF_table[,c(2,3)]+1), xlab = 'WT pilot', main = 'Raw Log2 NET-seq reads',
     col = rgb(0,0,0,.5), pch = 16, ylab = expression(paste(italic('ccr4'), Delta)))
points(log2(ORF_table[a,c(2,3)]+1), col = rgb(1,0,0,.5), pch = 16)
legend(x = 'bottomright', col = c(1,2), legend = c('Not used', paste0('Used ', '(', length(a), ')')),
       title = 'Normalization', pch = 16)

dhh1_norm_names = cdta$orf[cdta[,105] <= .25 & cdta[,105] >= -.25]
a = match(dhh1_norm_names, ORF_annot$name)
a = a[!is.na(a)]
n1 = estimateSizeFactorsForMatrix(ORF_table[a,c(2,4)])
estimateSizeFactorsForMatrix(ORF_table[,c(2,4)])
dhh1_norm_factor =  n1[2]/n1[1]
plot(log2(ORF_table[,c(2,4)]+1), xlab = 'WT pilot', main = 'Raw Log2 NET-seq reads',
     col = rgb(0,0,0,.5), pch = 16, ylab = expression(paste(italic('dhh1'), Delta)))
points(log2(ORF_table[a,c(2,4)]+1), col = rgb(1,0,0,.5), pch = 16)
legend(x = 'bottomright', col = c(1,2), legend = c('Not used', paste0('Used ', '(', length(a), ')')),
       title = 'Normalization', pch = 16)

lsm1_norm_names = cdta$orf[cdta[,189] <= .25 & cdta[,189] >= -.25]
a = match(lsm1_norm_names, ORF_annot$name)
a = a[!is.na(a)]
n1 = estimateSizeFactorsForMatrix(ORF_table[a,c(2,5)])
estimateSizeFactorsForMatrix(ORF_table[,c(2,5)])
lsm1_norm_factor =  n1[2]/n1[1]
plot(log2(ORF_table[,c(2,5)]+1), xlab = 'WT pilot', main = 'Raw Log2 NET-seq reads',
     col = rgb(0,0,0,.5), pch = 16, ylab = expression(paste(italic('lsm1'), Delta)))
points(log2(ORF_table[a,c(2,5)]+1), col = rgb(1,0,0,.5), pch = 16)
legend(x = 'bottomright', col = c(1,2), legend = c('Not used', paste0('Used ', '(', length(a), ')')),
       title = 'Normalization', pch = 16)

rpb4_norm_names = GRO_ref$ORF_Name[abs(log2(RPB4_GRO$VALUE+1)-log2(WT2_GRO$VALUE+1)) <= .25]
a = match(rpb4_norm_names, ORF_annot$name)
a = a[!is.na(a)]
n1 = estimateSizeFactorsForMatrix(ORF_table[a,c(2,6)])
estimateSizeFactorsForMatrix(ORF_table[,c(2,6)])
rpb4_norm_factor = n1[2]/n1[1]
plot(log2(ORF_table[,c(2,6)]+1), xlab = 'WT pilot', main = 'Raw Log2 NET-seq reads',
     col = rgb(0,0,0,.5), pch = 16, ylab = expression(paste(italic('rpb4'), Delta)))
points(log2(ORF_table[a,c(2,6)]+1), col = rgb(1,0,0,.5), pch = 16)
legend(x = 'bottomright', col = c(1,2), legend = c('Not used', paste0('Used ', '(', length(a), ')')),
       title = 'Normalization', pch = 16)

sfp1_norm_names = GRO_ref$ORF_Name[abs(log2(SFP1_GRO$VALUE+1)-log2(WT2_GRO$VALUE+1)) <= .25]
a = match(sfp1_norm_names, ORF_annot$name)
a = a[!is.na(a)]
n1 = estimateSizeFactorsForMatrix(ORF_table[a,c(2,7)])
estimateSizeFactorsForMatrix(ORF_table[,c(2,7)])
sfp1_norm_factor = n1[2]/n1[1]
#plot(log2(ORF_table[,c(2,7)]+1), xlab = 'WT pilot', main = 'Raw Log2 NET-seq reads',
#     col = rgb(0,0,0,.5), pch = 16, ylab = expression(paste(italic('sfp1'), Delta)))
#points(log2(ORF_table[a,c(2,7)]+1), col = rgb(1,0,0,.5), pch = 16)
#legend(x = 'bottomright', col = c(1,2), legend = c('Not used', paste0('Used ', '(', length(a), ')')),
#       title = 'Normalization', pch = 16)

xrn1_norm_names = GRO_ref$ORF_Name[abs(log2(XRN1_GRO$VALUE+1)-log2(WT_GRO$VALUE+1)) <= .25]
a = match(xrn1_norm_names, ORF_annot$name)
a = a[!is.na(a)]
n1 = estimateSizeFactorsForMatrix(ORF_table[a,c(2,8,9)])
xrn1_norm_factor = n1[2]/n1[1]
xrn1_pilot_norm_factor = n1[3]/n1[1]
plot(log2(ORF_table[,c(2,8)]+1), xlab = 'WT pilot', main = 'Raw Log2 NET-seq reads',
     col = rgb(0,0,0,.5), pch = 16, ylab = expression(paste(italic('xrn1'), Delta)))
points(log2(ORF_table[a,c(2,8)]+1), col = rgb(1,0,0,.5), pch = 16)
legend(x = 'bottomright', col = c(1,2), legend = c('Not used', paste0('Used ', '(', length(a), ')')),
       title = 'Normalization', pch = 16)
plot(log2(ORF_table[,c(2,9)]+1), xlab = 'WT pilot', main = 'Raw Log2 NET-seq reads',
     col = rgb(0,0,0,.5), pch = 16, ylab = expression(paste(italic('xrn1'), Delta, ' Pilot')))
points(log2(ORF_table[a,c(2,9)]+1), col = rgb(1,0,0,.5), pch = 16)
legend(x = 'bottomright', col = c(1,2), legend = c('Not used', paste0('Used ', '(', length(a), ')')),
       title = 'Normalization', pch = 16)

n1 = estimateSizeFactorsForMatrix(ORF_table[,1:2])
WT_norm_factor = n1[1]/n1[2]
WT_pilot_norm_factor = n1[2]/n1[2]

dst1_norm_names = GRO_ref$ORF_Name[abs(log2(DST1_GRO$VALUE+1)-log2(WT2_GRO$VALUE+1)) <= .25]
a = match(dst1_norm_names, ORF_annot_c$name)
a = a[!is.na(a)]
n1 = estimateSizeFactorsForMatrix(their_sense_reads[a,c(7,10)])
estimateSizeFactorsForMatrix(their_sense_reads[a,c(7,10)])
dst1_norm_factor = n1[2]/n1[1]

norm_factors = c(WT_norm_factor, WT_pilot_norm_factor, ccr4_norm_factor, dhh1_norm_factor, 
             lsm1_norm_factor, rpb4_norm_factor, sfp1_norm_factor, xrn1_norm_factor,
             xrn1_pilot_norm_factor)





z = -log2(dst1_norm_factor) + log2(their_sense_reads$DST1+1) - log2(their_sense_reads$WT_NC)

plot(-100:500, rowMeans(WT_NC_TSS[, z > .5]), typ = 'l', ylim = c(0, 8), main = 'TSS metagenes',
     xlab = 'Distance from TSS', ylab = 'Mean NET-seq reads (A.U.)')
lines(-100:500, rowMeans(DST1_TSS[, z >.5])/dst1_norm_factor, col = 2)
legend(x = 'topright', lwd = 1, col = 1:2, legend = c('WT', 'dst1D'))

plot(-100:500, rowMeans(WT_NC_TSS[, z <= -.5]), typ = 'l', ylim = c(0, 8), main = 'TSS metagenes',
     xlab = 'Distance from TSS', ylab = 'Mean NET-seq reads (A.U.)')
lines(-100:500, rowMeans(DST1_TSS[, z <= -.5])/dst1_norm_factor, col = 2)
legend(x = 'topright', lwd = 1, col = 1:2, legend = c('WT', 'dst1D'))

par(mfrow = c(2,3))
col_vec = 1 + (abs(log2(XRN1_GRO$VALUE+1) - log2(WT_GRO$VALUE+1)) <= .25)
plot(log2(WT_GRO$VALUE+1), log2(XRN1_GRO$VALUE+1), main = 'GRO', xlab = 'WT', ylab = expression(paste(italic('xrn1'), Delta)),
     col = col_vec)
legend(x = 'bottomright', col = c(1,2), legend = c('Not eligible', 'Eligible'),
       title = 'Normalization', pch = 1)

col_vec = 1 + (abs(cdta[,105]) <= .25)
plot(log2(cdta[,3]+1), log2(cdta[,101]+1), main = 'CDTA', xlab = 'WT', ylab = expression(paste(italic('dhh1'), Delta)),
     col = col_vec)
legend(x = 'bottomright', col = c(1,2), legend = c('Not eligible', 'Eligible'),
       title = 'Normalization', pch = 1)

col_vec = 1 + (abs(cdta[,189]) <= .25)
plot(log2(cdta[,3]+1), log2(cdta[,185]+1), main = 'CDTA', xlab = 'WT', ylab = expression(paste(italic('lsm1'), Delta)),
     col = col_vec)
legend(x = 'bottomright', col = c(1,2), legend = c('Not eligible', 'Eligible'),
       title = 'Normalization', pch = 1)

col_vec = 1 + (abs(cdta[,45]) <= .25)
plot(log2(cdta[,3]+1), log2(cdta[,41]+1), main = 'CDTA', xlab = 'WT', ylab = expression(paste(italic('ccr4'), Delta)),
     col = col_vec)
legend(x = 'bottomright', col = c(1,2), legend = c('Not eligible', 'Eligible'),
       title = 'Normalization', pch = 1)

col_vec = 1 + (abs(log2(RPB4_GRO$VALUE+1) - log2(WT2_GRO$VALUE+1)) <= .25)
plot(log2(WT2_GRO$VALUE+1), log2(RPB4_GRO$VALUE+1), main = 'GRO', xlab = 'WT', ylab = expression(paste(italic('rpb4'), Delta)),
     col = col_vec)
legend(x = 'bottomright', col = c(1,2), legend = c('Not eligible', 'Eligible'),
       title = 'Normalization', pch = 1)



t0 = read.delim('data/gro_auxin_t0.txt')
t30 = read.delim('data/gro_auxin_t30.txt')

hist(log2(t30$VALUE+1)-log2(t0$VALUE+1))

fc = cbind.data.frame(t0$ID_REF ,log2(t30$VALUE+1)-log2(t0$VALUE+1))
fc = fc[!is.na(fc[,2]),]

fc = cbind.data.frame(GRO_ref$ORF_Name[match(fc$`t0$ID_REF`, GRO_ref$ID)], GRO_ref$Gene_Name[match(fc$`t0$ID_REF`, GRO_ref$ID)], fc)
colnames(fc) = c('ORF_Name', 'Common_name', 'ID_REF', 'FC')

fc_down = fc$ORF_Name[order(fc$FC, decreasing = F)]
write.table(fc_down, 'auxin_down.txt', quote = F, col.names = F, row.names = F)

xrn1_sorted_down = rownames(ORF_fc)[order(ORF_fc$XRN1, decreasing = F)]
write.table(xrn1_sorted_down, 'xrn1_sorted_down.txt', quote = F, col.names = F, row.names = F)



#fc = cbind.data.frame(,fc)


fc_up = fc$ORF_Name[order(fc$FC, decreasing = T)]
write.table(fc_up, 'auxin_up.txt', quote = F, col.names = F, row.names = F)

xrn1_sorted_up = rownames(ORF_fc)[order(ORF_fc$XRN1, decreasing = T)]
write.table(xrn1_sorted_up, 'xrn1_sorted_up.txt', quote = F, col.names = F, row.names = F)



go_ko = read.csv(file = '~/Downloads/go_ko.csv', sep = '\t')[,1:2]
go_aux = read.csv(file = '~/Downloads/go_aux.csv', sep = '\t')

go_merge = merge(go_ko, go_aux, by = 'GO.Term')



pdf(file = '~/Documents/R_Projects/NETSEQ/paper/DF2/Figures/supplement/p_vals_scatter_go.pdf', width = 10, height = 8)
par(mar = c(5,5,5,3))
plot(-log10(go_merge$P.value.x), -log10(go_merge$P.value.y), pch = 16,
     xlab = 'NET-seq, XRN1 KO', ylab = 'GRO, Xrn1-AID', main = '-log10 p-values for genes stimulated by Xrn1',
     cex.main = 2, cex.axis = 2, cex.lab = 2, cex = 1.5)
dev.off()


fc_net = cbind.data.frame(ORF_annot$name, ORF_fc$XRN1)
colnames(fc_net) = c('ORF_Name', 'FC_NET')

a = merge(fc_net, fc)

plot(a$FC_NET, a$FC)
cor(a$FC_NET, a$FC, use = 'complete.obs')
cor(a$FC_NET, a$FC, use = 'complete.obs', method = 'spearman')


go_ko_up = read.csv(file = '~/Downloads/go_ko_up.csv', sep = '\t')
go_aux_up = read.csv(file = '~/Downloads/go_aux_up.csv', sep = '\t')

go_merge_up = merge(go_ko_up, go_aux_up, by = 'GO.Term')
plot(-log10(go_merge_up$P.value.x), -log10(go_merge_up$P.value.y), pch = 16,
     xlab = '-log10 P-value NET-seq KO', ylab = '-log10 P-value GRO Auxin')


mHG.test(1 * (fc$ORF_Name[order(fc$FC, decreasing = T)] %in% ar2))$p.value

mHG.test(1 * (fc$ORF_Name[order(fc$FC, decreasing = T)] %in% z$V2))$p.value
mHG.test(1 * (fc$ORF_Name[order(fc$FC, decreasing = F)] %in% z$V2))$p.value
#mHG.test(1 * (fc$Common_name[order(fc$FC, decreasing = T)] %in% gp2))$p.value
#mHG.test(1 * (fc$Common_name[order(fc$FC, decreasing = F)] %in% gp2))$p.value

mHG.test(1 * (fc$ORF_Name[order(fc$FC, decreasing = T)] %in% y$V2))$p.value
mHG.test(1 * (fc$ORF_Name[order(fc$FC, decreasing = F)] %in% y$V2))$p.value
mHG.test(1 * (fc$Common_name[order(fc$FC, decreasing = T)] %in% ar2))$p.value
mHG.test(1 * (fc$Common_name[order(fc$FC, decreasing = F)] %in% ar2))$p.value


mHG.test(1 * (fc$Common_name[order(fc$FC, decreasing = F)] %in% ar2))$p.value
mHG.test(1 * (fc$Common_name[order(fc$FC, decreasing = F)] %in% gp2))$p.value


zz = c(mHG.test(1 * (ORF_annot$commonName[order(ORF_fc[,5], decreasing = F)] %in% ar2))$p.value,
mHG.test(1 * (ORF_annot$commonName[order(ORF_fc[,5], decreasing = F)] %in% gp2))$p.value,
mHG.test(1 * (fc$Common_name[order(fc$FC, decreasing = F)] %in% ar2))$p.value,
mHG.test(1 * (fc$Common_name[order(fc$FC, decreasing = F)] %in% gp2))$p.value,

mHG.test(1 * (ORF_annot$commonName[order(ORF_fc[,5], decreasing = T)] %in% ar2))$p.value,
mHG.test(1 * (ORF_annot$commonName[order(ORF_fc[,5], decreasing = T)] %in% gp2))$p.value,
mHG.test(1 * (fc$Common_name[order(fc$FC, decreasing = T)] %in% ar2))$p.value,
mHG.test(1 * (fc$Common_name[order(fc$FC, decreasing = T)] %in% gp2))$p.value)

pdf(file = '~/Documents/R_Projects/NETSEQ/paper/DF2/Figures/supplement/aerobic_glycolysis_xrn1_2.pdf', width = 10, height = 8)
par(xpd=FALSE)
par(mar = c(5,5,5,3))
barplot(-log10(zz), space = c(1,0,1,0,2,0,1,0),
        ylab = '-log10 p-value', xlab = 'Gene response to Xrn1', main = 'Aerobic respiration and glycolysis enrichment',
        col = c('blue', 'red'), cex.axis = 2, cex.lab = 2, cex.main = 2)
par(xpd=TRUE)
text(x = c(2,5,9,12), y = -.25, c('NET-seq, XRN1 KO', 'GRO, Xrn1-AID', 'NET-seq, XRN1 KO', 'GRO, Xrn1-AID'), cex = 1.25)
text(x = c(3.5,10.5), y = -.65, c('Stimulated', 'Repressed'), cex = 1.5)
legend(x = .8, y = 8, fill = c('blue', 'red'), legend = c('aerobic respiration', 'glycolytic process'), cex = 1.5)
dev.off()