end_metagene_plotter = function(gene_ind, gene_title, which_end, which_strain = 'all', sm_val = .05){
  
  if(which_end == 'TSS'){
    WT_mean = lowess(rowMeans(WT_TSS[,gene_ind]+WT_pilot_TSS[,gene_ind])/(nf[1]+nf[2]), f = sm_val)$y
    CCR4_mean = lowess(rowMeans(CCR4_TSS[,gene_ind])/nf[3], f = sm_val)$y
    DHH1_mean = lowess(rowMeans(DHH1_TSS[,gene_ind])/nf[4], f = sm_val)$y
    LSM1_mean = lowess(rowMeans(LSM1_TSS[,gene_ind])/nf[5], f = sm_val)$y
    RPB4_mean = lowess(rowMeans(RPB4_TSS[,gene_ind])/nf[6], f = sm_val)$y
    XRN1_mean = lowess(rowMeans(XRN1_TSS[,gene_ind]+XRN1_pilot_TSS[,gene_ind])/(nf[8]+nf[9]), f = sm_val)$y
    
    y_max = max(WT_mean, CCR4_mean, DHH1_mean, LSM1_mean, RPB4_mean, XRN1_mean)
    
    if(which_strain == 'all'){
      par(mfrow = c(2,3))
      flags = rep(T,5)
    } else {
      #par(mfrow = c(1,1))
      flags = rep(0,5)
      if(which_strain == 'xrn1'){flags[1]=T}
      if(which_strain == 'dhh1'){flags[2]=T}
      if(which_strain == 'lsm1'){flags[3]=T}
      if(which_strain == 'ccr4'){flags[4]=T}
      if(which_strain == 'rpb4'){flags[5]=T}
    }
    
    par(mar = c(5,5,5,3))
    
    if(flags[1]){
      plot(-100:500, WT_mean, typ = 'l', ylab = 'Mean normalized reads',
           xlab = 'Distance from TSS', ylim = c(0, y_max), lwd = 2, main = paste('TSS metagenes,', gene_title), cex.lab = 1.75, cex.main = 1.75, cex.axis = 1.6)
      lines(-100:500, XRN1_mean, col = 2, lwd = 2)
      legend(x = 'topright', legend = c('WT', expression(paste(italic("xrn1"), Delta))), col = c(1,2), lwd = 2, cex = 1.5)
    }

    if(flags[2]){
    plot(-100:500, WT_mean, typ = 'l', ylab = 'Mean normalized reads',
         xlab = 'Distance from TSS', ylim = c(0, y_max), lwd = 2, main = paste('TSS metagenes,', gene_title), cex.lab = 1.75, cex.main = 1.75, cex.axis = 1.6)
    lines(-100:500, DHH1_mean, col = 2, lwd = 2)
    legend(x = 'topright', legend = c('WT', expression(paste(italic("dhh1"), Delta))), col = c(1,2), lwd = 2, cex = 1.5)
    }
    
    if(flags[3]){
    plot(-100:500, WT_mean, typ = 'l', ylab = 'Mean normalized reads',
         xlab = 'Distance from TSS', ylim = c(0, y_max), lwd = 2, main = paste('TSS metagenes,', gene_title), cex.lab = 1.75, cex.main = 1.75, cex.axis = 1.6)
    lines(-100:500, LSM1_mean, col = 2, lwd = 2)
    legend(x = 'topright', legend = c('WT', expression(paste(italic("lsm1"), Delta))), col = c(1,2), lwd = 2, cex = 1.5)
    }
    
    if(flags[4]){
    plot(-100:500, WT_mean, typ = 'l', ylab = 'Mean normalized reads',
         xlab = 'Distance from TSS', ylim = c(0, y_max), lwd = 2, main = paste('TSS metagenes,', gene_title), cex.lab = 1.75, cex.main = 1.75, cex.axis = 1.6)
    lines(-100:500, CCR4_mean, col = 2, lwd = 2)
    legend(x = 'topright', legend = c('WT', expression(paste(italic("ccr4"), Delta))), col = c(1,2), lwd = 2, cex = 1.5)
    }
    
    if(flags[5]){
    plot(-100:500, WT_mean, typ = 'l', ylab = 'Mean normalized reads',
         xlab = 'Distance from TSS', ylim = c(0, y_max), lwd = 2, main = paste('TSS metagenes,', gene_title), cex.lab = 1.75, cex.main = 1.75, cex.axis = 1.6)
    lines(-100:500, RPB4_mean, col = 2, lwd = 2)
    legend(x = 'topright', legend = c('WT', expression(paste(italic("rpb4"), Delta))), col = c(1,2), lwd = 2, cex = 1.5)
    }
    
  } else if(which_end == 'PAS'){
    
    WT_mean = lowess(rowMeans(WT_TES[,gene_ind]+WT_pilot_TES[,gene_ind])/(nf[1]+nf[2]), f = .05)$y
    CCR4_mean = lowess(rowMeans(CCR4_TES[,gene_ind])/nf[3], f = .05)$y
    DHH1_mean = lowess(rowMeans(DHH1_TES[,gene_ind])/nf[4], f = .05)$y
    LSM1_mean = lowess(rowMeans(LSM1_TES[,gene_ind])/nf[5], f = .05)$y
    RPB4_mean = lowess(rowMeans(RPB4_TES[,gene_ind])/nf[6], f = .05)$y
    XRN1_mean = lowess(rowMeans(XRN1_TES[,gene_ind]+XRN1_pilot_TES[,gene_ind])/(nf[8]+nf[9]), f = .05)$y
    
    y_max = max(WT_mean, CCR4_mean, DHH1_mean, LSM1_mean, RPB4_mean, XRN1_mean)
    
    if(which_strain == 'all'){
      par(mfrow = c(2,3))
      flags = rep(T,5)
    } else {
      #par(mfrow = c(1,1))
      flags = rep(0,5)
      if(which_strain == 'xrn1'){flags[1]=T}
      if(which_strain == 'dhh1'){flags[2]=T}
      if(which_strain == 'lsm1'){flags[3]=T}
      if(which_strain == 'ccr4'){flags[4]=T}
      if(which_strain == 'rpb4'){flags[5]=T}
    }
    
    par(mar = c(5,5,5,3))
    
    if(flags[1]){
    plot(-150:150, WT_mean, typ = 'l', ylab = 'Mean normalized reads',
         xlab = 'Distance from PAS', ylim = c(0, y_max), lwd = 2, main = paste('PAS metagenes,', gene_title), cex.lab = 1.75, cex.main = 1.75, cex.axis = 1.6)
    lines(-150:150, XRN1_mean, col = 2, lwd = 2)
    legend(x = 'topright', legend = c('WT', expression(paste(italic("xrn1"), Delta))), col = c(1,2), lwd = 2, cex = 1.5)
    }
    
    if(flags[2]){
    plot(-150:150, WT_mean, typ = 'l', ylab = 'Mean normalized reads',
         xlab = 'Distance from PAS', ylim = c(0, y_max), lwd = 2, main = paste('PAS metagenes,', gene_title), cex.lab = 1.75, cex.main = 1.75, cex.axis = 1.6)
    lines(-150:150, DHH1_mean, col = 2, lwd = 2)
    legend(x = 'topright', legend = c('WT', expression(paste(italic("dhh1"), Delta))), col = c(1,2), lwd = 2, cex = 1.5)
    }
  
    if(flags[3]){
    plot(-150:150, WT_mean, typ = 'l', ylab = 'Mean normalized reads',
         xlab = 'Distance from PAS', ylim = c(0, y_max), lwd = 2, main = paste('PAS metagenes,', gene_title), cex.lab = 1.75, cex.main = 1.75, cex.axis = 1.6)
    lines(-150:150, LSM1_mean, col = 2, lwd = 2)
    legend(x = 'topright', legend = c('WT', expression(paste(italic("lsm1"), Delta))), col = c(1,2), lwd = 2, cex = 1.5)
    }
    
    if(flags[4]){
    plot(-150:150, WT_mean, typ = 'l', ylab = 'Mean normalized reads',
         xlab = 'Distance from PAS', ylim = c(0, y_max), lwd = 2, main = paste('PAS metagenes,', gene_title), cex.lab = 1.75, cex.main = 1.75, cex.axis = 1.6)
    lines(-150:150, CCR4_mean, col = 2, lwd = 2)
    legend(x = 'topright', legend = c('WT', expression(paste(italic("ccr4"), Delta))), col = c(1,2), lwd = 2, cex = 1.5)
    }
    
    if(flags[5]){
    plot(-150:150, WT_mean, typ = 'l', ylab = 'Mean normalized reads',
         xlab = 'Distance from PAS', ylim = c(0, y_max), lwd = 2, main = paste('PAS metagenes,', gene_title), cex.lab = 1.75, cex.main = 1.75, cex.axis = 1.6)
    lines(-150:150, RPB4_mean, col = 2, lwd = 2)
    legend(x = 'topright', legend = c('WT', expression(paste(italic("rpb4"), Delta))), col = c(1,2), lwd = 2, cex = 1.5)
    }
  }
  
}


######
# Figure 1 - ORF FC distributions 

a1 = density(ORF_fc$CCR4, na.rm = T)
a2 = density(ORF_fc$DHH1, na.rm = T)
a3 = density(ORF_fc$LSM1, na.rm = T)
a4 = density(ORF_fc$RPB4, na.rm = T)
a5 = density(ORF_fc$XRN1, na.rm = T)

pdf(file = 'xrn1_fc_density.pdf')
par(mar = c(5,5,5,3))
plot(a5, xlim = c(-5,5), xlab = 'Standardized Log2 NET-seq FC', main = 'Gene fold change distribution(s)',
     col = 1, lwd = 2, ylim = c(0,.5), cex.lab = 2, cex.main = 2, cex.axis = 2)
abline(v = 0, col = 1, lty = 2, lwd = 2)
legend(x = 'topright', legend = c(expression(italic(paste("xrn1", Delta)))),
       col = 1, lwd = 2, cex = 2)
dev.off()
pdf(file = 'other_fc_density.pdf')
par(mar = c(5,5,5,3))
plot(a2, xlim = c(-5,5), xlab = 'Standardized Log2 NET-seq FC', main = 'Gene fold change distribution(s)',
     col = 2, lwd = 2, ylim = c(0,.5), cex.lab = 2, cex.main = 2, cex.axis = 2)
lines(a3,col = 3, lwd = 2)
lines(a4,col = 4, lwd = 2)
abline(v = 0, col = 1, lty = 2, lwd = 2)
legend(x = 'topright', legend = c(expression(italic(paste("ccr4", Delta))), expression(italic(paste("dhh1", Delta))), expression(italic(paste("lsm1", Delta)))),
       col = 1:3, lwd = 2, cex = 2)
dev.off()

pdf(file = 'transcriptional_efficiency.pdf')
par(mar = c(5,5,5,3))
a1 = log2(merged_mat_2$`XRN1_GRO$VALUE`) - log2(merged_mat_2$XRN1) - (log2(merged_mat_2$`WT_GRO$VALUE`) - log2(merged_mat_2$WT))
b1 = density(a1, na.rm = T)
a2 = log2(merged_mat_2$`RPB4_GRO$VALUE`) - log2(merged_mat_2$RPB4) - (log2(merged_mat_2$`WT2_GRO$VALUE`) - log2(merged_mat_2$WT))
b2 = density(a2, na.rm = T)
plot(b1, lwd = 2, col = 1, xlab = 'Log2 GRO/NET-seq FC',
     main = 'Changes in elongation efficiency', cex.lab = 2, cex.axis = 2, cex.main = 2)
lines(b2, lwd = 2, col = 6)
legend(x = 'topleft', col = c(1,6), lwd = 2, legend = c(expression(paste(italic('xrn1'), Delta)), expression(paste(italic('rpb4'), Delta))), cex = 1.5)
abline(v = 0, lty = 2, lwd = 2)
dev.off()

######
# Figure 2 - Comparison with Xrn1 responsiveness

xrn1_rating_list = lapply(2:10, function(x) which(ORF_annot$name %in% xrn1_synth_ratings$ORF.name[xrn1_synth_ratings$responsiveness == x]))

pdf(file = 'xrn1_box.pdf', width = 7, height = 6)
par(mar = c(5,5,5,3))
boxplot(ORF_fc[xrn1_rating_list[[1]],5], ORF_fc[xrn1_rating_list[[2]],5], ORF_fc[xrn1_rating_list[[3]],5], 
        ORF_fc[xrn1_rating_list[[4]],5], ORF_fc[xrn1_rating_list[[5]],5], ORF_fc[xrn1_rating_list[[6]],5], 
        ORF_fc[xrn1_rating_list[[7]],5], ORF_fc[xrn1_rating_list[[8]],5], ORF_fc[xrn1_rating_list[[9]],5],
        ylim = c(-5,5), ylab = 'Standardized NET-seq FC', xlab = 'Xrn1 responsiveness', names = 2:10, main = expression(italic(paste('xrn1', Delta))),
        cex.main = 3, cex.lab = 2, cex.axis = 1.95)
abline(h = 0, lwd = 2, lty = 2, col = 2)
dev.off()

pdf(file = 'dhh1_box.pdf', width = 7, height = 6)
par(mar = c(5,5,5,3))
boxplot(ORF_fc[xrn1_rating_list[[1]],2], ORF_fc[xrn1_rating_list[[2]],2], ORF_fc[xrn1_rating_list[[3]],2], 
        ORF_fc[xrn1_rating_list[[4]],2], ORF_fc[xrn1_rating_list[[5]],2], ORF_fc[xrn1_rating_list[[6]],2], 
        ORF_fc[xrn1_rating_list[[7]],2], ORF_fc[xrn1_rating_list[[8]],2], ORF_fc[xrn1_rating_list[[9]],2],
        ylim = c(-5,5), ylab = 'Standardized NET-seq FC', xlab = 'Xrn1 responsiveness', names = 2:10, main = expression(italic(paste('dhh1', Delta))),
        cex.main = 3, cex.lab = 2, cex.axis = 1.95)
abline(h = 0, lwd = 2, lty = 2, col = 2)
dev.off()

pdf(file = 'lsm1_box.pdf', width = 7, height = 6)
par(mar = c(5,5,5,3))
boxplot(ORF_fc[xrn1_rating_list[[1]],3], ORF_fc[xrn1_rating_list[[2]],3], ORF_fc[xrn1_rating_list[[3]],3], 
        ORF_fc[xrn1_rating_list[[4]],3], ORF_fc[xrn1_rating_list[[5]],3], ORF_fc[xrn1_rating_list[[6]],3], 
        ORF_fc[xrn1_rating_list[[7]],3], ORF_fc[xrn1_rating_list[[8]],3], ORF_fc[xrn1_rating_list[[9]],3],
        ylim = c(-5,5), ylab = 'Standardized NET-seq FC', xlab = 'Xrn1 responsiveness', names = 2:10, main = expression(italic(paste('lsm1', Delta))),
        cex.main = 3, cex.lab = 2, cex.axis = 1.95)
abline(h = 0, lwd = 2, lty = 2, col = 2)
dev.off()

pdf(file = 'ccr4_box.pdf', width = 7, height = 6)
par(mar = c(5,5,5,3))
boxplot(ORF_fc[xrn1_rating_list[[1]],1], ORF_fc[xrn1_rating_list[[2]],1], ORF_fc[xrn1_rating_list[[3]],1], 
        ORF_fc[xrn1_rating_list[[4]],1], ORF_fc[xrn1_rating_list[[5]],1], ORF_fc[xrn1_rating_list[[6]],1], 
        ORF_fc[xrn1_rating_list[[7]],1], ORF_fc[xrn1_rating_list[[8]],1], ORF_fc[xrn1_rating_list[[9]],1],
        ylim = c(-5,5), ylab = 'Standardized NET-seq FC', xlab = 'Xrn1 responsiveness', names = 2:10, main = expression(italic(paste('ccr4', Delta))),
        cex.main = 3, cex.lab = 2, cex.axis = 1.95)
abline(h = 0, lwd = 2, lty = 2, col = 2, cex = 2)
dev.off()

pdf(file = 'rpb4_box.pdf', width = 7, height = 6)
par(mar = c(5,5,5,3))
boxplot(ORF_fc[xrn1_rating_list[[1]],4], ORF_fc[xrn1_rating_list[[2]],4], ORF_fc[xrn1_rating_list[[3]],4], 
        ORF_fc[xrn1_rating_list[[4]],4], ORF_fc[xrn1_rating_list[[5]],4], ORF_fc[xrn1_rating_list[[6]],4], 
        ORF_fc[xrn1_rating_list[[7]],4], ORF_fc[xrn1_rating_list[[8]],4], ORF_fc[xrn1_rating_list[[9]],4],
        ylim = c(-5,5), ylab = 'Standardized NET-seq FC', xlab = 'Xrn1 responsiveness', names = 2:10, main = expression(italic(paste('rpb4', Delta))),
        cex.main = 3, cex.lab = 2, cex.axis = 1.95)
abline(h = 0, lwd = 2, lty = 2, col = 2)
dev.off()

######
### Figure 4 - Full metagenes ###

pdf(file = 'full_metagenes_1_smoothed.pdf', width = 14, height = 5)
par(mar = c(5,5,5,3))
plot(WT_smoothed$x, WT_smoothed$y/sum(nf[1:2]), col = 'gray', lwd = 3, typ = 'l', ylim = c(0, 1.2),
     xlab = 'Metaposition', ylab = 'Mean NET-seq reads', main = 'Full metagenes',
     cex.main = 2, cex.axis = 2, cex.lab = 2)
lines(xrn1_smoothed$x, xrn1_smoothed$y/sum(nf[8:9]), col = 1, lwd = 3)
legend(x = 'topright', col = c('gray', 1), lwd = 3, legend = c('WT', expression(paste(italic('xrn1'), Delta))), cex = 1.5)
dev.off()

pdf(file = 'full_metagenes_2_smoothed.pdf', width = 14, height = 5)
par(mar = c(5,5,5,3))
plot(WT_smoothed$x, WT_smoothed$y/sum(nf[1:2]), col = 'gray', lwd = 3, typ = 'l', ylim = c(0, 1.2),
     xlab = 'Metaposition', ylab = 'Mean NET-seq reads', main = 'Full metagenes',
     cex.main = 2, cex.axis = 2, cex.lab = 2)
lines(dhh1_smoothed$x, dhh1_smoothed$y/nf[4], col = 2, lwd = 3)
lines(lsm1_smoothed$x, lsm1_smoothed$y/nf[5], col = 3, lwd = 3)
legend(x = 'topright', col = c('gray', 2, 3), lwd = 3, legend = c('WT', expression(paste(italic('dhh1'), Delta)), expression(paste(italic('lsm1'), Delta))), cex = 1.5)
dev.off()

pdf(file = 'full_metagenes_3_smoothed.pdf', width = 14, height = 5)
par(mar = c(5,5,5,3))
plot(WT_smoothed$x, WT_smoothed$y/sum(nf[1:2]), col = 'gray', lwd = 3, typ = 'l', ylim = c(0, 1.2),
     xlab = 'Metaposition', ylab = 'Mean NET-seq reads', main = 'Full metagenes',
     cex.main = 2, cex.axis = 2, cex.lab = 2)
lines(ccr4_smoothed$x, ccr4_smoothed$y/nf[3], col = 4, lwd = 3)
lines(rpb4_smoothed$x, rpb4_smoothed$y/nf[6], col = 6, lwd = 3)
legend(x = 'topright', col = c('gray', 4, 6), lwd = 3, legend = c('WT', expression(paste(italic('ccr4'), Delta)), expression(paste(italic('rpb4'), Delta))), cex = 1.5)
dev.off()

######
### Figure 5 - TSS and PAS metagenes ###

pdf(file = 'TSS_metagenes_all.pdf', width = 14)
end_metagene_plotter(gene_ind = 1:4973, gene_title = 'All genes', which_end = 'TSS')
dev.off()

pdf(file = 'PAS_metagenes_all.pdf', width = 14)
end_metagene_plotter(gene_ind = 1:4973, gene_title = 'All genes', which_end = 'PAS')
dev.off()

######
### Figure 7 - Convergent/Divergent Pol II ###
pdf(file = 'convergent_pol_ii.pdf', width = 14)
par(mar = c(5,5,5,3))
par(mfrow = c(2,3))
plot(-500:500, lowess(wt_mean, f = .1)$y, typ = 'l', ylim = c(0,.35), lwd = 3, 
     xlab = 'Distance from midpoint', ylab = 'Average reads', main = expression(paste(italic('xrn1'), Delta)),
     cex.main = 2, cex.axis = 1.75, cex.lab = 2)
lines(-500:500, lowess(xrn1_mean, f = .1)$y, col = 'red', lwd = 3)
legend(x = 'bottomleft', col = 1:2, lwd = 3, legend = c('WT', expression(paste(italic('xrn1'), Delta))), cex = 1.5, ncol = 2)

plot(-500:500, lowess(wt_mean, f = .1)$y, typ = 'l', ylim = c(0,.35), lwd = 3, 
     xlab = 'Distance from midpoint', ylab = 'Average reads', main = expression(paste(italic('dhh1'), Delta)),
     cex.main = 2, cex.axis = 1.75, cex.lab = 2)
lines(-500:500, lowess(dhh1_mean, f = .1)$y, col = 'red', lwd = 3)
legend(x = 'bottomleft', col = 1:2, lwd = 3, legend = c('WT', expression(paste(italic('dhh1'), Delta))), cex = 1.5, ncol = 2)

plot(-500:500, lowess(wt_mean, f = .1)$y, typ = 'l', ylim = c(0,.35), lwd = 3, 
     xlab = 'Distance from midpoint', ylab = 'Average reads', main = expression(paste(italic('lsm1'), Delta)),
     cex.main = 2, cex.axis = 1.75, cex.lab = 2)
lines(-500:500, lowess(lsm1_mean, f = .1)$y, col = 'red', lwd = 3)
legend(x = 'bottomleft', col = 1:2, lwd = 3, legend = c('WT', expression(paste(italic('lsm1'), Delta))), cex = 1.5, ncol = 2)

plot(-500:500, lowess(wt_mean, f = .1)$y, typ = 'l', ylim = c(0,.35), lwd = 3, 
     xlab = 'Distance from midpoint', ylab = 'Average reads', main = expression(paste(italic('ccr4'), Delta)),
     cex.main = 2, cex.axis = 1.75, cex.lab = 2)
lines(-500:500, lowess(ccr4_mean, f = .1)$y, col = 'red', lwd = 3)
legend(x = 'bottomleft', col = 1:2, lwd = 3, legend = c('WT', expression(paste(italic('ccr4'), Delta))), cex = 1.5, ncol = 2)

plot(-500:500, lowess(wt_mean, f = .1)$y, typ = 'l', ylim = c(0,.35), lwd = 3, 
     xlab = 'Distance from midpoint', ylab = 'Average reads', main = expression(paste(italic('rpb4'), Delta)),
     cex.main = 2, cex.axis = 1.75, cex.lab = 2)
lines(-500:500, lowess(rpb4_mean, f = .1)$y, col = 'red', lwd = 3)
legend(x = 'bottomleft', col = 1:2, lwd = 3, legend = c('WT', expression(paste(italic('rpb4'), Delta))), cex = 1.5, ncol = 2)

plot(x = 0, y = 0, col = 'white', axes = F, xlab = '', ylab = '', main = 'Convergent gene orientation', ylim = c(-1,1))
abline(h = c(.5, -.5))
#lines(c(-1,-.35), c(.5,.5), lwd = 10, col = 2)
arrows(-1, .5, -.25, .5, lwd = 5, col = 2, length = .2)
arrows(1, -.5, .25, -.5, lwd = 5, col = 2, length = .2)
lines(c(-1, -1), c(.5, .75))
arrows(-1, .75, -.75, .75, length = .2)
lines(c(1, 1), c(-.5, -.75))
arrows(1, -.75, .75, -.75, length = .2)
dev.off()

pdf(file = 'divergent_pol_ii.pdf', width = 14)
par(mar = c(5,5,5,3))
par(mfrow = c(2,3))

plot(-500:500, lowess(rowMeans(WT_div+WT_pilot_div)/sum(norm_factors[1:2]), f = .1)$y, ylim = c(0,.4), main = expression(paste(italic('xrn1'), Delta)),
     xlab = 'Distance from midpoint', ylab = 'Normalized mean reads', typ = 'l', lwd = 2, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
lines(-500:500, lowess(rowMeans(XRN1_div+XRN1_pilot_div)/sum(norm_factors[8:9]), f = .1)$y, col = 'red', lwd = 2)
legend(x = 'topleft', col = 1:2, lwd = 2, legend = c('WT', expression(paste(italic('xrn1'), Delta))), cex = 1.5, ncol = 2)

plot(-500:500, lowess(rowMeans(WT_div+WT_pilot_div)/sum(norm_factors[1:2]), f = .1)$y, ylim = c(0,.4), main = expression(paste(italic('dhh1'), Delta)),
     xlab = 'Distance from midpoint', ylab = 'Normalized mean reads', typ = 'l', lwd = 2, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
lines(-500:500, lowess(rowMeans(DHH1_div)/norm_factors[4], f = .1)$y, col = 'red', lwd = 2)
legend(x = 'topleft', col = 1:2, lwd = 2, legend = c('WT', expression(paste(italic('dhh1'), Delta))), cex = 1.5, ncol = 2)

plot(-500:500, lowess(rowMeans(WT_div+WT_pilot_div)/sum(norm_factors[1:2]), f = .1)$y, ylim = c(0,.4), main = expression(paste(italic('lsm1'), Delta)),
     xlab = 'Distance from midpoint', ylab = 'Normalized mean reads', typ = 'l', lwd = 2, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
lines(-500:500, lowess(rowMeans(LSM1_div)/norm_factors[5], f = .1)$y, col = 'red', lwd = 2)
legend(x = 'topleft', col = 1:2, lwd = 2, legend = c('WT', expression(paste(italic('lsm1'), Delta))), cex = 1.5, ncol = 2)

plot(-500:500, lowess(rowMeans(WT_div+WT_pilot_div)/sum(norm_factors[1:2]), f = .1)$y, ylim = c(0,.4), main = expression(paste(italic('ccr4'), Delta)),
     xlab = 'Distance from midpoint', ylab = 'Normalized mean reads', typ = 'l', lwd = 2, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
lines(-500:500, lowess(rowMeans(CCR4_div)/norm_factors[3], f = .1)$y, col = 'red', lwd = 2)
legend(x = 'topleft', col = 1:2, lwd = 2, legend = c('WT', expression(paste(italic('ccr4'), Delta))), cex = 1.5, ncol = 2)

plot(-500:500, lowess(rowMeans(WT_div+WT_pilot_div)/sum(norm_factors[1:2]), f = .1)$y, ylim = c(0,.4), main = expression(paste(italic('rpb4'), Delta)),
     xlab = 'Distance from midpoint', ylab = 'Normalized mean reads', typ = 'l', lwd = 2, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
lines(-500:500, lowess(rowMeans(RPB4_div)/norm_factors[6], f = .1)$y, col = 'red', lwd = 2)
legend(x = 'topleft', col = 1:2, lwd = 2, legend = c('WT', expression(paste(italic('rpb4'), Delta))), cex = 1.5, ncol = 2)

plot(x = 0, y = 0, col = 'white', axes = F, xlab = '', ylab = '', main = 'Divergent gene orientation', ylim = c(-1,1))
abline(h = c(.5, -.5))
#lines(c(-1,-.35), c(.5,.5), lwd = 10, col = 2)
arrows(-.25, .5, -1, .5, lwd = 5, col = 2, length = .2)
arrows(.25, -.5, 1, -.5, lwd = 5, col = 2, length = .2)
lines(c(-.25, -.25), c(.5, .75))
arrows(-.25, .75, -.5, .75, length = .2)
lines(c(.25, .25), c(-.5, -.75))
arrows(.25, -.75, .5, -.75, length = .2)
dev.off()

######
### Figure 8 - Comparing response in SAGA and TFIID genes ###

pdf(file = '~/Documents/R_Projects/NETSEQ/paper/DF2/Figures/main/FC_boxplot_SAGA_TFIID.pdf', width = 10)
par(mar = c(5,5,5,3))
boxplot(ORF_fc[SAGA_ind,5], ORF_fc[TFIID_ind,5],
        ORF_fc[SAGA_ind,2], ORF_fc[TFIID_ind,2],
        ORF_fc[SAGA_ind,3], ORF_fc[TFIID_ind,3],
        ORF_fc[SAGA_ind,1], ORF_fc[TFIID_ind,1],
        ORF_fc[SAGA_ind,4], ORF_fc[TFIID_ind,4],
        at = c(1.1,1.9,3.1,3.9,5.1,5.9,7.1,7.9,9.1,9.9),
        ylim = c(-5,5), ylab = 'Standardized NET-seq FC',
        main = 'Fold changes for SAGA- and TFIID-dominated genes',
        notch = F, col = c(4,2,4,2,4,2,4,2,4,2), xlab = 'Mutant strain',
        xaxt = "n",
        cex.axis = 2, cex.main = 2, cex.lab = 2)
axis(1, cex.axis = 2, at=c(1.5,3.5,5.5,7.5,9.5), labels=c(expression(paste(italic('xrn1'), Delta)), expression(paste(italic('dhh1'), Delta)), expression(paste(italic('lsm1'), Delta)), expression(paste(italic('ccr4'), Delta)), expression(paste(italic('rpb4'), Delta))))
abline(h = 0, lty = 2, lwd = 2)
legend(x = 'topright', fill = c(4,2), 
       legend = c('SAGA-dominated', 'TFIID-dominated'),
       horiz = T, cex = 1.2)
dev.off()

pdf(file = 'TSS_xrn1_SAGA.pdf', width = 10)
end_metagene_plotter(gene_ind = SAGA_ind, gene_title = 'SAGA-dominated genes', which_end = 'TSS', which_strain = 'xrn1')
dev.off()

pdf(file = 'PAS_xrn1_SAGA.pdf', width = 10)
end_metagene_plotter(gene_ind = SAGA_ind, gene_title = 'SAGA-dominated genes', which_end = 'PAS', which_strain = 'xrn1')
dev.off()

pdf(file = 'TSS_xrn1_TFIID.pdf', width = 10)
end_metagene_plotter(gene_ind = TFIID_ind, gene_title = 'TFIID-dominated genes', which_end = 'TSS', which_strain = 'xrn1')
dev.off()

pdf(file = 'PAS_xrn1_TFIID.pdf', width = 10)
end_metagene_plotter(gene_ind = TFIID_ind, gene_title = 'TFIID-dominated genes', which_end = 'PAS', which_strain = 'xrn1')
dev.off()

par(mar = c(5,5,5,3))
boxplot(rowMeans(log2(SKI2_CRAC_mat_norm[SAGA_ind,]+1))-rowMeans(log2(XRN1_CRAC_mat_norm[SAGA_ind,]+1)),
        rowMeans(log2(SKI2_CRAC_mat_norm[TFIID_ind,]+1))-rowMeans(log2(XRN1_CRAC_mat_norm[TFIID_ind,]+1)),
        log2(saga_decay$HL._xrn1) - log2(saga_decay$HL.WT),
        log2(tfiid_decay$HL._xrn1) - log2(tfiid_decay$HL.WT),
        at = c(1.1,1.9,3.1,3.9),
        ylim = c(-5,5), ylab = 'Log2 values',
        main = 'DF binding and mRNA HL FCs',
        notch = F, col = c(4,2,4,2), xlab = '',
        xaxt = "n",
        cex.axis = 2, cex.main = 2, cex.lab = 2)
axis(1, cex.axis = 2, at=c(1.5,3.5), labels=c('SKI2/XRN1 CRAC Ratio', expression(paste('mRNA HL FC in ', italic('xrn1'), Delta))))
abline(v = 2.5, lty = 1, lwd = 2)
legend(x = 'bottomright', fill = c(4,2), 
       legend = c('SAGA-dominated', 'TFIID-dominated'),
       horiz = T)


#######
### Figure 9 - ncRNA FCs
pdf(file = 'xrn1_fc_dens.pdf', width = 9)
i = 5
a1 = density(ORF_fc[,i], na.rm = T)
a2 = density(CUT_fc[,i], na.rm = T)
a3 = density(NUT_fc[,i], na.rm = T)
a4 = density(SUT_fc[,i], na.rm = T)
a5 = density(XUT_fc[,i], na.rm = T)
plot(a1, typ = 'l', xlim = c(-5,5), ylim = c(0,.5), lwd = 2,
     main = expression(italic(paste('xrn1', Delta))), xlab = 'Standardized NET-seq FC',
     cex.axis = 2, cex.lab = 2, cex.main = 3)
lines(a2, col = 2, lwd = 2)
lines(a3, col = 3, lwd = 2)
lines(a4, col = 4, lwd = 2)
lines(a5, col = 5, lwd = 2)
legend(x = 'topright', legend = c('Genes', 'CUTs', 'NUTs', 'SUTs', 'XUTs'),
       col = 1:5, lwd = 2, cex = 2)
abline(v = 0, lty = 2, lwd = 2)
dev.off()

pdf(file = 'dhh1_fc_dens.pdf', width = 9)
i = 2
a1 = density(ORF_fc[,i], na.rm = T)
a2 = density(CUT_fc[,i], na.rm = T)
a3 = density(NUT_fc[,i], na.rm = T)
a4 = density(SUT_fc[,i], na.rm = T)
a5 = density(XUT_fc[,i], na.rm = T)
plot(a1, typ = 'l', xlim = c(-5,5), ylim = c(0,.5), lwd = 2,
     main = expression(italic(paste('dhh1', Delta))), xlab = 'Standardized NET-seq FC',
     cex.axis = 2, cex.lab = 2, cex.main = 3)
lines(a2, col = 2, lwd = 2)
lines(a3, col = 3, lwd = 2)
lines(a4, col = 4, lwd = 2)
lines(a5, col = 5, lwd = 2)
legend(x = 'topright', legend = c('Genes', 'CUTs', 'NUTs', 'SUTs', 'XUTs'),
       col = 1:5, lwd = 2, cex = 2)
abline(v = 0, lty = 2, lwd = 2)
dev.off()

pdf(file = 'lsm1_fc_dens.pdf', width = 9)
i = 3
a1 = density(ORF_fc[,i], na.rm = T)
a2 = density(CUT_fc[,i], na.rm = T)
a3 = density(NUT_fc[,i], na.rm = T)
a4 = density(SUT_fc[,i], na.rm = T)
a5 = density(XUT_fc[,i], na.rm = T)
plot(a1, typ = 'l', xlim = c(-5,5), ylim = c(0,.5), lwd = 2,
     main = expression(italic(paste('lsm1', Delta))), xlab = 'Standardized NET-seq FC',
     cex.axis = 2, cex.lab = 2, cex.main = 3)
lines(a2, col = 2, lwd = 2)
lines(a3, col = 3, lwd = 2)
lines(a4, col = 4, lwd = 2)
lines(a5, col = 5, lwd = 2)
legend(x = 'topright', legend = c('Genes', 'CUTs', 'NUTs', 'SUTs', 'XUTs'),
       col = 1:5, lwd = 2, cex = 2)
abline(v = 0, lty = 2, lwd = 2)
dev.off()

pdf(file = 'ccr4_fc_dens.pdf', width = 9)
i = 1
a1 = density(ORF_fc[,i], na.rm = T)
a2 = density(CUT_fc[,i], na.rm = T)
a3 = density(NUT_fc[,i], na.rm = T)
a4 = density(SUT_fc[,i], na.rm = T)
a5 = density(XUT_fc[,i], na.rm = T)
plot(a1, typ = 'l', xlim = c(-5,5), ylim = c(0,.5), lwd = 2,
     main = expression(italic(paste('ccr4', Delta))), xlab = 'Standardized NET-seq FC',
     cex.axis = 2, cex.lab = 2, cex.main = 3)
lines(a2, col = 2, lwd = 2)
lines(a3, col = 3, lwd = 2)
lines(a4, col = 4, lwd = 2)
lines(a5, col = 5, lwd = 2)
legend(x = 'topright', legend = c('Genes', 'CUTs', 'NUTs', 'SUTs', 'XUTs'),
       col = 1:5, lwd = 2, cex = 2)
abline(v = 0, lty = 2, lwd = 2)
dev.off()

pdf(file = 'rpb4_fc_dens.pdf', width = 9)
i = 4
a1 = density(ORF_fc[,i], na.rm = T)
a2 = density(CUT_fc[,i], na.rm = T)
a3 = density(NUT_fc[,i], na.rm = T)
a4 = density(SUT_fc[,i], na.rm = T)
a5 = density(XUT_fc[,i], na.rm = T)
plot(a1, typ = 'l', xlim = c(-5,5), ylim = c(0,.5), lwd = 2,
     main = expression(italic(paste('rpb4', Delta))), xlab = 'Standardized NET-seq FC',
     cex.axis = 2, cex.lab = 2, cex.main = 3)
lines(a2, col = 2, lwd = 2)
lines(a3, col = 3, lwd = 2)
lines(a4, col = 4, lwd = 2)
lines(a5, col = 5, lwd = 2)
legend(x = 'topright', legend = c('Genes', 'CUTs', 'NUTs', 'SUTs', 'XUTs'),
       col = 1:5, lwd = 2, cex = 2)
abline(v = 0, lty = 2, lwd = 2)
dev.off()

######
### Figure 10 - ncRNA metagenes ###
pdf(file = 'CUT_meta.pdf')
par(mar = c(5,5,5,3))
plot(WT_smoothed_CUT$x, WT_smoothed_CUT$y/sum(nf[1:2]), typ = 'l', col = 'gray', lwd = 2,
     xlab = 'Metaposition', ylab = 'Mean NET-seq reads', main = 'CUT metagenes',
     cex.main = 2, cex.axis = 2, cex.lab = 2, ylim = c(0, .14))
lines(xrn1_smoothed_CUT$x, xrn1_smoothed_CUT$y/sum(nf[8:9]), col = 1, lwd = 2)
lines(dhh1_smoothed_CUT$x, dhh1_smoothed_CUT$y/nf[4], col = 2, lwd = 2)
lines(lsm1_smoothed_CUT$x, lsm1_smoothed_CUT$y/nf[5], col = 3, lwd = 2)
lines(ccr4_smoothed_CUT$x, ccr4_smoothed_CUT$y/nf[3], col = 4, lwd = 2)
lines(rpb4_smoothed_CUT$x, rpb4_smoothed_CUT$y/nf[6], col = 6, lwd = 2)
legend(x = 'topright', col = c('gray', 1, 2, 3, 4, 6), cex = 1.5, lwd = 2,
       legend = c('WT',expression(paste(italic("xrn1"), Delta)), expression(paste(italic("dhh1"), Delta)), expression(paste(italic("lsm1"), Delta)), expression(paste(italic("ccr4"), Delta)), expression(paste(italic("rpb4"), Delta))))
dev.off()

pdf(file = 'NUT_meta.pdf')
par(mar = c(5,5,5,3))
plot(WT_smoothed_NUT$x, WT_smoothed_NUT$y/sum(nf[1:2]), typ = 'l', col = 'gray', lwd = 2,
     xlab = 'Metaposition', ylab = 'Mean NET-seq reads', main = 'NUT metagenes',
     cex.main = 2, cex.axis = 2, cex.lab = 2, ylim = c(0, .5))
lines(xrn1_smoothed_NUT$x, xrn1_smoothed_NUT$y/sum(nf[8:9]), col = 1, lwd = 2)
lines(dhh1_smoothed_NUT$x, dhh1_smoothed_NUT$y/nf[4], col = 2, lwd = 2)
lines(lsm1_smoothed_NUT$x, lsm1_smoothed_NUT$y/nf[5], col = 3, lwd = 2)
lines(ccr4_smoothed_NUT$x, ccr4_smoothed_NUT$y/nf[3], col = 4, lwd = 2)
lines(rpb4_smoothed_NUT$x, rpb4_smoothed_NUT$y/nf[6], col = 6, lwd = 2)
legend(x = 'topright', col = c('gray', 1, 2, 3, 4, 6), cex = 1.5, lwd = 2,
       legend = c('WT',expression(paste(italic("xrn1"), Delta)), expression(paste(italic("dhh1"), Delta)), expression(paste(italic("lsm1"), Delta)), expression(paste(italic("ccr4"), Delta)), expression(paste(italic("rpb4"), Delta))))
dev.off()

pdf(file = 'SUT_meta.pdf')
par(mar = c(5,5,5,3))
plot(WT_smoothed_SUT$x, WT_smoothed_SUT$y/sum(nf[1:2]), typ = 'l', col = 'gray', lwd = 2,
     xlab = 'Metaposition', ylab = 'Mean NET-seq reads', main = 'SUT metagenes',
     cex.main = 2, cex.axis = 2, cex.lab = 2, ylim = c(0, .16))
lines(xrn1_smoothed_SUT$x, xrn1_smoothed_SUT$y/sum(nf[8:9]), col = 1, lwd = 2)
lines(dhh1_smoothed_SUT$x, dhh1_smoothed_SUT$y/nf[4], col = 2, lwd = 2)
lines(lsm1_smoothed_SUT$x, lsm1_smoothed_SUT$y/nf[5], col = 3, lwd = 2)
lines(ccr4_smoothed_SUT$x, ccr4_smoothed_SUT$y/nf[3], col = 4, lwd = 2)
lines(rpb4_smoothed_SUT$x, rpb4_smoothed_SUT$y/nf[6], col = 6, lwd = 2)
legend(x = 'topright', col = c('gray', 1, 2, 3, 4, 6), cex = 1.5, lwd = 2,
       legend = c('WT',expression(paste(italic("xrn1"), Delta)), expression(paste(italic("dhh1"), Delta)), expression(paste(italic("lsm1"), Delta)), expression(paste(italic("ccr4"), Delta)), expression(paste(italic("rpb4"), Delta))))
dev.off()

pdf(file = 'XUT_meta.pdf')
par(mar = c(5,5,5,3))
plot(WT_smoothed_XUT$x, WT_smoothed_XUT$y/sum(nf[1:2]), typ = 'l', col = 'gray', lwd = 2,
     xlab = 'Metaposition', ylab = 'Mean NET-seq reads', main = 'XUT metagenes',
     cex.main = 2, cex.axis = 2, cex.lab = 2, ylim = c(0, .12))
lines(xrn1_smoothed_XUT$x, xrn1_smoothed_XUT$y/sum(nf[8:9]), col = 1, lwd = 2)
lines(dhh1_smoothed_XUT$x, dhh1_smoothed_XUT$y/nf[4], col = 2, lwd = 2)
lines(lsm1_smoothed_XUT$x, lsm1_smoothed_XUT$y/nf[5], col = 3, lwd = 2)
lines(ccr4_smoothed_XUT$x, ccr4_smoothed_XUT$y/nf[3], col = 4, lwd = 2)
lines(rpb4_smoothed_XUT$x, rpb4_smoothed_XUT$y/nf[6], col = 6, lwd = 2)
legend(x = 'topright', col = c('gray', 1, 2, 3, 4, 6), cex = 1.5, lwd = 2,
       legend = c('WT',expression(paste(italic("xrn1"), Delta)), expression(paste(italic("dhh1"), Delta)), expression(paste(italic("lsm1"), Delta)), expression(paste(italic("ccr4"), Delta)), expression(paste(italic("rpb4"), Delta))))
dev.off()

pdf(file = 'xrn1_fc_metagenes.pdf', width = 14, height = 9)
par(mfrow = c(2,3))
end_metagene_plotter(gene_ind = 1:4973, gene_title = 'All genes', which_end = 'TSS', which_strain = 'xrn1')
end_metagene_plotter(gene_ind = which(ORF_fc[,5] < -2), gene_title = 'FC < -2', which_end = 'TSS', which_strain = 'xrn1')
end_metagene_plotter(gene_ind = which(ORF_fc[,5] > 0), gene_title = 'FC > 0', which_end = 'TSS', which_strain = 'xrn1')
end_metagene_plotter(gene_ind = 1:4973, gene_title = 'All genes', which_end = 'PAS', which_strain = 'xrn1')
end_metagene_plotter(gene_ind = which(ORF_fc[,5] < -2), gene_title = 'FC < -2', which_end = 'PAS', which_strain = 'xrn1')
end_metagene_plotter(gene_ind = which(ORF_fc[,5] > 0), gene_title = 'FC > 0', which_end = 'PAS', which_strain = 'xrn1')
dev.off()


pdf(file = 'partitioned_metagenes.pdf', width = 12)
par(mfrow = c(3,2))
end_metagene_plotter(gene_ind = which(ORF_fc[,5] < -2), gene_title = 'FC < -2', which_end = 'TSS', which_strain = 'xrn1')
end_metagene_plotter(gene_ind = which(ORF_fc[,5] < -2), gene_title = 'FC < -2', which_end = 'PAS', which_strain = 'xrn1')
a = which(ORF_fc[,5] > -2 & ORF_fc[,5] < 0)
end_metagene_plotter(gene_ind = a, gene_title = '-2 < FC < 0', which_end = 'TSS', which_strain = 'xrn1')
end_metagene_plotter(gene_ind = a, gene_title = '-2 < FC < 0', which_end = 'PAS', which_strain = 'xrn1')
end_metagene_plotter(gene_ind = which(ORF_fc[,5] > 0), gene_title = 'FC > 0', which_end = 'TSS', which_strain = 'xrn1')
end_metagene_plotter(gene_ind = which(ORF_fc[,5] > 0), gene_title = 'FC > 0', which_end = 'PAS', which_strain = 'xrn1')

end_metagene_plotter(gene_ind = which(ORF_fc[,2] < -2), gene_title = 'FC < -2', which_end = 'TSS', which_strain = 'dhh1')
end_metagene_plotter(gene_ind = which(ORF_fc[,2] < -2), gene_title = 'FC < -2', which_end = 'PAS', which_strain = 'dhh1')
a = which(ORF_fc[,2] > -2 & ORF_fc[,2] < 0)
end_metagene_plotter(gene_ind = a, gene_title = '-2 < FC < 0', which_end = 'TSS', which_strain = 'dhh1')
end_metagene_plotter(gene_ind = a, gene_title = '-2 < FC < 0', which_end = 'PAS', which_strain = 'dhh1')
end_metagene_plotter(gene_ind = which(ORF_fc[,2] > 0), gene_title = 'FC > 0', which_end = 'TSS', which_strain = 'dhh1')
end_metagene_plotter(gene_ind = which(ORF_fc[,2] > 0), gene_title = 'FC > 0', which_end = 'PAS', which_strain = 'dhh1')

end_metagene_plotter(gene_ind = which(ORF_fc[,3] < -2), gene_title = 'FC < -2', which_end = 'TSS', which_strain = 'lsm1')
end_metagene_plotter(gene_ind = which(ORF_fc[,3] < -2), gene_title = 'FC < -2', which_end = 'PAS', which_strain = 'lsm1')
a = which(ORF_fc[,3] > -2 & ORF_fc[,3] < 0)
end_metagene_plotter(gene_ind = a, gene_title = '-2 < FC < 0', which_end = 'TSS', which_strain = 'lsm1')
end_metagene_plotter(gene_ind = a, gene_title = '-2 < FC < 0', which_end = 'PAS', which_strain = 'lsm1')
end_metagene_plotter(gene_ind = which(ORF_fc[,3] > 0), gene_title = 'FC > 0', which_end = 'TSS', which_strain = 'lsm1')
end_metagene_plotter(gene_ind = which(ORF_fc[,3] > 0), gene_title = 'FC > 0', which_end = 'PAS', which_strain = 'lsm1')

end_metagene_plotter(gene_ind = which(ORF_fc[,1] < -2), gene_title = 'FC < -2', which_end = 'TSS', which_strain = 'ccr4')
end_metagene_plotter(gene_ind = which(ORF_fc[,1] < -2), gene_title = 'FC < -2', which_end = 'PAS', which_strain = 'ccr4')
a = which(ORF_fc[,1] > -2 & ORF_fc[,1] < 0)
end_metagene_plotter(gene_ind = a, gene_title = '-2 < FC < 0', which_end = 'TSS', which_strain = 'ccr4')
end_metagene_plotter(gene_ind = a, gene_title = '-2 < FC < 0', which_end = 'PAS', which_strain = 'ccr4')
end_metagene_plotter(gene_ind = which(ORF_fc[,1] > 0), gene_title = 'FC > 0', which_end = 'TSS', which_strain = 'ccr4')
end_metagene_plotter(gene_ind = which(ORF_fc[,1] > 0), gene_title = 'FC > 0', which_end = 'PAS', which_strain = 'ccr4')

end_metagene_plotter(gene_ind = which(ORF_fc[,4] < -2), gene_title = 'FC < -2', which_end = 'TSS', which_strain = 'rpb4')
end_metagene_plotter(gene_ind = which(ORF_fc[,4] < -2), gene_title = 'FC < -2', which_end = 'PAS', which_strain = 'rpb4')
a = which(ORF_fc[,4] > -2 & ORF_fc[,4] < 0)
end_metagene_plotter(gene_ind = a, gene_title = '-2 < FC < 0', which_end = 'TSS', which_strain = 'rpb4')
end_metagene_plotter(gene_ind = a, gene_title = '-2 < FC < 0', which_end = 'PAS', which_strain = 'rpb4')
end_metagene_plotter(gene_ind = which(ORF_fc[,4] > 0), gene_title = 'FC > 0', which_end = 'TSS', which_strain = 'rpb4')
end_metagene_plotter(gene_ind = which(ORF_fc[,4] > 0), gene_title = 'FC > 0', which_end = 'PAS', which_strain = 'rpb4')
dev.off()

a = WT_TSS + WT_pilot_TSS
a2 = WT_TES + WT_pilot_TES
b = XRN1_TSS + XRN1_TSS
b2 = XRN1_TES + XRN1_TES
#b = CCR4_TSS
#b2 = CCR4_TES

i1 = which(ORF_fc[,5] <= -2)
i2 = which(ORF_fc[,5] < 0 & ORF_fc[,5] > -2)
i3 = which(ORF_fc[,5] >= 0)

z = rep(0, 4973)
z[i1] = 1
z[i2] = 2
z[i3] = 3

plot(log2(colSums(b)/colSums(a)), log2(colSums(b2)/colSums(a2)),
     col = z, xlab = 'Raw Log2 TSS FC', ylab = 'Raw Log2 PAS FC',
     main = expression(paste('Comparison of TSS and PAS FC, ', italic("xrn1"), Delta)))
legend(x = 'bottomright', legend = c('Gene FC < -2', '-2 < Gene FC < 0', '0 < Gene FC'), pch = 1, col = 1:3)
cor(log2(colSums(b)/colSums(a)), log2(colSums(b2)/colSums(a2)), method = 'spearman', use = 'complete.obs')
cor(log2(colSums(b)/colSums(a))[i2], log2(colSums(b2)/colSums(a2))[i2], method = 'spearman', use = 'complete.obs')


hist(ORF_fc$XRN1, breaks = seq(-10, 10, by = .5), xlab = 'Standardized Log2 NET-seq FC', main = expression(paste(italic(xrn1), Delta)))
hist(ORF_fc$XRN1[ORF_fc$XRN1<=-2], breaks = seq(-10, 10, by = .5), col = 3, add = T)
hist(ORF_fc$XRN1[ORF_fc$XRN1>=0], breaks = seq(-10, 10, by = .5), col = 2, add = T)
legend(x = 'topright', fill = c(3,2), legend = c('Stimulated by Xrn1', 'Repressed by Xrn1'))

pdf('~/Documents/R_Projects/NETSEQ/paper/DF2/Figures/main/FC_SAGA_TFIID.pdf')
hist(ORF_fc$XRN1[SAGA_ind], breaks = 60, freq = F, ylim = c(0, .4), col = rgb(0,0,1,.5),
     main = 'FCs by association with SAGA or TFIID', xlab = 'Standardized NET-seq Pol II FCs')
hist(ORF_fc$XRN1[TFIID_ind], breaks = 75, add = T, freq = F, col = rgb(1,0,0,.5))
legend(x = 'topright', legend = c('SAGA-dominated', 'TFIID-dominated'), fill = c(rgb(0,0,1,.5), rgb(1,0,0,.5)))
abline(v = median(ORF_fc$XRN1[SAGA_ind], na.rm = T), col = 'blue', lwd = 2, lty = 2)
abline(v = median(ORF_fc$XRN1[TFIID_ind], na.rm = T), col = 'red', lwd = 2, lty = 2)
legend(x = 'topleft', lty = 2, lwd = 2, legend = 'median value')
dev.off()

break_vec = seq(-4, 4, by = .1)
pdf('~/Documents/R_Projects/NETSEQ/paper/DF2/Figures/review/TSS_FC_SAGA_TFIID.pdf')
par(mar = c(5,5,5,3))
hist(z[SAGA_ind], breaks = break_vec, freq = F, ylim = c(0, 1), col = rgb(0,0,1,.5),
     main = 'FCs near TSS', xlab = 'Log2 NET-seq Pol II FCs', cex.main = 2, 
     cex.lab = 2, cex.axis = 2)
hist(z[TFIID_ind], breaks = break_vec, add = T, freq = F, col = rgb(1,0,0,.5))
legend(x = 'topright', legend = c('SAGA-dominated', 'TFIID-dominated'), fill = c(rgb(0,0,1,.5), rgb(1,0,0,.5)), cex = 1.35)
abline(v = median(z[SAGA_ind], na.rm = T), col = 'blue', lwd = 2, lty = 2)
abline(v = median(z[TFIID_ind], na.rm = T), col = 'red', lwd = 2, lty = 2)
legend(x = 'topleft', lty = 2, lwd = 2, legend = 'median value', cex = 1.25)
dev.off()

pdf('~/Documents/R_Projects/NETSEQ/paper/DF2/Figures/review/PAS_FC_SAGA_TFIID.pdf')
par(mar = c(5,5,5,3))
hist(y[SAGA_ind], breaks = break_vec, freq = F, ylim = c(0, 1), col = rgb(0,0,1,.5),
     main = 'FCs near PAS', xlab = 'Log2 NET-seq Pol II FCs', cex.main = 2, 
     cex.lab = 2, cex.axis = 2)
hist(y[TFIID_ind], breaks = break_vec, add = T, freq = F, col = rgb(1,0,0,.5))
legend(x = 'topright', legend = c('SAGA-dominated', 'TFIID-dominated'), fill = c(rgb(0,0,1,.5), rgb(1,0,0,.5)), cex = 1.35)
abline(v = median(y[SAGA_ind], na.rm = T), col = 'blue', lwd = 2, lty = 2)
abline(v = median(y[TFIID_ind], na.rm = T), col = 'red', lwd = 2, lty = 2)
legend(x = 'topleft', lty = 2, lwd = 2, legend = 'median value', cex = 1.25)
dev.off()

break_vec2 = seq(from = -8, to = 8, by = .1)
pdf('~/Documents/R_Projects/NETSEQ/paper/DF2/Figures/review/DF_binding_SAGA_TFIID.pdf')
par(mar = c(5,5,5,3))
hist(zz[SAGA_ind], breaks = break_vec2, freq = F, ylim = c(0, 1), col = rgb(0,0,1,.5),
     main = 'Relative binding of DFs to mRNAs', xlab = 'Log2 XRN1/SKI2 CRAC ratio', cex.main = 2, 
     cex.lab = 2, cex.axis = 2, xlim = c(-5, 6))
hist(zz[TFIID_ind], breaks = break_vec2, add = T, freq = F, col = rgb(1,0,0,.5))
legend(x = 'topright', legend = c('SAGA-dominated', 'TFIID-dominated'), fill = c(rgb(0,0,1,.5), rgb(1,0,0,.5)), cex = 1.35)
abline(v = median(zz[SAGA_ind], na.rm = T), col = 'blue', lwd = 2, lty = 2)
abline(v = median(zz[TFIID_ind], na.rm = T), col = 'red', lwd = 2, lty = 2)
legend(x = 'topleft', lty = 2, lwd = 2, legend = 'median value', cex = 1.25)
dev.off()

zz = rowMeans(log2(XRN1_CRAC_mat_norm+3))-rowMeans(log2(SKI2_CRAC_mat_norm+3))


plot(-100:500, rowMeans(log2((BG_WT_TSS_2[-c(1:15,617:631),b]+s)/(WT_NC_TSS[,a]+s))), typ = 'l', main = 'WT', xlab = 'Distance to TSS', ylab = 'Log2 BioGRO/NET-seq', ylim = c(0.5, 1.45))
plot(-100:500, rowMeans(log2((BG_xrn1_TSS_2[-c(1:15,617:631),b]+s)/(XRN1_TSS[,a]+XRN1_pilot_TSS[,a]+s))), typ = 'l', main = expression(paste(italic('xrn1'), Delta)), xlab = 'Distance to TSS', ylab = 'Log2 BioGRO/NET-seq', ylim = c(0.5, 1.45))

plot(-150:150, rowMeans(log2((BG_WT_TES_2[-c(1:15,317:331),b]+s)/(WT_NC_TES[,a]+s))), typ = 'l', main = 'WT', xlab = 'Distance to PAS', ylab = 'Log2 BioGRO/NET-seq', ylim = c(0.85, 1.4))
plot(-150:150, rowMeans(log2((BG_xrn1_TES_2[-c(1:15,317:331),b]+s)/(XRN1_TES[,a]+XRN1_pilot_TES[,a]+s))), typ = 'l', main = expression(paste(italic('xrn1'), Delta)), xlab = 'Distance to PAS', ylab = 'Log2 BioGRO/NET-seq', ylim = c(0.85, 1.4))


