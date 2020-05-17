nn_sd_est_2 = function(exp_table, n_neighbor = nrow(exp_table/100), cutoff = 1, smoother = 2/3){
  l = nrow(exp_table)
  exp_sorted = arrange(exp_table, A)
  k = ceiling(n_neighbor)
  below_cutoff = which(exp_sorted$A < cutoff)
  if(length(below_cutoff > 0)){
    exp_nz = exp_sorted[-below_cutoff,]
  } else {
    exp_nz = exp_sorted
  }
  k_nn = lapply(1:nrow(exp_nz), function(x) (x-k+(x<k)*(k-x)+order(abs(exp_nz[x,1] - exp_nz[max(1,x-k+1):min(nrow(exp_nz),x+k-1),1])))[1:k])
  sd_est = unlist(lapply(k_nn, function(x) mad(exp_nz$M[x], center = 0)))
  sd_df = data.frame(cbind(sd_est, exp_nz$A))
  colnames(sd_df) = c('SD', 'A')
  #plot(sd_df$A, sd_df$SD, xlab = 'Avg Log2 Expression', ylab = 'SD in Neighborhood', main = 'Gene kNN SD Estimates, with smoothing')
  smoothed_sd = data.frame(lowess(unique(sd_df)$A, unique(sd_df)$SD, f = smoother))
  colnames(smoothed_sd) = c('A', 'SD')
  to_return = list(smoothed_sd, unique(sd_df))
  return(to_return)
}

differential_expression = function(expression_table, num_neighbors = 25, min_reads = 5, sd_cut = 1, sm_val = 1/4, norm_factors = NULL, WT_ind = c(1,2), xrn1_ind = c(8,9)){
  
  if(is.null(norm_factors)){
    norm_factors = estimateSizeFactorsForMatrix(expression_table)
  }
  normalized_table = sweep(expression_table, 2, norm_factors, '/')
  norm_log2_table = log2(normalized_table+1)
  
  M.WT = norm_log2_table[,WT_ind[1]]-norm_log2_table[,WT_ind[2]]
  A.WT = .5*(norm_log2_table[,WT_ind[1]]+norm_log2_table[,WT_ind[2]])
  MA.WT = cbind.data.frame(A.WT,M.WT)
  colnames(MA.WT) = c('A', 'M')
  
  M.xrn = norm_log2_table[,xrn1_ind[1]]-norm_log2_table[,xrn1_ind[2]]
  A.xrn = .5*(norm_log2_table[,xrn1_ind[1]]+norm_log2_table[,xrn1_ind[2]])
  MA.xrn = cbind.data.frame(A.xrn,M.xrn)
  colnames(MA.xrn) = c('A', 'M')
  
  MA.total = rbind.data.frame(MA.WT, MA.xrn)
  
  WT_too_low = which(rowSums(expression_table[,WT_ind] < min_reads) > 0)
  xrn_enough_reads = (expression_table[,xrn1_ind] >= min_reads)
  xrn_too_low = which(rowSums(xrn_enough_reads) < 2)
  low_counts = cbind(expression_table[,-c(WT_ind, xrn1_ind)] < min_reads, rowSums(xrn_enough_reads) == 0)
  
  rep_too_low = union(WT_too_low, xrn_too_low)
  MA_to_est_SD = MA.total[-rep_too_low,]
  
  #plot(MA.total$A, MA.total$M, xlab = 'Avg Log2 Expression (A)', ylab = 'Log2 Fold Change (M)', main = 'Normalized Replicate MA Plot')
  ggplot(data = MA.total)+geom_point(aes(x = A, y = M, alpha = .1)) + 
    ggtitle('Normalized Replicate MA Plot') + xlab('Avg Log2 Expression') + 
    ylab('Log2 Fold Change')+ theme(legend.position="none")
  
  sd_est = nn_sd_est_2(MA_to_est_SD, n_neighbor = num_neighbors, cutoff = sd_cut, smoother = sm_val)
  sd_smoothed = sd_est[[1]]
  sd_NN = sd_est[[2]]
  sd_scatter = ggplot()+geom_point(data = sd_NN, aes(x = A, y = SD))+
    xlab('Avg Log 2 Expression (A)')+ylab('Estimated SD')+ggtitle('Smoothed SD Estimate, Genes')+
    geom_line(data = sd_smoothed, aes(x = A, y = SD, colour = 'red'))+ theme(legend.position="none")
  #lines(sd_est$x, sd_est$y, type = 'l', col = 'red')
  #plot(sd_est$x, sd_est$y, type = 'l', xlab = 'Avg Log2 Expression', ylab = 'Est SD', 
  #main = 'Smoothed SD Estimates, Genes')
  
  sd_df = data.frame(sd_smoothed)
  colnames(sd_df) = c('A', 'SD')
  
  # Separates mutants and wild types
  mutant_table = norm_log2_table
  #mutant_table[,xrn1_ind[1]] = .5*rowSums(mutant_table[,xrn1_ind])
  mutant_table[,xrn1_ind[1]] = rowSums(mutant_table[,xrn1_ind]*xrn_enough_reads) / rowSums(xrn_enough_reads)
  mutant_table = mutant_table[,-c(WT_ind, xrn1_ind[2])]
  WT_table = norm_log2_table[,WT_ind]
  
  # Subtracts wild type values (equivalent to log fold change)
  mutant_table_1 <- mutant_table - WT_table[,1]
  mutant_table_2 <- mutant_table - WT_table[,2]
  
  gene_sd = sapply(MA.WT$A, function(x) unique(sd_df$SD[which.min(abs(sd_df$A - x))]))
  gene_sd_2 = sapply(MA.xrn$A, function(x) unique(sd_df$SD[which.min(abs(sd_df$A - x))]))
  
  rep_vec = c(MA.WT[,2], MA.xrn[,2]) / c(gene_sd, gene_sd_2)
  null_sd = mad(rep_vec)
  
  # Creates tables for standardized values
  mutant_table_stand_1 <- sweep(mutant_table_1 / gene_sd, 2, c(rep(1, 5), sqrt(4/3)), '*')
  mutant_table_stand_2 <- sweep(mutant_table_2 / gene_sd, 2, c(rep(1, 5), sqrt(4/3)), '*')
  
  mutant_table_stand_1[rowSums(xrn_enough_reads)==1,6] = mutant_table_stand_1[rowSums(xrn_enough_reads)==1,6]*sqrt(3/4)
  mutant_table_stand_2[rowSums(xrn_enough_reads)==1,6] = mutant_table_stand_2[rowSums(xrn_enough_reads)==1,6]*sqrt(3/4)
  
  p1 = matrix(2*pnorm(-abs(unlist(mutant_table_stand_1)), 0, null_sd), ncol = 6)
  p2 = matrix(2*pnorm(-abs(unlist(mutant_table_stand_2)), 0, null_sd), ncol = 6)
  
  q1 = p1
  q2 = p2
  
  for(i in 1:ncol(p1)){
    q1[,i] <- p.adjust(p1[,i], "fdr")
    q2[,i] <- p.adjust(p2[,i], "fdr")
  }
  
  max_q_values <- (q1>q2) * q1 + (q2>=q1) * q2
  
  WT_avg = .5*rowSums(WT_table)
  
  fold_changes = mutant_table - WT_avg
  SFC = sweep(fold_changes / gene_sd, 2, c(rep(sqrt(4/3), 5), sqrt(2)), '*')
  SFC[rowSums(xrn_enough_reads)==1,6] = SFC[rowSums(xrn_enough_reads)==1,6]*sqrt(2/3)
  P_test = matrix(2*pnorm(-abs(unlist(SFC)), 0, null_sd), ncol = 6)
  
  P_test[low_counts] = NA
  P_test[WT_too_low, ] = rep(NA, 6)
  Q_test = P_test
  
  for(i in 1:ncol(Q_test)){
    Q_test[,i] <- p.adjust(Q_test[,i], "fdr")
  }
  
  p = .05
  q.cut = .15
  
  DE_up = lapply(1:6, function(x) which(Q_test[,x] < p & SFC[,x] > 0 & !is.na(Q_test[,x]) & max_q_values[,x] < q.cut))
  DE_down = lapply(1:6, function(x) which(Q_test[,x] < p & SFC[,x] < 0 & !is.na(Q_test[,x]) & max_q_values[,x] < q.cut))
  DE_list = list(DE_up, DE_down)
  
  to_return = list(DE_list, Q_test, max_q_values, sd_df, normalized_table, P_test, SFC, list(p1,p2))
  names(to_return) = c('DE_list', 'Q_test', 'max_q_values', 'sd_df', 'normalized_table', 'P_test', 'SFC', 'p_list')
  
  return(to_return)
  
}

true_false_positives = function(set1, set2, true_set, list_or_vec = 'vec'){
  t_all = intersect(true_set, intersect(set1, set2))
  t_neither = setdiff(true_set, union(set1, set2))
  t_us_only = setdiff(intersect(set1, true_set), set2)
  t_them_only = setdiff(intersect(set2, true_set), set1)
  
  f_all = setdiff(intersect(set1, set2), true_set)
  f_neither = integer(0)
  f_us_only = setdiff(setdiff(set1, true_set), set2)
  f_them_only = setdiff(setdiff(set2, true_set), set1)
  
  my_list = list(t_all, t_neither, t_us_only, t_them_only, f_all, f_neither, f_us_only, f_them_only)
  
  my_vec = unlist(lapply(my_list, function(x) length(x)))
  
  if(list_or_vec == 'vec'){
    return(my_vec)
  } else {
    return(my_list)
  }
  
}

netseq_plotter = function(gene_name, plus_list = WT_NC_plus, minus_list = WT_NC_minus, annot = ORF_annot, sense_flag = 0, plot_flag = 1){
  # gene_name: gene common name, locus name, or row index
  # plus_list and minus_list: lists containing number of (nonzero) reads at sites along the plus and minus strands
  # annot: annotation of gene names and locations
  # sense_val: 1 for sense, -1 for antisense, 0 for both
  # plot_flag: should read profile be plotted
  
  # Finds requested gene in annotation file
  if(is.numeric(gene_name)){
    i = gene_name
  } else if(is.character(gene_name)){
    i = which(annot[,6] == gene_name)
    if(length(i) == 0){
      i = which(annot[,5] == gene_name)
    }
  }
  
  # Extracts chromosome and strand information for gene
  gene_chr = as.numeric(annot[i,1])
  gene_str = annot[i,2]
  
  # Get reads for sense and antisense strands
  if(gene_str == '+'){
    s_dummy = matrix(plus_list[[gene_chr]][which(plus_list[[gene_chr]][,1] %in% annot[i,3]:annot[i,4]),], ncol = 2)
    as_dummy = matrix(minus_list[[gene_chr]][which(minus_list[[gene_chr]][,1] %in% annot[i,3]:annot[i,4]),], ncol = 2)
  } else if(gene_str == '-'){
    s_dummy = matrix(minus_list[[gene_chr]][which(minus_list[[gene_chr]][,1] %in% annot[i,3]:annot[i,4]),], ncol = 2)
    as_dummy = matrix(plus_list[[gene_chr]][which(plus_list[[gene_chr]][,1] %in% annot[i,3]:annot[i,4]),], ncol = 2)
  }
  sense_reads = cbind(annot[i,3]:annot[i,4], rep(0, length(annot[i,3]:annot[i,4])))
  sense_reads[which(sense_reads[,1] %in% s_dummy[,1]), 2] = s_dummy[,2]
  
  antisense_reads = cbind(annot[i,3]:annot[i,4], rep(0, length(annot[i,3]:annot[i,4])))
  antisense_reads[which(antisense_reads[,1] %in% as_dummy[,1]), 2] = as_dummy[,2]
  
  if(plot_flag){
    if(sense_flag == 1){
      if(gene_str == '+'){
        plot(1:nrow(sense_reads), sense_reads[,2], typ = 'l', main = annot[i,6], xlab = 'Position', ylab = 'Counts')
      } else if(gene_str == '-'){
        plot(1:nrow(sense_reads), rev(sense_reads[,2]), typ = 'l', main = annot[i,6], xlab = 'Position', ylab = 'Counts')
      }
    } else if(sense_flag == -1){
      if(gene_str == '+'){
        plot(1:nrow(antisense_reads), antisense_reads[,2], typ = 'l', main = annot[i,6], xlab = 'Position', ylab = 'Antisense Counts')
      } else if(gene_str == '-'){
        plot(1:nrow(antisense_reads), rev(antisense_reads[,2]), typ = 'l', main = annot[i,6], xlab = 'Position', ylab = 'Antisense Counts')
      }
    } else if(sense_flag == 0){
      if(gene_str == '+'){
        plot(annot[i,3]:annot[i,4], sense_reads[,2], typ = 'l', main = annot[i,6], xlab = "Position (5' -> 3')", ylab = 'Counts', ylim = range(c(sense_reads[,2], -antisense_reads[,2])))
        lines(annot[i,3]:annot[i,4], -antisense_reads[,2], col = 'red')
      } else if(gene_str == '-'){
        plot(annot[i,3]:annot[i,4], -sense_reads[,2], typ = 'l', main = annot[i,6], xlab = "Position (3' <- 5')", ylab = 'Counts', ylim = range(c(-sense_reads[,2], antisense_reads[,2])))
        lines(annot[i,3]:annot[i,4], antisense_reads[,2], col = 'red')
      }
    }
    
  }
  
  return(c(sum(sense_reads[,2]), sum(antisense_reads[,2])))
  
}

five_vs_three_reads = function(gene_ind = 1, annot = ORF_annot, p = .15, plus_list = WT_NC_plus, minus_list = WT_NC_minus){
  chr = ORF_annot$chr[gene_ind]
  strand = ORF_annot$strand[gene_ind]
  l = abs(ORF_annot$start[gene_ind] - ORF_annot$end[gene_ind])
  if(strand == '+'){
    five_start = ORF_annot$start[gene_ind]
    five_end = ceiling(five_start+p*l)
    three_end = ORF_annot$end[gene_ind]
    three_start = floor(three_end-p*l)
  } else if(strand == '-'){
    five_start = ORF_annot$end[gene_ind]
    five_end = floor(five_start - p*l)
    three_end = ORF_annot$start[gene_ind]
    three_start = ceiling(three_end+p*l)
  }
  
  five_reads = interval_reads(start = five_start, end = five_end, chrom = chr, strand = strand, 
                              plus_list = plus_list, minus_list = minus_list, plot_flag = F)[1]
  three_reads = interval_reads(start = three_start, end = three_end, chrom = chr, strand = strand, 
                               plus_list = plus_list, minus_list = minus_list, plot_flag = F)[1]
  
  read_vec = c(five_reads, three_reads)
  names(read_vec) = c("5'", "3'")
  
  return(read_vec)
}

interval_reads = function(start, end, chrom, strand = '+', plus_list = WT_NC_plus, minus_list = WT_NC_minus, sense_flag = 0, plot_flag = 1, sums = T){
  # start: starting genomic position
  # end: ending genomic position
  # chrom: chromosome of interest
  # strand: strand of interest if not both
  # plus_list and minus_list: lists containing number of (nonzero) reads at sites along the plus and minus strands
  # sense_val: 1 for sense, -1 for antisense, 0 for both
  # plot_flag: should read profile be plotted
  
  if(chrom == 'M'){
    chrom = 17
  }
  chrom = as.numeric(as.character(chrom))
  # Get reads for sense and antisense strands
  if(strand == '+'){
    site_range = min(start, end):max(start, end)
    s_dummy = matrix(plus_list[[chrom]][which(plus_list[[chrom]][,1] %in% site_range),], ncol = 2)
    as_dummy = matrix(minus_list[[chrom]][which(minus_list[[chrom]][,1] %in% site_range),], ncol = 2)
  } else if(strand == '-'){
    site_range = min(start, end):max(start, end)
    s_dummy = matrix(minus_list[[chrom]][which(minus_list[[chrom]][,1] %in% site_range),], ncol = 2)
    as_dummy = matrix(plus_list[[chrom]][which(plus_list[[chrom]][,1] %in% site_range),], ncol = 2)
  }
  sense_reads = cbind(site_range, rep(0, length(site_range)))
  sense_reads[which(sense_reads[,1] %in% s_dummy[,1]), 2] = s_dummy[,2]
  
  antisense_reads = cbind(site_range, rep(0, length(site_range)))
  antisense_reads[which(antisense_reads[,1] %in% as_dummy[,1]), 2] = as_dummy[,2]
  
  if(plot_flag){
    if(sense_flag == 1){
      if(strand == '+'){
        plot(1:nrow(sense_reads), sense_reads[,2], typ = 'l', main = paste('Chromosome', chrom), xlab = 'Position', ylab = 'Counts')
      } else if(strand == '-'){
        plot(1:nrow(sense_reads), rev(sense_reads[,2]), typ = 'l', main = paste('Chromosome', chrom), xlab = 'Position', ylab = 'Counts')
      }
    } else if(sense_flag == -1){
      if(strand == '+'){
        plot(1:nrow(antisense_reads), antisense_reads[,2], typ = 'l', main = paste('Chromosome', chrom), xlab = 'Position', ylab = 'Antisense Counts')
      } else if(strand == '-'){
        plot(1:nrow(antisense_reads), rev(antisense_reads[,2]), typ = 'l', main = paste('Chromosome', chrom), xlab = 'Position', ylab = 'Antisense Counts')
      }
    } else if(sense_flag == 0){
      if(strand == '+'){
        plot(site_range, sense_reads[,2], typ = 'l', main = paste('Chromosome', chrom), xlab = "Position", ylab = 'Counts', ylim = range(c(sense_reads[,2], -antisense_reads[,2])))
        lines(site_range, -antisense_reads[,2], col = 'red')
      } else if(strand == '-'){
        plot(site_range, -sense_reads[,2], typ = 'l', main = paste('Chromosome', chrom), xlab = "Position", ylab = 'Counts', ylim = range(c(-sense_reads[,2], antisense_reads[,2])))
        lines(site_range, antisense_reads[,2], col = 'red')
      }
    }
    
  }
  
  if(sums){
    to_return = c(sum(sense_reads[,2]), sum(antisense_reads[,2]))
  }
  else {
    to_return = cbind(sense_reads[,2], antisense_reads[,2])
  }
  return(to_return)
  
}

complex_enrichment <- function(de_list, IC_annot, genes_table, complex_id = 'SAGA-dominated'){
  # Function to compute hypergeometric enrichment scores for initiation complexes
  
  de_table = genes_table[genes_table$commonName %in% de_list,]
  dt_hit = sum(de_table$complex == complex_id)
  hit = sum(genes_table$complex == complex_id)
  total = nrow(genes_table)
  dt_total = nrow(de_table)
  #f = matrix(c(dt_hit, dt_total-dt_hit, hit-dt_hit, total-dt_total-(hit-dt_hit)), ncol = 2)
  #fisher.test(f, alternative = 'greater')
  
  # Compute the upper tail p-value using hypergeometric test
  hyper_value = phyper(dt_hit-1, hit, total-hit, dt_total, lower.tail = F)
  
  to_be_returned <- data.frame(hyper_value, dt_hit, dt_total, dt_hit/dt_total, hit, total, hit/total)
  
  colnames(to_be_returned) <- c('P-value', 'DT genes w/ complex', 
                                'DT genes', 'Fraction in DT', 'Genes w/ complex', 'All genes', 'All fraction')
  
  return(to_be_returned)
  
}

ncRNA_function = function(my_table, norm_factors = SizeFactors, ncRNA = 'XUT', data_flag = 1){
  
  norm_table = sweep(my_table, 2, norm_factors, '/')
  log_norm_table = log2(norm_table+1)
  
  if(data_flag == 1){
    suff_reads = which(norm_table[,1]+norm_table[,2] >= 40)
    
    DT = differential_expression(expression_table = my_table, num_neighbors = 25, min_reads = 5, sd_cut = 3, sm_val = 1/4, WT_ind = c(1,2), xrn1_ind = c(8,9))
    
    if(ncRNA == 'XUT'){
      my_lim = c(-.7,.2)
    } else if(ncRNA == 'SUT'){
      my_lim = c(-.8, .2)
    } else if(ncRNA == 'CUT'){
      my_lim = c(-1.1, .4)
    } else if(ncRNA == 'NUT'){
      my_lim = c(-.8,.1)
    }
    
    barplot(rbind(colMedians( (log_norm_table - log_norm_table[,1])[suff_reads,-c(1,2)]),
                  colMedians( (log_norm_table - log_norm_table[,2])[suff_reads,-c(1,2)])), beside = T,
            col = c('white', 'gray'), legend = c('vs WT', 'vs WT pilot'), ylab = 'Median Log2 FC', main = paste(ncRNA, 'RNAPII levels'), args.legend = list(x = 22, y = .27, cex = .75),
            xlab = 'Sample',
            ylim = my_lim, names = c(expression(paste('Ccr4', Delta)), expression(paste('Dhh1', Delta)), expression(paste('Lsm1', Delta)), expression(paste('Rpb4', Delta)), expression(paste('Sfp1', Delta)), expression(paste('Xrn1', Delta)), expression(paste('Xrn1', Delta, 'pilot'))))
    
    abline(h = 0)
    ###############################################################################################
    w = wilcox.test(log_norm_table[suff_reads,3], log_norm_table[suff_reads,1], conf.int = T, paired = T)
    arrows(1.5, w$conf.int[1], 1.5, w$conf.int[2], angle=90, code=3, length=.05, col = 'black', lwd = 1)
    
    w = wilcox.test(log_norm_table[suff_reads,4], log_norm_table[suff_reads,1], conf.int = T, paired = T)
    arrows(4.5, w$conf.int[1], 4.5, w$conf.int[2], angle=90, code=3, length=.05, col = 'black', lwd = 1)
    
    w = wilcox.test(log_norm_table[suff_reads,5], log_norm_table[suff_reads,1], conf.int = T, paired = T)
    arrows(7.5, w$conf.int[1], 7.5, w$conf.int[2], angle=90, code=3, length=.05, col = 'black', lwd = 1)
    
    w = wilcox.test(log_norm_table[suff_reads,6], log_norm_table[suff_reads,1], conf.int = T, paired = T)
    arrows(10.5, w$conf.int[1], 10.5, w$conf.int[2], angle=90, code=3, length=.05, col = 'black', lwd = 1)
    
    w = wilcox.test(log_norm_table[suff_reads,7], log_norm_table[suff_reads,1], conf.int = T, paired = T)
    arrows(13.5, w$conf.int[1], 13.5, w$conf.int[2], angle=90, code=3, length=.05, col = 'black', lwd = 1)
    
    w = wilcox.test(log_norm_table[suff_reads,8], log_norm_table[suff_reads,1], conf.int = T, paired = T)
    arrows(16.5, w$conf.int[1], 16.5, w$conf.int[2], angle=90, code=3, length=.05, col = 'black', lwd = 1)
    
    w = wilcox.test(log_norm_table[suff_reads,9], log_norm_table[suff_reads,1], conf.int = T, paired = T)
    arrows(19.5, w$conf.int[1], 19.5, w$conf.int[2], angle=90, code=3, length=.05, col = 'black', lwd = 1)
    
    
    w = wilcox.test(log_norm_table[suff_reads,3], log_norm_table[suff_reads,2], conf.int = T, paired = T)
    arrows(2.5, w$conf.int[1], 2.5, w$conf.int[2], angle=90, code=3, length=.05, col = 'black', lwd = 1)
    
    w = wilcox.test(log_norm_table[suff_reads,4], log_norm_table[suff_reads,2], conf.int = T, paired = T)
    arrows(5.5, w$conf.int[1], 5.5, w$conf.int[2], angle=90, code=3, length=.05, col = 'black', lwd = 1)
    
    w = wilcox.test(log_norm_table[suff_reads,5], log_norm_table[suff_reads,2], conf.int = T, paired = T)
    arrows(8.5, w$conf.int[1], 8.5, w$conf.int[2], angle=90, code=3, length=.05, col = 'black', lwd = 1)
    
    w = wilcox.test(log_norm_table[suff_reads,6], log_norm_table[suff_reads,2], conf.int = T, paired = T)
    arrows(11.5, w$conf.int[1], 11.5, w$conf.int[2], angle=90, code=3, length=.05, col = 'black', lwd = 1)
    
    w = wilcox.test(log_norm_table[suff_reads,7], log_norm_table[suff_reads,2], conf.int = T, paired = T)
    arrows(14.5, w$conf.int[1], 14.5, w$conf.int[2], angle=90, code=3, length=.05, col = 'black', lwd = 1)
    
    w = wilcox.test(log_norm_table[suff_reads,8], log_norm_table[suff_reads,2], conf.int = T, paired = T)
    arrows(17.5, w$conf.int[1], 17.5, w$conf.int[2], angle=90, code=3, length=.05, col = 'black', lwd = 1)
    
    w = wilcox.test(log_norm_table[suff_reads,9], log_norm_table[suff_reads,2], conf.int = T, paired = T)
    arrows(20.5, w$conf.int[1], 20.5, w$conf.int[2], angle=90, code=3, length=.05, col = 'black', lwd = 1)
    ###############################################################################################
    
    log_norm_avg_table = cbind(rowMeans(log_norm_table[,1:2]), log_norm_table[,-c(1:2,8:9)], rowMeans(log_norm_table[,8:9]))
    colnames(log_norm_avg_table) = c('WT', 'CCR4', 'DHH1', 'LSM1', 'RPB4', 'SFP1', 'XRN1') 
    
    FC_mat = log_norm_avg_table[,-1]-log_norm_avg_table[,1]
    b = max(abs(median(log_norm_table[suff_reads,8] - log_norm_table[suff_reads,9])),
            abs(median(log_norm_table[suff_reads,1] - log_norm_table[suff_reads,2])))
            
    
    vioplot(FC_mat[suff_reads,1], FC_mat[suff_reads,2], FC_mat[suff_reads,3],
            FC_mat[suff_reads,4], FC_mat[suff_reads,5], FC_mat[suff_reads,6],
            names = c('CCR4', 'DHH1', 'LSM1', 'RPB4', 'SFP1', 'XRN1'), col = 'gray')
    abline(h = 0, col='red')
    abline(h = c(-b,b), col='green', lty = 2)
    title(main = paste('Fold changes in', ncRNA ,'transcription'), ylab = 'Log2 FC',
          xlab = 'Deletion')
    
    median_FCs = colMedians(FC_mat[suff_reads,])
    
    p = matrix(rep(NA, ncol(norm_table)^2), ncol = ncol(norm_table))
    colnames(p) = colnames(norm_table)
    rownames(p) = colnames(p)
    for(i in 1:(ncol(norm_table)-1)){
      for(j in (i+1):ncol(norm_table)){
        w_test = wilcox.test(x = log2(norm_table[suff_reads,i]+1), y = log2(norm_table[suff_reads,j]+1), paired = T, conf.int = T)
        p[i,j] = w_test$p.value
      }
    }
    to_return = list(p, DT, median_FCs, FC_mat, suff_reads)
    names(to_return) = c('p_values', 'DT_transcripts', 'median_FCs', 'FC_mat', 'sufficient_reads')
    
  } else if(data_flag == 2){
    suff_reads = which(norm_table[,1] >= 20)
    
    barplot(colMedians( (log_norm_table - log_norm_table[,1])[suff_reads,])[-c(1,2)],
            col = 'red', ylab = 'Median Log2 FC', main = paste(ncRNA, 'transcription'), args.legend = list(x = 22, y = .27, cex = .75))
    
    FC_mat = log_norm_table[,-c(1,2)]-log_norm_table[,1]
    vioplot(FC_mat[suff_reads,1], FC_mat[suff_reads,2], FC_mat[suff_reads,3],
            FC_mat[suff_reads,4], FC_mat[suff_reads,5],
            names = c('RCO1', 'DST1', 'EAF3', 'SET1', 'SET2'), col = 'gray')
    abline(h = 0, col='red')
    title(main = paste('Fold changes in', ncRNA ,'transcription'), ylab = 'Log2 FC',
          xlab = 'Deletion')
    
    median_FCs = colMedians(FC_mat[suff_reads,])
    
    p = matrix(rep(NA, ncol(norm_table)^2), ncol = ncol(norm_table))
    colnames(p) = colnames(norm_table)
    rownames(p) = colnames(p)
    for(i in 1:(ncol(norm_table)-1)){
      for(j in (i+1):ncol(norm_table)){
        p[i,j] = wilcox.test(x = norm_table[suff_reads,i], y = norm_table[suff_reads,j], paired = T)$p.value
      }
    }
    to_return = list(p, median_FCs, FC_mat, suff_reads)
    names(to_return) = c('p_values', 'median_FCs', 'FC_mat', 'sufficient_reads')
  }
  
  return(to_return)
  
}

chromosome_enrichment <- function(de_list, gene_annot){
  # Function to compute hypergeometric enrichment scores for chromosomal locations
  
  de_ind = which(ORF_annot$commonName %in% de_list)
  
  chrom = unique(gene_annot$chr)
  
  to_be_returned = NULL
  for(i in chrom){
    dt_hit = sum(gene_annot$chr[de_ind] == i)
    hit = sum(ORF_annot$chr == i)
    total = nrow(gene_annot)
    dt_total = length(de_list)
    
    # Compute the upper tail p-value using hypergeometric test
    hyper_value = phyper(dt_hit-1, hit, total-hit, dt_total, lower.tail = F)
    
    test_results = c(i, hyper_value, dt_hit, dt_total, dt_hit/dt_total, hit, total, hit/total)
    to_be_returned = rbind(to_be_returned, test_results)
    
  }
  
  colnames(to_be_returned) = c('Chromosome', 'P-value', 'DT genes on chrom', 
                                'DT genes', 'Fraction in DT', 'Genes on chrom', 'All genes', 'All fraction')
  rownames(to_be_returned) = NULL
  return(to_be_returned)
  
}

aneuploidy_checker = function(ref_ind = 12, mut_ind = 13, reads_table = all_reads_table){
  med_fc = NULL
  for(i in unique(reads_table$chr)){
    chrom_table = all_reads_table[all_reads_table$chr == i, c(ref_ind, mut_ind)]
    med_fc = c(med_fc, median(log2(chrom_table[,2]+1) - log2(chrom_table[,1]+1)))
  }
  plot(log2(reads_table[,mut_ind]+1) - log2(reads_table[,ref_ind]+1), xlab = 'Location',
       ylab = 'FC wrt WT')
  return(med_fc)
}

reads_scrubber = function(annot = tRNA_annot, plus_list = WT_plus, minus_list = WT_minus){
  
  for(i in 1:nrow(annot)){
    if(annot$chrom[i] == 'M'){
      chr = 17
    } else {
      chr = as.numeric(annot$chrom[i])
    }
    strand = annot$strand[i]
    if(strand == '+'){
      a = which(!is.na(match(plus_list[[chr]][,1], annot$start[i]:annot$end[i])))
      if(length(a) > 0){
        plus_list[[chr]][a,2] = 0
      }
    } else if (strand == '-'){
      a = which(!is.na(match(minus_list[[chr]][,1], annot$start[i]:annot$end[i])))
      if(length(a) > 0){
        minus_list[[chr]][a,2] = 0
      }
    }
  }
  
  to_return = list(plus_list, minus_list)
  names(to_return) = c('plus', 'minus')
  
  return(to_return)
  
}

profile_reads = function(start = ORF_annot$start[i], end = ORF_annot$end[i], chrom = ORF_annot$chr[i], strand = ORF_annot$strand[i], plus_list = WT_NC_plus, minus_list = WT_NC_minus, rad = c(100, 100), TSS_or_TES = 'TSS'){
  
  if(TSS_or_TES == 'TSS'){
    if(strand == '+'){
      site_range = (start-rad[1]) : (start+rad[2])
      s_dummy = as.matrix(plus_list[[chrom]][which(plus_list[[chrom]][,1] %in% site_range),])
      s_dummy = s_dummy[order(s_dummy[,1], decreasing = F),]
    } else if(strand == '-'){
      site_range = (end+rad[1]) : (end-rad[2])
      s_dummy = as.matrix(minus_list[[chrom]][which(minus_list[[chrom]][,1] %in% site_range),])
      s_dummy = s_dummy[order(s_dummy[,1], decreasing = T),]
      }
    sense_reads = cbind(site_range, rep(0, length(site_range)))
    if(!is.null(nrow(s_dummy))){
      sense_reads[which(sense_reads[,1] %in% s_dummy[,1]), 2] = s_dummy[,2]
      }
    }
  else if(TSS_or_TES == 'TES'){
      if(strand == '+'){
        site_range = (end-rad[1]) : (end+rad[2])
        s_dummy = as.matrix(plus_list[[chrom]][which(plus_list[[chrom]][,1] %in% site_range),])
        s_dummy = s_dummy[order(s_dummy[,1], decreasing = F),]
      } else if(strand == '-'){
        site_range = (start+rad[1]) : (start-rad[2])
        s_dummy = as.matrix(minus_list[[chrom]][which(minus_list[[chrom]][,1] %in% site_range),])
        s_dummy = s_dummy[order(s_dummy[,1], decreasing = T),]
      }
      sense_reads = cbind(site_range, rep(0, length(site_range)))
      if(!is.null(nrow(s_dummy))){
        sense_reads[which(sense_reads[,1] %in% s_dummy[,1]), 2] = s_dummy[,2]
      }
  }
  
  return(sense_reads[,2])
}

# function to produce matrix of gene profiles
gene_profile = function(annot = ORF_annot, plus_list = WT_NC_plus, minus_list = WT_NC_minus, rad = c(500, 500), TSS_or_TES = 'TSS', sense_or_antisense = 'sense'){
  if(sense_or_antisense == 'sense'){
    profile_mat = sapply(1:nrow(annot), function(i) profile_reads(start = annot$start[i], end = annot$end[i], chrom = annot$chr[i], strand = annot$strand[i], plus_list = plus_list, minus_list = minus_list, rad = rad, TSS_or_TES = TSS_or_TES))
  } else if (sense_or_antisense == 'antisense'){
    profile_mat = sapply(1:nrow(annot), function(i) profile_reads(start = annot$start[i], end = annot$end[i], chrom = annot$chr[i], strand = annot$strand[i], plus_list = minus_list, minus_list = plus_list, rad = rad, TSS_or_TES = TSS_or_TES))
  }
  colnames(profile_mat) = annot$name
  return(profile_mat)
}

# function to produce metagene plots
metagene = function(profile_mat = WT_pilot_TSS_prof, genes = NULL, annot = ORF_annot, coverage_mat = inclusion_mat){
  if(is.null(genes)){
    ind = 1:ncol(profile_mat)
  } else {
    ind = which(ORF_annot$name %in% genes | ORF_annot$commonName %in% genes)
  }
  if(is.null(coverage_mat)){
    avg_prof = apply(profile_mat[,ind], 1, function(x) mean(x, trim = .001))
  } else {
    avg_prof = rowSums((profile_mat[,ind] * coverage_mat[,ind]))/rowSums(coverage_mat[,ind])
  }
  x = 1:nrow(profile_mat)
  smoothed_prof = loess(avg_prof ~ x, span = .1)
  to_return = cbind(x-(nrow(profile_mat)-1)/2, smoothed_prof$fitted)
  return(to_return)
}

neighbor_flag = function(gene, annot = ORF_annot, rad = 500){
  ind = which(annot$name == gene)
  strand = annot$strand[ind]
  tss = (strand == '+')*annot$start[ind] + (strand == '-')*annot$end[ind]
  tes = (strand == '+')*annot$end[ind] + (strand == '-')*annot$start[ind]
  chrom = annot$chr[ind] 
  tss_sites = (strand == '+')*((tss-rad):(tss-1)) + (strand == '-')*((tss+1):(tss+rad))
  tes_sites = (strand == '+')*((tes+1):(tes+rad)) + (strand == '-')*((tss-rad):(tss-1))
  if(strand == '+'){
    tss_same_chrom_and_strand_start = annot$start[annot$chr == chrom & annot$strand == strand]
    tss_same_chrom_and_strand_end = annot$end[annot$chr == chrom & annot$strand == strand]
  }  else if(strand == '-'){
    tss_same_chrom_and_strand_start = annot$end[annot$chr == chrom & annot$strand == strand]
    tss_same_chrom_and_strand_end = annot$start[annot$chr == chrom & annot$strand == strand]
  }
  if(strand == '+'){
    tes_same_chrom_and_strand_start = annot$start[annot$chr == chrom & annot$strand == strand]
    tes_same_chrom_and_strand_end = annot$end[annot$chr == chrom & annot$strand == strand]
  }  else if(strand == '-'){
    tes_same_chrom_and_strand_start = annot$end[annot$chr == chrom & annot$strand == strand]
    tes_same_chrom_and_strand_end = annot$start[annot$chr == chrom & annot$strand == strand]
  }
  tss_flag = sum( (tss_same_chrom_and_strand_start %in% tss_sites) | (tss_same_chrom_and_strand_end %in% tss_sites) ) > 0
  tes_flag = sum( (tes_same_chrom_and_strand_start %in% tes_sites) | (tes_same_chrom_and_strand_end %in% tes_sites) ) > 0
  to_return = c(tss_flag, tes_flag)
  return(to_return)
}

non_orf_neighbor_flag = function(gene, annot = ORF_annot, rad = 500, non_orf_annot = trna_annot){
  ind = which(annot$name == gene)
  strand = annot$strand[ind]
  tss = (strand == '+')*annot$start[ind] + (strand == '-')*annot$end[ind]
  tes = (strand == '+')*annot$end[ind] + (strand == '-')*annot$start[ind]
  chrom = annot$chr[ind] 
  tss_sites = (strand == '+')*((tss-rad):(tss+rad)) + (strand == '-')*((tss-rad):(tss+rad))
  tes_sites = (strand == '+')*((tes-rad):(tes+rad)) + (strand == '-')*((tss-rad):(tss+rad))
  
  a = non_orf_annot[non_orf_annot$chr == chrom & non_orf_annot$strand == strand,]
  if(nrow(a) > 0){
    overlap_flag_tss = 0
    overlap_flag_tes = 0
    for(i in 1:nrow(a)){
      overlap_flag_tss = overlap_flag_tss + length(intersect((a$start[i]:a$end[i]), tss_sites)) > 0
      overlap_flag_tes = overlap_flag_tes + length(intersect((a$start[i]:a$end[i]), tes_sites)) > 0
    }
  } else {
    overlap_flag_tss = 0
    overlap_flag_tes = 0
  }
  
  to_return = c(overlap_flag_tss > 0, overlap_flag_tes > 0)
  return(to_return)
}

sum_collapser = function(start, stop, data){
  a = which(data[,1] > start & data[,1] < stop)
  s = sum(data[a,2])
  to_return = c(stop, s)
  names(to_return) = c('Coord', 'Sum')
  return(to_return)
}

reads_formatter = function(wig_file, norm = T){
  
  start_stop_points = c(which(wig_file[,1] == 'variableStep'), nrow(wig_file)+1)
  a = list()
  chrom_list = c('chrom=chrI', 'chrom=chrII', 'chrom=chrIII', 'chrom=chrIV', 'chrom=chrV', 'chrom=chrVI',
                 'chrom=chrVII', 'chrom=chrVIII', 'chrom=chrIX', 'chrom=chrX', 'chrom=chrXI', 'chrom=chrXII',
                 'chrom=chrXIII', 'chrom=chrXIV', 'chrom=chrXV', 'chrom=chrXVI', 'chrom=chrM')
  for(i in 1:(length(start_stop_points)-1)){
    j = which(chrom_list == wig_file[(start_stop_points[i]),2])
    a[[j]] = wig_file[(start_stop_points[i]+1):(start_stop_points[i+1]-1),]
  }
  
  a = lapply(a, function(x) cbind(as.numeric(as.character(x[,1])), as.numeric(as.character(x[,2]))))
  if(norm){
    norm_factor = min(unlist(lapply(a, function(x) min(x[,2]))))
    a = lapply(a, function(x) cbind(x[,1], x[,2]/norm_factor))
  }

  
  return(a)
}

profile_extractor = function(start, stop, chrom, strand, plus_list = WT_NC_plus, minus_list = WT_NC_minus){
  if(strand == '+'){
    a = which(plus_list[[chrom]][,1] >= start & plus_list[[chrom]][,1] <= stop)
    b = plus_list[[chrom]][a,]
  } else if (strand == '-'){
    a = which(minus_list[[chrom]][,1] >= start & minus_list[[chrom]][,1] <= stop)
    b = minus_list[[chrom]][a,]
  }
  
  if(!is.null(dim(b))){
    b[,1] = ((b[,1] - start) / (stop - start))*1000
    dilated_prof = t(sapply(0:999, function(x) sum_collapser(start = x, stop = x+1, data = b)))
    if(strand == '-'){
      dilated_prof[,2] = rev(dilated_prof[,2])
    }
  } else {
    dilated_prof = cbind(1:1000, 0)
  }
  
  return(dilated_prof)
}

my_fpca = function(mutant_mat, wt_mat, k1 = 1, k2 = NULL, to_log = T, make_dens = T, min_reads = 1){
  if(is.null(k2)){
    k2 = k1
  }
  if(make_dens){
    ind = colSums(mutant_mat) >= min_reads & colSums(wt_mat) >= min_reads
    mutant_mat = mutant_mat[,ind]
    mutant_mat = sweep(mutant_mat, 2, colSums(mutant_mat), '/')
    wt_mat = wt_mat[,ind]
    wt_mat = sweep(wt_mat, 2, colSums(wt_mat), '/')
  }
  if(to_log){
    new_mat = log2(mutant_mat+k2) - log2(wt_mat+k1)
  } else {
    new_mat = mutant_mat - wt_mat
  }
  N = ncol(new_mat)
  L3 = MakeFPCAInputs(IDs = rep(1:N, each = nrow(mutant_mat)), tVec = rep(1:nrow(mutant_mat), N), new_mat)
  FPCAdense <- FPCA(L3$Ly, L3$Lt)
  #plot(FPCAdense)
  return(FPCAdense)
}

fpca_prof_plots = function(fpca_output, input1, input2, name1, name2, prof_typ, norm_factor){
  
  if(prof_typ == 'TES'){
    leg_loc = 'topright'
  } else if(prof_typ == 'TSS'){
    leg_loc = 'topleft'
  }
  
  a = cut(fpca_output$xiEst[,1], breaks = quantile(fpca_output$xiEst[,1], c(0:5)/5), labels = F, include.lowest = T)
  
  par(mfrow = c(1,3))
  x1 = lapply(1:5, function(x) lowess(rowMeans(input1[,which(a == x)]), f = .05)$y)
  
  plot((-500:500), x1[[1]], typ = 'l', ylim = c(min(unlist(x1)), max(unlist(x1))),
       main = name1, xlab = paste('Distance from', prof_typ), ylab = 'Mean reads', col = 1)
  lines((-500:500), x1[[2]], col = 2)
  lines((-500:500), x1[[3]], col = 3)
  lines((-500:500), x1[[4]], col = 4)
  lines((-500:500), x1[[5]], col = 5)
  legend(x = leg_loc, lty = 1, col = c(1:5), legend = c(1:5), title = 'PC 1 Quintile', cex = .65)
  
  x2 = lapply(1:5, function(x) lowess(rowMeans(input2[,which(a == x)]/norm_factor), f = .05)$y)
  
  plot((-500:500), x2[[1]], typ = 'l', ylim = c(min(unlist(x2)), max(unlist(x2))),
       main = name2, xlab = paste('Distance from', prof_typ), ylab = 'Mean reads', col = 1)
  lines((-500:500), x2[[2]], col = 2)
  lines((-500:500), x2[[3]], col = 3)
  lines((-500:500), x2[[4]], col = 4)
  lines((-500:500), x2[[5]], col = 5)
  legend(x = leg_loc, lty = 1, col = c(1:5), legend = c(1:5), title = 'PC 1 Quintile', cex = .65)
  
  h = log2(input2/norm_factor+.001) - log2(input1+.001)
  h = lapply(1:5, function(x) lowess(rowMeans(h[,which(a == x)]), f = .05)$y)
  
  plot((-500:500), h[[1]], typ = 'l', ylim = c(min(unlist(h)), max(unlist(h))),
       main = paste(name2, 'vs', name1, 'FC'), xlab = paste('Distance from', prof_typ), ylab = 'Mean log2 FC', col = 1)
  lines((-500:500), h[[2]], col = 2)
  lines((-500:500), h[[3]], col = 3)
  lines((-500:500), h[[4]], col = 4)
  lines((-500:500), h[[5]], col = 5)
  legend(x = leg_loc, lty = 1, col = c(1:5), legend = c(1:5), title = 'PC 1 Quintile', cex = .65)
  
}

marginal_fpca = function(input_mat, to_dens = T){
  dens_mat = sweep(input_mat, 2, colSums(input_mat), '/')
  new_mat = dens_mat
  new_mat = new_mat[,which(colSums(input_mat) > 0)]
  new_mat2 = input_mat[,which(colSums(input_mat)>0)]
  N = ncol(new_mat)
  
  if(to_dens){
    L3 = MakeFPCAInputs(IDs = rep(1:N, each = 1001), tVec = rep(1:1001, N), new_mat)
  } else {
    L3 = MakeFPCAInputs(IDs = rep(1:N, each = 1001), tVec = rep(1:1001, N), new_mat2)
  }

  FPCAdense = FPCA(L3$Ly, L3$Lt)
  
  plot(FPCAdense)
  
  plot(rowMeans(new_mat[,FPCAdense$xiEst[,1] > 0 & FPCAdense$xiEst[,2] > 0]), typ = 'l')
  lines(rowMeans(new_mat[,FPCAdense$xiEst[,1] > 0 & FPCAdense$xiEst[,2] < 0]), col = 2)
  lines(rowMeans(new_mat[,FPCAdense$xiEst[,1] < 0 & FPCAdense$xiEst[,2] > 0]), col = 3)
  lines(rowMeans(new_mat[,FPCAdense$xiEst[,1] < 0 & FPCAdense$xiEst[,2] < 0]), col = 4)
  
  plot(rowMeans(new_mat2[,FPCAdense$xiEst[,1] > 0 & FPCAdense$xiEst[,2] > 0]), typ = 'l')
  lines(rowMeans(new_mat2[,FPCAdense$xiEst[,1] > 0 & FPCAdense$xiEst[,2] < 0]), col = 2)
  lines(rowMeans(new_mat2[,FPCAdense$xiEst[,1] < 0 & FPCAdense$xiEst[,2] > 0]), col = 3)
  lines(rowMeans(new_mat2[,FPCAdense$xiEst[,1] < 0 & FPCAdense$xiEst[,2] < 0]), col = 4)
  
  return(FPCAdense)
  
}

nearest_feature_finder = function(gene_vec = upstream_reads[1,], annot = master_annot, sense_or_antisense = 'sense', stream = 'up'){
  if(stream == 'up'){
    if(sense_or_antisense == 'sense'){
      strand = gene_vec$strand
      chr = gene_vec$chr
      sub_annot = annot[annot$chr == chr & annot$strand == strand, ]
      if( gene_vec$name %in% sub_annot$name ){
        matches = which(as.character(sub_annot$name) == as.character(gene_vec$name))
        sub_annot = sub_annot[-matches,]
      }
      if(strand == '+'){
        pt = gene_vec$start
        cand_pts = sub_annot$end
        cand_pts = cand_pts[cand_pts <= pt]
      } else if(strand == '-'){
        pt = gene_vec$end
        cand_pts = sub_annot$start
        cand_pts = cand_pts[cand_pts >= pt]
      }
      if(length(cand_pts) > 0){
        d = min(abs(cand_pts - pt))
      } else {
        d = NA
      }
    } else if (sense_or_antisense == 'antisense'){
      strand = setdiff(c('+', '-'), gene_vec$strand)
      chr = gene_vec$chr
      sub_annot = annot[annot$chr == chr & annot$strand == strand, ]
      if(gene_vec$strand == '+'){
        pt = gene_vec$start
        cand_pts = sub_annot$end
        cand_pts = cand_pts[cand_pts <= pt]
      } else if(gene_vec$strand == '-'){
        pt = gene_vec$end
        cand_pts = sub_annot$start
        cand_pts = cand_pts[cand_pts >= pt]
      }
      if(length(cand_pts) > 0){
        d = min(abs(cand_pts - pt))
      } else {
        d = NA
      }
    }
  } else if(stream == 'down'){
    if(sense_or_antisense == 'sense'){
      strand = gene_vec$strand
      chr = gene_vec$chr
      sub_annot = annot[annot$chr == chr & annot$strand == strand, ]
      if( gene_vec$name %in% sub_annot$name ){
        matches = which(as.character(sub_annot$name) == as.character(gene_vec$name))
        sub_annot = sub_annot[-matches,]
      }
      if(strand == '+'){
        pt = gene_vec$end
        cand_pts = sub_annot$start
        cand_pts = cand_pts[cand_pts >= pt]
      } else if(strand == '-'){
        pt = gene_vec$start
        cand_pts = sub_annot$end
        cand_pts = cand_pts[cand_pts <= pt]
      }
      if(length(cand_pts) > 0){
        d = min(abs(cand_pts - pt))
      } else {
        d = NA
      }
    } else if (sense_or_antisense == 'antisense'){
      strand = setdiff(c('+', '-'), gene_vec$strand)
      chr = gene_vec$chr
      sub_annot = annot[annot$chr == chr & annot$strand == strand, ]
      if(gene_vec$strand == '+'){
        pt = gene_vec$start
        cand_pts = c(sub_annot$start, sub_annot$end)
        cand_pts = cand_pts[cand_pts >= pt]
      } else if(gene_vec$strand == '-'){
        pt = gene_vec$end
        cand_pts = c(sub_annot$start, sub_annot$end)
        cand_pts = cand_pts[cand_pts <= pt]
      }
      if(length(cand_pts) > 0){
        d = min(abs(cand_pts - pt))
      } else {
        d = NA
      }
    }
  }
  
  return(d)
}

data_shifter = function(track, loc = 200701, shift_size = 355){
  track[,1][track[,1] >= loc] = track[,1][track[,1] >= loc]-shift_size
  dups = track[duplicated(track[,1]),1]
  for(x in dups){
    track[track[,1] == x, 2] = sum(track[track[,1] == x, 2])
  }
  track = track[!duplicated(track[,1]),]
  track = track[order(track[,1]),]
  return(track)
}

dens_vs_length = function(reads, column_ind){
  lengths = reads[,4]-reads[,3]
  a = reads[,column_ind]/lengths
  b = which(a <= quantile(a, .95))
  plot(lengths[b], a[b], xlab = 'Gene Length', ylab = 'Reads per bp', main = colnames(reads)[column_ind])
}

netseq_reads = function(gene_name, plus_list = WT_NC_plus, minus_list = WT_NC_minus, annot = ORF_annot){
  # gene_name: gene common name, locus name, or row index
  # plus_list and minus_list: lists containing number of (nonzero) reads at sites along the plus and minus strands
  # annot: annotation of gene names and locations
  
  # Finds requested gene in annotation file
  if(is.numeric(gene_name)){
    i = gene_name
  } else if(is.character(gene_name)){
    i = which(annot[,6] == toupper(gene_name))
    if(length(i) == 0){
      i = which(annot[,5] == toupper(gene_name))
    }
  }
  
  # Extracts chromosome and strand information for gene
  gene_chr = annot[i,1]
  gene_str = annot[i,2]
  
  # Get reads for sense and antisense strands
  if(gene_str == '+'){
    s_dummy = matrix(plus_list[[gene_chr]][which(plus_list[[gene_chr]][,1] %in% annot[i,3]:annot[i,4]),], ncol = 2)
    as_dummy = matrix(minus_list[[gene_chr]][which(minus_list[[gene_chr]][,1] %in% annot[i,3]:annot[i,4]),], ncol = 2)
  } else if(gene_str == '-'){
    s_dummy = matrix(minus_list[[gene_chr]][which(minus_list[[gene_chr]][,1] %in% annot[i,3]:annot[i,4]),], ncol = 2)
    as_dummy = matrix(plus_list[[gene_chr]][which(plus_list[[gene_chr]][,1] %in% annot[i,3]:annot[i,4]),], ncol = 2)
  }
  sense_reads = cbind(annot[i,3]:annot[i,4], rep(0, length(annot[i,3]:annot[i,4])))
  sense_reads[which(sense_reads[,1] %in% s_dummy[,1]), 2] = s_dummy[,2]
  
  antisense_reads = cbind(annot[i,3]:annot[i,4], rep(0, length(annot[i,3]:annot[i,4])))
  antisense_reads[which(antisense_reads[,1] %in% as_dummy[,1]), 2] = as_dummy[,2]
  
  if(gene_str == '-'){
    sense_reads[,1] = rev(sense_reads[,1])
    sense_reads[,2] = rev(sense_reads[,2])
    
    antisense_reads[,1] = rev(antisense_reads[,1])
    antisense_reads[,2] = rev(antisense_reads[,2])
  }
  
  return(list(sense_reads, antisense_reads))
}

first_k_reads = function(gene_ind = 1, read_prof = WT_NC_profiles, k = 1000, s_or_as = 1, direction = 1){
  # s_or_as: sense or antisense; 1 for s, 2 for as
  # direction: 1 for 5' to 3', -1 for 3' to 5'
  n = nrow(read_prof[[gene_ind]][[s_or_as]])
  if(direction == 1){
    a = read_prof[[gene_ind]][[s_or_as]][1:min(k,n),2]
    b = rep(1, k)
    if(n < k){
      a = c(a, rep(0, k-n))
      b[(n+1):k] = 0
      #print(paste('Gene is only', n, 'bases long'))
    }
  } else if(direction == -1){
    a = read_prof[[gene_ind]][[s_or_as]][max(1, n-k+1):n,2]
    b = rep(1, k)
    if(n < k){
      a = c(rep(0, k-n), a)
      b[1:(k-n)] = 0
      #print(paste('Gene is only', n, 'bases long'))
    }
  }
  return(cbind(b,a))
}

# be careful, have since changed profiles lists
metagene_plotter = function(annot = ORF_annot, prof_list = WT_NC_profiles, n_bp = 1000, s_flag = 1, direct = 1, strain_name = 'WT (nascent)'){
  first_k_list = lapply(1:nrow(annot), function(x) first_k_reads(gene_ind = x, read_prof = prof_list, k = n_bp, s_or_as = s_flag, direction = direct))
  list_sum = Reduce('+', first_k_list)
  dens_vec = cbind(1:n_bp, list_sum[,2]/list_sum[,1])
  smoothed_dens = lowess(dens_vec, f = 1/10)
  if(direct == 1){
    x_label = 'Distance downstream of TSS'
  } else if(direct == -1){
    x_label = paste('Final', n_bp, 'bases of the gene')
  }
  plot(dens_vec, xlab = x_label, ylab = 'Reads per bp', main = strain_name, pch = 1)
  lines(smoothed_dens, col = 'red')
  plot(smoothed_dens, typ = 'l', col = 'red', xlab = x_label, ylab = 'Reads per bp', main = paste(strain_name, 'peak at', which.max(smoothed_dens$y)))
  return(which.max(smoothed_dens$y))
}

compute_entropy = function(x){
  y = x[x>0]
  z = y/sum(y)
  h = -sum(z*log(z))
  return(h)
}

FPCA_plotter = function(input_mat, read_limit = 50, ent_limit = 4){
  ent_vec = apply(input_mat, 2, compute_entropy)
  dens_mat = sweep(input_mat, 2, colSums(input_mat), '/')
  indices = which(colSums(input_mat) > read_limit & ent_vec > ent_limit)
  new_mat = dens_mat[,indices]
  N = ncol(new_mat)
  M = nrow(new_mat)
  fpca_input = MakeFPCAInputs(IDs = rep(1:N, each = M), tVec = rep(1:M, N), new_mat)
  fpca_run = FPCA(fpca_input$Ly, fpca_input$Lt)
  b = list(new_mat, input_mat, indices)
  names(b) = c('Density', 'Reads', 'Indices')
  fpca_run = list(fpca_run, b)
  names(fpca_run) = c('FPCA', 'Inputs')
  return(fpca_run)
}

BDP_counter = function(gene = 'URA1', plus_list = WT_NC_plus, minus_list = WT_NC_minus, annot = ORF_annot){
  strand = annot$strand[annot$name == gene]
  chrom = annot$chr[annot$name == gene]
  
  plus_mat = plus_list[[chrom]]
  minus_mat = minus_list[[chrom]]
  
  if(strand == '+'){
    start_pt = annot[annot$name == gene,]$start
    end_pt = annot[annot$name == gene,]$end
    
    b = plus_mat[plus_mat[,1] >= start_pt & plus_mat[,1] < (start_pt + 700),]
    if(length(b) > 2){
      region = start_pt:(start_pt+699)
      d = region[!(region %in% b[,1])]
      d = cbind(d, rep(0, length(d)))
      b = rbind(b,d)
      b = b[order(b[,1]),]
      sum_vec = rep(0, 201)
      sum_vec[1] = sum(b[1:500, 2])
      for(i in 2:201){
        sum_vec[i] = sum_vec[i-1]-b[i-1,2]+b[i+499,2]
      }
      sense_sum = max(sum_vec)
    } else if(length(b) == 2){
      sense_sum = b[2]
    } else if(length(b) == 0){
      sense_sum = 0
    }
    
    b = minus_mat[minus_mat[,1] < (start_pt - 100) & minus_mat[,1] > (start_pt - 1100),]
    if(length(b) > 2){
      region = (start_pt-101):(start_pt-1100)
      d = region[!(region %in% b[,1])]
      d = cbind(d, rep(0, length(d)))
      b = rbind(b,d)
      b = b[order(b[,1]),]
      sum_vec = rep(0, 501)
      sum_vec[1] = sum(b[1:500, 2])
      for(i in 2:501){
        sum_vec[i] = sum_vec[i-1]-b[i-1,2]+b[i+499,2]
      }
      antisense_sum = max(sum_vec)
    } else if(length(b) == 2){
      antisense_sum = b[2]
    } else if(length(b) == 0){
      antisense_sum = 0
    }
    
  } else if(strand == '-'){
    start_pt = annot[annot$name == gene,]$end
    end_pt = annot[annot$name == gene,]$start
    
    b = minus_mat[minus_mat[,1] <= start_pt & minus_mat[,1] > (start_pt - 700),]
    if(length(b) > 2){
      region = start_pt:(start_pt-699)
      d = region[!(region %in% b[,1])]
      d = cbind(d, rep(0, length(d)))
      b = rbind(b,d)
      b = b[order(b[,1]),]
      sum_vec = rep(0, 201)
      sum_vec[1] = sum(b[1:500, 2])
      for(i in 2:201){
        sum_vec[i] = sum_vec[i-1]-b[i-1,2]+b[i+499,2]
      }
      sense_sum = max(sum_vec)
    } else if(length(b) == 2){
      sense_sum = b[2]
    } else if(length(b) == 0){
      sense_sum = 0
    }
    
    b = plus_mat[plus_mat[,1] > (start_pt + 100) & plus_mat[,1] < (start_pt + 1100),]
    if(length(b) > 2){
      region = (start_pt+101):(start_pt+1100)
      d = region[!(region %in% b[,1])]
      d = cbind(d, rep(0, length(d)))
      b = rbind(b,d)
      b = b[order(b[,1]),]
      sum_vec = rep(0, 501)
      sum_vec[1] = sum(b[1:500, 2])
      for(i in 2:501){
        sum_vec[i] = sum_vec[i-1]-b[i-1,2]+b[i+499,2]
      }
      antisense_sum = max(sum_vec)
    } else if(length(b) == 2){
      antisense_sum = b[2]
    } else if(length(b) == 0){
      antisense_sum = 0
    }
    
  }
  
  to_return = c(sense_sum, antisense_sum)
  names(to_return) = c('S_DS', 'AS_US')
  return(to_return)
}

as_neighbor = function(gene = 'GDH3', annot = ORF_annot, d = 1100, div_flag = T){
  
  chrom = annot$chr[annot$commonName == gene]
  strand = annot$strand[annot$commonName == gene]
  
  if(div_flag){
    if(strand == '+'){
      TSS = annot$start[annot$commonName == gene]
      AS_start_pts = annot$end[annot$chr == chrom & annot$strand == '-']
      num_neighbors = sum(AS_start_pts %in% (TSS):(TSS-d))
    } else if(strand == '-'){
      TSS = annot$end[annot$commonName == gene]
      AS_start_pts = annot$start[annot$chr == chrom & annot$strand == '+']
      num_neighbors = sum(AS_start_pts %in% (TSS):(TSS+d))
    }
  } else {
    if(strand == '+'){
      TSS = annot$start[annot$commonName == gene]
      TES = annot$end[annot$commonName == gene]
      AS_start_pts = annot$end[annot$chr == chrom & annot$strand == '-']
      AS_end_pts = annot$start[annot$chr == chrom & annot$strand == '-']
      num_neighbors = ceiling(.5*(sum(AS_start_pts %in% TSS:TES) + sum(AS_end_pts %in% TSS:TES)))
    } else if(strand == '-'){
      TSS = annot$end[annot$commonName == gene]
      TES = annot$start[annot$commonName == gene]
      AS_start_pts = annot$start[annot$chr == chrom & annot$strand == '+']
      AS_end_pts = annot$end[annot$chr == chrom & annot$strand == '+']
      num_neighbors = ceiling(.5*(sum(AS_start_pts %in% TSS:TES) + sum(AS_end_pts %in% TSS:TES)))
    }
  }
  
  return(num_neighbors)
  
}

overlap_finder = function(start = XUT_annot$start[1], end = XUT_annot$end[1], chrom = XUT_annot$chr[1], strand = XUT_annot$strand[1], comp_annot = update_annot){
  if(strand == '+'){
    comp_strand = '-'
  } else if(strand == '-'){
    comp_strand = '+'
  }
  comp_ind = which(comp_annot$strand == comp_strand & comp_annot$chr == chrom)
  a = lapply(1:length(comp_ind), function(x) intersect(comp_annot$start[comp_ind[x]]:comp_annot$end[comp_ind[x]], start:end)) 
  l = sapply(a, length)
  overlap_ind = comp_ind[which(l > 0)]
  overlap_sites = a[which(l>0)]
  return(overlap_ind)
}

overlap_full = function(my_annot = XUT_annot, ref_annot = update_annot){
  my_overlaps = sapply(1:nrow(my_annot), function(x) overlap_finder(start = my_annot$start[x], end = my_annot$end[x], chrom = my_annot$chr[x], strand = my_annot$strand[x], comp_annot = ref_annot))
  dd = unlist(my_overlaps)
  b = sapply(my_overlaps, length)
  d = unlist(sapply(1:length(my_overlaps), function(x) rep(x, b[x])))
  a = cbind(dd, d)
  colnames(a) = c('Overlap', 'Ind')
  return(a)
}

gene_neighbors = function(start = ORF_annot$start[1], end = ORF_annot$end[1], chrom = ORF_annot$chr[1], strand = ORF_annot$strand[1], radius = 500, cov_list = site_coverage_list){
  if(strand == '+'){
    my_strand_ind = 1
    as_strand_ind = 2
  } else if(strand == '-'){
    my_strand_ind = 2
    as_strand_ind = 1
  }
  my_strand = cov_list[[chrom]][[my_strand_ind]]
  as_strand = cov_list[[chrom]][[as_strand_ind]]
  flag_vec = rep(0,6)
  names(flag_vec) = c('S_upstream', 'S_downstream', 'AS_upstream', 'AS_downstream', 'AS_overlap', 'divergent')
  
  if(strand == '+'){
    flag_vec[1] = sum( (start-1):(start-501) %in% my_strand ) > 0
    flag_vec[2] = sum( (end+1):(end+501) %in% my_strand ) > 0
    flag_vec[3] = sum( (start-1):(start-501) %in% as_strand ) > 0
    flag_vec[4] = sum( (end+1):(end+501) %in% as_strand ) > 0
    flag_vec[5] = sum( start:end %in% as_strand ) > 0
    flag_vec[6] = (!start %in% as_strand) & flag_vec[3]
  } else if(strand == '-'){
    flag_vec[1] = sum( (end+1):(end+501) %in% my_strand ) > 0
    flag_vec[2] = sum( (start-1):(start-501) %in% my_strand ) > 0
    flag_vec[3] = sum( (end+1):(end+501) %in% as_strand ) > 0
    flag_vec[4] = sum( (start-1):(start-501) %in% as_strand ) > 0
    flag_vec[5] = sum( start:end %in% as_strand ) > 0
    flag_vec[6] = (!end %in% as_strand) & flag_vec[3]
  }
  flag_vec = flag_vec*1
  
  return(flag_vec)
  
}

new_fct = function(fpca_output = WT_NC_FPCA_TSS, prof_mat = WT_NC_TSS_prof[,which(colSums(WT_NC_TSS_prof) > 0)], comp_num = 2, n_q = 5, input_samp = ''){
  r = cut(fpca_output$xiEst[,comp_num], quantile(fpca_output$xiEst[,comp_num], c(0:n_q)/n_q), include.lowest = T, labels = F)
  lowess_list = lapply(1:n_q, function(x) lowess(-500:500, rowMeans(prof_mat[1:1001, r == x]), f = .05))
  y_max = max(unlist(lapply(lowess_list, function(x) x[[2]])))
  
  plot(lowess_list[[1]], typ = 'l', ylim = c(0, y_max), xlab = 'Distance from TSS', 
       ylab = 'Mean reads', main = paste(input_samp, 'Metagenes by loading bins'))
  for(j in 2:n_q){
    lines(lowess_list[[j]], col = j)
  }
  legend(x = 'topleft', col = 1:n_q, legend = paste('Bin', 1:n_q), lwd = 1, cex = .5,
         title = 'Loading bin')
  return(lowess_list)
}

coverage_list = function(annot = big_annot){
  a = lapply(1:16, function(x) lapply(c('+', '-'), function(y) annot[annot$chr == x & annot$strand == y, 3:4]))
  b = lapply(a, function(x) lapply(x, function(y) unlist(sapply(1:nrow(y), function(z) y[z,1]:y[z,2]))))
  return(b)
}

# smooth, rescale, add, decompose
prof_maker = function(prof, strand = NULL){
  prof = prof[,2]
  L = length(prof)
  if(L >= 1000){
    #if(strand == '-'){
    #  prof = rev(prof)
    #}
    five_prime = prof[1:500]
    three_prime = prof[(L-199):L]
    middle = prof[-c(1:500, (L-199):L)]
    x = c(0, 1:length(middle), length(middle)+1)
    z = approx(x = x/length(middle)*500, y = c(prof[500], middle, prof[L-199]), xout = 1:500)
    new_prof = c(five_prime, z$y, three_prime)
  } else {
    new_prof = NULL
  }
  return(new_prof)
}

# compute Jaccard similarity
jaccard = function(x,y){
  j = length(intersect(x,y))/length(union(x,y))
  return(j)
}


augmented_annot = ORF_annot
augmented_annot$start = ORF_annot$start - 201
augmented_annot$end = ORF_annot$end + 201

WT_list = sapply(1:4973, function(i) netseq_reads(i, plus_list = WT_plus, minus_list = WT_minus, annot = augmented_annot)[[1]][,2])
WT_pilot_list = sapply(1:4973, function(i) netseq_reads(i, plus_list = WT_pilot_plus, minus_list = WT_pilot_minus, annot = augmented_annot)[[1]][,2])
CCR4_list = sapply(1:4973, function(i) netseq_reads(i, plus_list = CCR4_plus, minus_list = CCR4_minus, annot = augmented_annot)[[1]][,2])
DHH1_list = sapply(1:4973, function(i) netseq_reads(i, plus_list = DHH1_plus, minus_list = DHH1_minus, annot = augmented_annot)[[1]][,2])
LSM1_list = sapply(1:4973, function(i) netseq_reads(i, plus_list = LSM1_plus, minus_list = LSM1_minus, annot = augmented_annot)[[1]][,2])
RPB4_list = sapply(1:4973, function(i) netseq_reads(i, plus_list = RPB4_plus, minus_list = RPB4_minus, annot = augmented_annot)[[1]][,2])
XRN1_list = sapply(1:4973, function(i) netseq_reads(i, plus_list = XRN1_plus, minus_list = XRN1_minus, annot = augmented_annot)[[1]][,2])
XRN1_pilot_list = sapply(1:4973, function(i) netseq_reads(i, plus_list = XRN1_pilot_plus, minus_list = XRN1_pilot_minus, annot = augmented_annot)[[1]][,2])


dilated_metagenes = function(gene_list, cutoff = 1000){
  l = sapply(gene_list, length)
  a = gene_list[l >= 1500]
  fronts = sapply(a, function(x) x[1:701])
  middles = lapply(a, function(x) x[702:(length(x)-701)])
  ends = sapply(a, function(x) x[(length(x)-700):length(x)])
  
  compressed_meta_mat = matrix(0, nrow = 1902, ncol = ncol(fronts))
  for(j in 1:ncol(fronts)){
    x = c(0, 1:length(middles[[j]]), length(middles[[j]])+1)
    z = approx(x = x/length(middles[[j]])*500, y = c(fronts[701,j], middles[[j]], ends[1,j]), xout = 1:500)
    compressed_meta_mat[,j] = c(fronts[,j], z$y, ends[,j])
  }
  
  compressed_meta_mat[compressed_meta_mat >= cutoff] = 0
  
  return(compressed_meta_mat)
}

WT_dil = dilated_metagenes(WT_list, cutoff = 100)
WT_pilot_dil = dilated_metagenes(WT_pilot_list, cutoff = 100)
CCR4_dil = dilated_metagenes(CCR4_list, cutoff = 500)
DHH1_dil = dilated_metagenes(DHH1_list, cutoff = 500)
LSM1_dil = dilated_metagenes(LSM1_list, cutoff = 500)
RPB4_dil = dilated_metagenes(RPB4_list, cutoff = 500)
XRN1_dil = dilated_metagenes(XRN1_list, cutoff = 100)
XRN1_pilot_dil = dilated_metagenes(XRN1_pilot_list, cutoff = 100)

pdf(file = 'full_metagenes_1.pdf', width = 12, height = 5)
par(mar = c(5,5,5,3))
plot(rowMeans(WT_dil + WT_pilot_dil)/(nf[1]+nf[2]), typ = 'l', col = 'gray',
     xlab = 'Metaposition', ylab = 'Mean NET-seq reads', main = 'Full metagenes',
     ylim = c(0,1), xaxt = 'n', cex.main = 2, cex.axis = 1.5, cex.lab = 1.5, lwd = 2)
lines(rowMeans(XRN1_dil + XRN1_pilot_dil)/(nf[8]+nf[9]), typ = 'l', col = 'black', lwd = 3)
legend(x = 'topright', col = c('gray', 'black'), lwd = 2, legend = c('WT', expression(paste(italic('xrn1'), Delta))), cex = 1.5, ncol = 1)
axis(1, at = c(1, 201, 701, 1201, 1702, 1902), labels = c('TSS-200', 'TSS', 'TSS+500', 'PAS-500', 'PAS', 'PAS+200'), cex.axis = 1.5)
dev.off()

pdf(file = 'full_metagenes_2.pdf', width = 12, height = 5)
par(mar = c(5,5,5,3))
plot(rowMeans(WT_dil + WT_pilot_dil)/(nf[1]+nf[2]), typ = 'l', col = 'gray',
     xlab = 'Metaposition', ylab = 'Mean NET-seq reads', main = 'Full metagenes',
     ylim = c(0,1), xaxt = 'n', cex.main = 2, cex.axis = 1.5, cex.lab = 1.5, lwd = 2)
lines(rowMeans(DHH1_dil)/(nf[4]), typ = 'l', col = 'red')
lines(rowMeans(LSM1_dil)/(nf[5]), typ = 'l', col = 'green')
legend(x = 'topright', col = c('gray', 'red', 'green'), lwd = 2, legend = c('WT', expression(paste(italic('dhh1'), Delta)), expression(paste(italic('lsm1'), Delta))), cex = 1.5, ncol = 1)
axis(1, at = c(1, 201, 701, 1201, 1702, 1902), labels = c('TSS-200', 'TSS', 'TSS+500', 'PAS-500', 'PAS', 'PAS+200'), cex.axis = 1.5)
dev.off()

pdf(file = 'full_metagenes_3.pdf', width = 12, height = 5)
par(mar = c(5,5,5,3))
plot(rowMeans(WT_dil + WT_pilot_dil)/(nf[1]+nf[2]), typ = 'l', col = 'gray',
     xlab = 'Metaposition', ylab = 'Mean NET-seq reads', main = 'Full metagenes',
     ylim = c(0,1), xaxt = 'n', cex.main = 2, cex.axis = 1.5, cex.lab = 1.5, lwd = 2)
lines(rowMeans(CCR4_dil)/(nf[3]), typ = 'l', col = 'blue')
lines(rowMeans(RPB4_dil)/(nf[6]), typ = 'l', col = 'magenta')
legend(x = 'topright', col = c('gray', 'blue', 'magenta'), lwd = 2, legend = c('WT', expression(paste(italic('ccr4'), Delta)), expression(paste(italic('rpb4'), Delta))), cex = 1.5, ncol = 1)
axis(1, at = c(1, 201, 701, 1201, 1702, 1902), labels = c('TSS-200', 'TSS', 'TSS+500', 'PAS-500', 'PAS', 'PAS+200'), cex.axis = 1.5)
dev.off()














augmented_annot = XUT_annot
augmented_annot$start = XUT_annot$start - 201
augmented_annot$end = XUT_annot$end + 201

WT_list = sapply(1:4973, function(i) netseq_reads(i, plus_list = WT_plus, minus_list = WT_minus, annot = augmented_annot)[[1]][,2])
WT_pilot_list = sapply(1:4973, function(i) netseq_reads(i, plus_list = WT_pilot_plus, minus_list = WT_pilot_minus, annot = augmented_annot)[[1]][,2])
CCR4_list = sapply(1:4973, function(i) netseq_reads(i, plus_list = CCR4_plus, minus_list = CCR4_minus, annot = augmented_annot)[[1]][,2])
DHH1_list = sapply(1:4973, function(i) netseq_reads(i, plus_list = DHH1_plus, minus_list = DHH1_minus, annot = augmented_annot)[[1]][,2])
LSM1_list = sapply(1:4973, function(i) netseq_reads(i, plus_list = LSM1_plus, minus_list = LSM1_minus, annot = augmented_annot)[[1]][,2])
RPB4_list = sapply(1:4973, function(i) netseq_reads(i, plus_list = RPB4_plus, minus_list = RPB4_minus, annot = augmented_annot)[[1]][,2])
XRN1_list = sapply(1:4973, function(i) netseq_reads(i, plus_list = XRN1_plus, minus_list = XRN1_minus, annot = augmented_annot)[[1]][,2])
XRN1_pilot_list = sapply(1:4973, function(i) netseq_reads(i, plus_list = XRN1_pilot_plus, minus_list = XRN1_pilot_minus, annot = augmented_annot)[[1]][,2])

WT_dil = dilated_metagenes(WT_list, cutoff = 100)
WT_pilot_dil = dilated_metagenes(WT_pilot_list, cutoff = 100)
CCR4_dil = dilated_metagenes(CCR4_list, cutoff = 500)
DHH1_dil = dilated_metagenes(DHH1_list, cutoff = 500)
LSM1_dil = dilated_metagenes(LSM1_list, cutoff = 500)
RPB4_dil = dilated_metagenes(RPB4_list, cutoff = 500)
XRN1_dil = dilated_metagenes(XRN1_list, cutoff = 100)
XRN1_pilot_dil = dilated_metagenes(XRN1_pilot_list, cutoff = 100)

pdf(file = 'XUT_metagenes.pdf', width = 12, height = 5)
par(mar = c(5,5,5,3))
plot(rowMeans(WT_dil + WT_pilot_dil)/(nf[1]+nf[2]), typ = 'l', col = 'gray',
     xlab = 'Metaposition', ylab = 'Mean NET-seq reads', main = 'Full metagenes',
     ylim = c(0,1), xaxt = 'n', cex.main = 2, cex.axis = 1.5, cex.lab = 1.5, lwd = 2)
lines(rowMeans(XRN1_dil + XRN1_pilot_dil)/(nf[8]+nf[9]), typ = 'l', col = 'black', lwd = 3)
lines(rowMeans(DHH1_dil)/(nf[4]), typ = 'l', col = 'red')
lines(rowMeans(LSM1_dil)/(nf[5]), typ = 'l', col = 'green')
lines(rowMeans(CCR4_dil)/(nf[3]), typ = 'l', col = 'blue')
lines(rowMeans(RPB4_dil)/(nf[6]), typ = 'l', col = 'magenta')
axis(1, at = c(1, 201, 701, 1201, 1702, 1902), labels = c('TSS-200', 'TSS', 'TSS+500', 'PAS-500', 'PAS', 'PAS+200'), cex.axis = 1.5)
legend(x = 'topright', col = c('gray', 'black'), lwd = 2, legend = c('WT', expression(paste(italic('xrn1'), Delta))), cex = 1.5, ncol = 1)
legend(x = 'topright', col = c('gray', 'red', 'green'), lwd = 2, legend = c('WT', expression(paste(italic('dhh1'), Delta)), expression(paste(italic('lsm1'), Delta))), cex = 1.5, ncol = 1)
legend(x = 'topright', col = c('gray', 'blue', 'magenta'), lwd = 2, legend = c('WT', expression(paste(italic('ccr4'), Delta)), expression(paste(italic('rpb4'), Delta))), cex = 1.5, ncol = 1)
dev.off()


