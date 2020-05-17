library(plyr)
source('code/clean_code/netseq_functions.R')

# annotation for Churchman data
ORF_annot_c = read.delim("~/Documents/R_Projects/NETSEQ/data/ORF_annot.txt")

# combine and reformat annotations for our data
TSS_annot = read.delim('data/GSE49026_S-TSS.txt')
PAS_annot = read.delim('data/GSE49026_S-PAS.txt')
joint_annot = merge(TSS_annot, PAS_annot, by = c('chr', 'ORF'))[,c(2,3,6)]
colnames(joint_annot) = c('ORF', 'TSS', 'PAS')
a = sapply(as.character(joint_annot[,1]), function(x) substr(x, 2, 2))
b = as.numeric(mapvalues(a, c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P'), 1:16))
joint_annot = cbind(joint_annot, b)
colnames(joint_annot)[4] = 'chr'
s = sapply(1:nrow(joint_annot), function(x) if(joint_annot$TSS[x] < joint_annot$PAS[x]){'+'}else{'-'})
joint_annot = cbind(joint_annot, s)
colnames(joint_annot)[5] = 'strand'
joint_annot = joint_annot[,c(4,5,2,3,1)]
h = sapply(1:nrow(joint_annot), function(x) min(joint_annot$TSS[x], joint_annot$PAS[x]))
j = sapply(1:nrow(joint_annot), function(x) max(joint_annot$TSS[x], joint_annot$PAS[x]))
ORF_annot = cbind(joint_annot[,1:2], h, j, joint_annot[,5])
colnames(ORF_annot) = c('chr', 'strand', 'start', 'end', 'name')
ORF_annot$commonName = as.character(ORF_annot$name)
ORF_annot$commonName[which(ORF_annot$name %in% ORF_annot_c$name)] = as.character(ORF_annot_c$commonName[match(ORF_annot$name[which(ORF_annot$name %in% ORF_annot_c$name)], ORF_annot_c$name)])

# load annotations for regions other than ORFs
tRNA_annot = read.table('~/Documents/R_Projects/NETSEQ/data/tRNA_annot.txt', header = T)
snRNA_annot = read.table('~/Documents/R_Projects/NETSEQ/data/snRNA_annot.txt', header = T)
ARS_annot = read.table('~/Documents/R_Projects/NETSEQ/data/ARS_annot.txt', header = T)

CUT_annot = read.table('~/Documents/R_Projects/NETSEQ/data/CUT_annot.txt', header = T)[,c(2,3,4,5)]
NUT_annot = read.table('~/Documents/R_Projects/NETSEQ/data/NUT_annot.txt', header = T)
SUT_annot = read.table('~/Documents/R_Projects/NETSEQ/data/SUT_annot.txt', header = T)[,c(2,3,4,5)]
XUT_annot = read.table('~/Documents/R_Projects/NETSEQ/data/XUT_annot.txt', header = T)

ncRNA_annot = rbind(CUT_annot, NUT_annot, SUT_annot, XUT_annot)

# load NET-seq data from Churchman paper and convert it to a workable format
WT_NC_plus = read.csv("~/Documents/R_Projects/NETSEQ/data/GSE25107_RAW/GSM617027_WT_NC_plus.wig", sep="")
WT_mRNA_plus = read.csv("~/Documents/R_Projects/NETSEQ/data/GSE25107_RAW/GSM617028_WT_mRNA_plus.wig", sep="")
RCO1_plus = read.csv("~/Documents/R_Projects/NETSEQ/data/GSE25107_RAW/GSM617029_RCO1D_plus.wig", sep="")
DST1_plus = read.csv("~/Documents/R_Projects/NETSEQ/data/GSE25107_RAW/GSM617030_DST1D_plus.wig", sep="")
EAF3_plus = read.csv("~/Documents/R_Projects/NETSEQ/data/GSE25107_RAW/GSM617031_EAF3D_plus.wig", sep="")
SET1_plus = read.csv("~/Documents/R_Projects/NETSEQ/data/GSE25107_RAW/GSM617033_SET1D_plus.wig", sep="")
SET2_plus = read.csv("~/Documents/R_Projects/NETSEQ/data/GSE25107_RAW/GSM617032_SET2D_plus.wig", sep="")

WT_NC_minus = read.csv("~/Documents/R_Projects/NETSEQ/data/GSE25107_RAW/GSM617027_WT_NC_minus.wig", sep="")
WT_mRNA_minus = read.csv("~/Documents/R_Projects/NETSEQ/data/GSE25107_RAW/GSM617028_WT_mRNA_minus.wig", sep="")
RCO1_minus = read.csv("~/Documents/R_Projects/NETSEQ/data/GSE25107_RAW/GSM617029_RCO1D_minus.wig", sep="")
DST1_minus = read.csv("~/Documents/R_Projects/NETSEQ/data/GSE25107_RAW/GSM617030_DST1D_minus.wig", sep="")
EAF3_minus = read.csv("~/Documents/R_Projects/NETSEQ/data/GSE25107_RAW/GSM617031_EAF3D_minus.wig", sep="")
SET1_minus = read.csv("~/Documents/R_Projects/NETSEQ/data/GSE25107_RAW/GSM617033_SET1D_minus.wig", sep="")
SET2_minus = read.csv("~/Documents/R_Projects/NETSEQ/data/GSE25107_RAW/GSM617032_SET2D_minus.wig", sep="")


WT_NC_plus = reads_formatter(WT_NC_plus)
WT_NC_minus = reads_formatter(WT_NC_minus)

WT_mRNA_plus = reads_formatter(WT_mRNA_plus)
WT_mRNA_minus = reads_formatter(WT_mRNA_minus)

RCO1_plus = reads_formatter(RCO1_plus)
RCO1_minus = reads_formatter(RCO1_minus)

DST1_plus = reads_formatter(DST1_plus)
DST1_minus = reads_formatter(DST1_minus)

EAF3_plus = reads_formatter(EAF3_plus)
EAF3_minus = reads_formatter(EAF3_minus)

SET1_plus = reads_formatter(SET1_plus)
SET1_minus = reads_formatter(SET1_minus)

SET2_plus = reads_formatter(SET2_plus)
SET2_minus = reads_formatter(SET2_minus)

# load our data and convert it to a workable format
CCR4_plus = read.csv("~/Documents/R_Projects/NETSEQ/data/wigs/ccr4.plus.wig", sep="")
CCR4_minus = read.csv("~/Documents/R_Projects/NETSEQ/data/wigs/ccr4.minus.wig", sep="")
DHH1_plus = read.csv("~/Documents/R_Projects/NETSEQ/data/wigs/dhh1.plus.wig", sep="")
DHH1_minus = read.csv("~/Documents/R_Projects/NETSEQ/data/wigs/dhh1.minus.wig", sep="")
LSM1_plus = read.csv("~/Documents/R_Projects/NETSEQ/data/wigs/lsm1.plus.wig", sep="")
LSM1_minus = read.csv("~/Documents/R_Projects/NETSEQ/data/wigs/lsm1.minus.wig", sep="")
RPB4_plus = read.csv("~/Documents/R_Projects/NETSEQ/data/wigs/rpb4.plus.wig", sep="")
RPB4_minus = read.csv("~/Documents/R_Projects/NETSEQ/data/wigs/rpb4.minus.wig", sep="")
WT_plus = read.csv("~/Documents/R_Projects/NETSEQ/data/wigs/WT.plus.wig", sep="")
WT_minus = read.csv("~/Documents/R_Projects/NETSEQ/data/wigs/WT.minus.wig", sep="")
WT_pilot_plus = read.csv("~/Documents/R_Projects/NETSEQ/data/wigs/WTPilot.plus.wig", sep="")
WT_pilot_minus = read.csv("~/Documents/R_Projects/NETSEQ/data/wigs/WTPilot.minus.wig", sep="")
XRN1_plus = read.csv("~/Documents/R_Projects/NETSEQ/data/wigs/xrn1.plus.wig", sep="")
XRN1_minus = read.csv("~/Documents/R_Projects/NETSEQ/data/wigs/xrn1.minus.wig", sep="")
XRN1_pilot_plus = read.csv("~/Documents/R_Projects/NETSEQ/data/wigs/xrn1Pilot.plus.wig", sep="")
XRN1_pilot_minus = read.csv("~/Documents/R_Projects/NETSEQ/data/wigs/xrn1Pilot.minus.wig", sep="")


CCR4_plus = reads_formatter(CCR4_plus)
CCR4_minus = reads_formatter(CCR4_minus)
DHH1_plus = reads_formatter(DHH1_plus)
DHH1_minus = reads_formatter(DHH1_minus)
LSM1_plus = reads_formatter(LSM1_plus)
LSM1_minus = reads_formatter(LSM1_minus)
RPB4_plus = reads_formatter(RPB4_plus)
RPB4_minus = reads_formatter(RPB4_minus)
WT_plus = reads_formatter(WT_plus)
WT_minus = reads_formatter(WT_minus)
WT_pilot_plus = reads_formatter(WT_pilot_plus)
WT_pilot_minus = reads_formatter(WT_pilot_minus)
XRN1_plus = reads_formatter(XRN1_plus)
XRN1_minus = reads_formatter(XRN1_minus)
XRN1_pilot_plus = reads_formatter(XRN1_pilot_plus)
XRN1_pilot_minus = reads_formatter(XRN1_pilot_minus)

CCR4_plus[[8]][,2] = CCR4_plus[[8]][,2]/2
CCR4_minus[[8]][,2] = CCR4_minus[[8]][,2]/2
XRN1_plus[[11]][,2] = XRN1_plus[[11]][,2]/2
XRN1_minus[[11]][,2] = XRN1_minus[[11]][,2]/2
XRN1_pilot_plus[[11]][,2] = XRN1_pilot_plus[[11]][,2]/2
XRN1_pilot_minus[[11]][,2] = XRN1_pilot_minus[[11]][,2]/2

# Remove reads from problematic regions
CCR4_clean = reads_scrubber(annot = tRNA_annot, plus_list = CCR4_plus, minus_list = CCR4_minus)
DHH1_clean = reads_scrubber(annot = tRNA_annot, plus_list = DHH1_plus, minus_list = DHH1_minus)
LSM1_clean = reads_scrubber(annot = tRNA_annot, plus_list = LSM1_plus, minus_list = LSM1_minus)
RPB4_clean = reads_scrubber(annot = tRNA_annot, plus_list = RPB4_plus, minus_list = RPB4_minus)
SFP1_clean = reads_scrubber(annot = tRNA_annot, plus_list = SFP1_plus, minus_list = SFP1_minus)
WT_clean = reads_scrubber(annot = tRNA_annot, plus_list = WT_plus, minus_list = WT_minus)
WT_pilot_clean = reads_scrubber(annot = tRNA_annot, plus_list = WT_pilot_plus, minus_list = WT_pilot_minus)
XRN1_clean = reads_scrubber(annot = tRNA_annot, plus_list = XRN1_plus, minus_list = XRN1_minus)
XRN1_pilot_clean = reads_scrubber(annot = tRNA_annot, plus_list = XRN1_pilot_plus, minus_list = XRN1_pilot_minus)

CCR4_clean = reads_scrubber(annot = snRNA_annot, plus_list = CCR4_clean[[1]], minus_list = CCR4_clean[[2]])
DHH1_clean = reads_scrubber(annot = snRNA_annot, plus_list = DHH1_clean[[1]], minus_list = DHH1_clean[[2]])
LSM1_clean = reads_scrubber(annot = snRNA_annot, plus_list = LSM1_clean[[1]], minus_list = LSM1_clean[[2]])
RPB4_clean = reads_scrubber(annot = snRNA_annot, plus_list = RPB4_clean[[1]], minus_list = RPB4_clean[[2]])
SFP1_clean = reads_scrubber(annot = snRNA_annot, plus_list = SFP1_clean[[1]], minus_list = SFP1_clean[[2]])
WT_clean = reads_scrubber(annot = snRNA_annot, plus_list = WT_clean[[1]], minus_list = WT_clean[[2]])
WT_pilot_clean = reads_scrubber(annot = snRNA_annot, plus_list = WT_pilot_clean[[1]], minus_list = WT_pilot_clean[[2]])
XRN1_clean = reads_scrubber(annot = snRNA_annot, plus_list = XRN1_clean[[1]], minus_list = XRN1_clean[[2]])
XRN1_pilot_clean = reads_scrubber(annot = snRNA_annot, plus_list = XRN1_pilot_clean[[1]], minus_list = XRN1_pilot_clean[[2]])

CCR4_clean = reads_scrubber(annot = ARS_annot, plus_list = CCR4_clean[[1]], minus_list = CCR4_clean[[2]])
DHH1_clean = reads_scrubber(annot = ARS_annot, plus_list = DHH1_clean[[1]], minus_list = DHH1_clean[[2]])
LSM1_clean = reads_scrubber(annot = ARS_annot, plus_list = LSM1_clean[[1]], minus_list = LSM1_clean[[2]])
RPB4_clean = reads_scrubber(annot = ARS_annot, plus_list = RPB4_clean[[1]], minus_list = RPB4_clean[[2]])
SFP1_clean = reads_scrubber(annot = ARS_annot, plus_list = SFP1_clean[[1]], minus_list = SFP1_clean[[2]])
WT_clean = reads_scrubber(annot = ARS_annot, plus_list = WT_clean[[1]], minus_list = WT_clean[[2]])
WT_pilot_clean = reads_scrubber(annot = ARS_annot, plus_list = WT_pilot_clean[[1]], minus_list = WT_pilot_clean[[2]])
XRN1_clean = reads_scrubber(annot = ARS_annot, plus_list = XRN1_clean[[1]], minus_list = XRN1_clean[[2]])
XRN1_pilot_clean = reads_scrubber(annot = ARS_annot, plus_list = XRN1_pilot_clean[[1]], minus_list = XRN1_pilot_clean[[2]])

# aggregate counts in ORF regions for Churchman data
WT_NC_counts = t(sapply(1:nrow(ORF_annot_c), function(x) netseq_plotter(x, plot_flag = 0, annot = ORF_annot_c)))
WT_mRNA_counts = t(sapply(1:nrow(ORF_annot_c), function(x) netseq_plotter(x, plus_list = WT_mRNA_plus, minus_list = WT_mRNA_minus, plot_flag = 0, annot = ORF_annot_c)))
RCO1_counts = t(sapply(1:nrow(ORF_annot_c), function(x) netseq_plotter(x, plus_list = RCO1_plus, minus_list = RCO1_minus, plot_flag = 0, annot = ORF_annot_c)))
DST1_counts = t(sapply(1:nrow(ORF_annot_c), function(x) netseq_plotter(x, plus_list = DST1_plus, minus_list = DST1_minus, plot_flag = 0, annot = ORF_annot_c)))
EAF3_counts = t(sapply(1:nrow(ORF_annot_c), function(x) netseq_plotter(x, plus_list = EAF3_plus, minus_list = EAF3_minus, plot_flag = 0, annot = ORF_annot_c)))
SET1_counts = t(sapply(1:nrow(ORF_annot_c), function(x) netseq_plotter(x, plus_list = SET1_plus, minus_list = SET1_minus, plot_flag = 0, annot = ORF_annot_c)))
SET2_counts = t(sapply(1:nrow(ORF_annot_c), function(x) netseq_plotter(x, plus_list = SET2_plus, minus_list = SET2_minus, plot_flag = 0, annot = ORF_annot_c)))

their_sense_reads = cbind.data.frame(ORF_annot_c, WT_NC_counts[,1], WT_mRNA_counts[,1], RCO1_counts[,1], DST1_counts[,1], EAF3_counts[,1], SET1_counts[,1], SET2_counts[,1])
colnames(their_sense_reads) = c( colnames(their_sense_reads)[1:6], 'WT_NC', 'WT_mRNA', 'RCO1', 'DST1', 'EAF3', 'SET1', 'SET2')

their_antisense_reads = cbind.data.frame(ORF_annot_c, WT_NC_counts[,2], WT_mRNA_counts[,2], RCO1_counts[,2], DST1_counts[,2], EAF3_counts[,2], SET1_counts[,2], SET2_counts[,2])
colnames(their_antisense_reads) = c( colnames(their_antisense_reads)[1:6], 'WT_NC', 'WT_mRNA', 'RCO1', 'DST1', 'EAF3', 'SET1', 'SET2')

# aggregate counts in ORF regions for our data
WT_counts = t(sapply(1:nrow(ORF_annot), function(x) netseq_plotter(x, plus_list = WT_clean[[1]], minus_list = WT_clean[[2]], plot_flag = 0, annot = ORF_annot)))
WT_pilot_counts = t(sapply(1:nrow(ORF_annot), function(x) netseq_plotter(x, plus_list = WT_pilot_clean[[1]], minus_list = WT_pilot_clean[[2]], plot_flag = 0, annot = ORF_annot)))
CCR4_counts = t(sapply(1:nrow(ORF_annot), function(x) netseq_plotter(x, plus_list = CCR4_clean[[1]], minus_list = CCR4_clean[[2]], plot_flag = 0, annot = ORF_annot)))
DHH1_counts = t(sapply(1:nrow(ORF_annot), function(x) netseq_plotter(x, plus_list = DHH1_clean[[1]], minus_list = DHH1_clean[[2]], plot_flag = 0, annot = ORF_annot)))
LSM1_counts = t(sapply(1:nrow(ORF_annot), function(x) netseq_plotter(x, plus_list = LSM1_clean[[1]], minus_list = LSM1_clean[[2]], plot_flag = 0, annot = ORF_annot)))
RPB4_counts = t(sapply(1:nrow(ORF_annot), function(x) netseq_plotter(x, plus_list = RPB4_clean[[1]], minus_list = RPB4_clean[[2]], plot_flag = 0, annot = ORF_annot)))
XRN1_counts = t(sapply(1:nrow(ORF_annot), function(x) netseq_plotter(x, plus_list = XRN1_clean[[1]], minus_list = XRN1_clean[[2]], plot_flag = 0, annot = ORF_annot)))
XRN1_pilot_counts = t(sapply(1:nrow(ORF_annot), function(x) netseq_plotter(x, plus_list = XRN1_pilot_clean[[1]], minus_list = XRN1_pilot_clean[[2]], plot_flag = 0, annot = ORF_annot)))

our_sense_reads = cbind.data.frame(ORF_annot, WT_counts[,1], WT_pilot_counts[,1], CCR4_counts[,1], DHH1_counts[,1], LSM1_counts[,1], RPB4_counts[,1], 
                                   XRN1_counts[,1], XRN1_pilot_counts[,1])
colnames(our_sense_reads) = c( colnames(our_sense_reads)[1:6], 'WT', 'WT_pilot', 'CCR4', 'DHH1', 'LSM1', 'RPB4', 'XRN1', 'XRN1_pilot')
rownames(our_sense_reads) = ORF_annot$commonName

our_antisense_reads = cbind.data.frame(ORF_annot, WT_counts[,2], WT_pilot_counts[,2], CCR4_counts[,2], DHH1_counts[,2], LSM1_counts[,2], RPB4_counts[,2], 
                                       XRN1_counts[,2], XRN1_pilot_counts[,2])
colnames(our_antisense_reads) = c( colnames(our_antisense_reads)[1:5], 'WT', 'WT_pilot', 'CCR4', 'DHH1', 'LSM1', 'RPB4', 'XRN1', 'XRN1_pilot')

ORF_table = our_sense_reads[,-c(1:6)]
CUT_annot_c = CUT_annot
NUT_annot_c = NUT_annot
SUT_annot_c = SUT_annot
XUT_annot_c = XUT_annot

# Account for shift identified between annotations
CUT_annot$start[CUT_annot$chr == 11 & CUT_annot$start > 201700] = CUT_annot$start[CUT_annot$chr == 11 & CUT_annot$start > 201700]+355
CUT_annot$start[CUT_annot$chr == 11 & CUT_annot$end > 201700] = CUT_annot$end[CUT_annot$chr == 11 & CUT_annot$end > 201700]+355
NUT_annot$start[NUT_annot$chr == 11 & NUT_annot$start > 201700] = NUT_annot$start[NUT_annot$chr == 11 & NUT_annot$start > 201700]+355
NUT_annot$start[NUT_annot$chr == 11 & NUT_annot$end > 201700] = NUT_annot$end[NUT_annot$chr == 11 & NUT_annot$end > 201700]+355
SUT_annot$start[SUT_annot$chr == 11 & SUT_annot$start > 201700] = SUT_annot$start[SUT_annot$chr == 11 & SUT_annot$start > 201700]+355
SUT_annot$start[SUT_annot$chr == 11 & SUT_annot$end > 201700] = SUT_annot$end[SUT_annot$chr == 11 & SUT_annot$end > 201700]+355
XUT_annot$start[XUT_annot$chr == 11 & XUT_annot$start > 201700] = XUT_annot$start[XUT_annot$chr == 11 & XUT_annot$start > 201700]+355
XUT_annot$start[XUT_annot$chr == 11 & XUT_annot$end > 201700] = XUT_annot$end[XUT_annot$chr == 11 & XUT_annot$end > 201700]+355

# Get counts for ncRNA regions
###################################################################################################################################################################################################################################################
CUT_reads_WT = t(sapply(1:nrow(CUT_annot), function(x) interval_reads(start = CUT_annot$start[x], end = CUT_annot$end[x], chrom = CUT_annot$chr[x], strand = CUT_annot$strand[x], plus_list = WT_clean[[1]], minus_list = WT_clean[[2]], plot_flag = 0)))[,1]
CUT_reads_WT_pilot = t(sapply(1:nrow(CUT_annot), function(x) interval_reads(start = CUT_annot$start[x], end = CUT_annot$end[x], chrom = CUT_annot$chr[x], strand = CUT_annot$strand[x], plus_list = WT_pilot_clean[[1]], minus_list = WT_pilot_clean[[2]], plot_flag = 0)))[,1]
CUT_reads_CCR4 = t(sapply(1:nrow(CUT_annot), function(x) interval_reads(start = CUT_annot$start[x], end = CUT_annot$end[x], chrom = CUT_annot$chr[x], strand = CUT_annot$strand[x], plus_list = CCR4_clean[[1]], minus_list = CCR4_clean[[2]], plot_flag = 0)))[,1]
CUT_reads_DHH1 = t(sapply(1:nrow(CUT_annot), function(x) interval_reads(start = CUT_annot$start[x], end = CUT_annot$end[x], chrom = CUT_annot$chr[x], strand = CUT_annot$strand[x], plus_list = DHH1_clean[[1]], minus_list = DHH1_clean[[2]], plot_flag = 0)))[,1]
CUT_reads_LSM1 = t(sapply(1:nrow(CUT_annot), function(x) interval_reads(start = CUT_annot$start[x], end = CUT_annot$end[x], chrom = CUT_annot$chr[x], strand = CUT_annot$strand[x], plus_list = LSM1_clean[[1]], minus_list = LSM1_clean[[2]], plot_flag = 0)))[,1]
CUT_reads_RPB4 = t(sapply(1:nrow(CUT_annot), function(x) interval_reads(start = CUT_annot$start[x], end = CUT_annot$end[x], chrom = CUT_annot$chr[x], strand = CUT_annot$strand[x], plus_list = RPB4_clean[[1]], minus_list = RPB4_clean[[2]], plot_flag = 0)))[,1]
CUT_reads_XRN1 = t(sapply(1:nrow(CUT_annot), function(x) interval_reads(start = CUT_annot$start[x], end = CUT_annot$end[x], chrom = CUT_annot$chr[x], strand = CUT_annot$strand[x], plus_list = XRN1_clean[[1]], minus_list = XRN1_clean[[2]], plot_flag = 0)))[,1]
CUT_reads_XRN1_pilot = t(sapply(1:nrow(CUT_annot), function(x) interval_reads(start = CUT_annot$start[x], end = CUT_annot$end[x], chrom = CUT_annot$chr[x], strand = CUT_annot$strand[x], plus_list = XRN1_pilot_clean[[1]], minus_list = XRN1_pilot_clean[[2]], plot_flag = 0)))[,1]
CUT_reads = cbind(CUT_reads_WT, CUT_reads_WT_pilot, CUT_reads_CCR4, CUT_reads_DHH1, CUT_reads_LSM1, CUT_reads_RPB4, CUT_reads_SFP1, CUT_reads_XRN1, CUT_reads_XRN1_pilot)
colnames(CUT_reads) = colnames(our_sense_reads)[-c(1:5)]
###################################################################################################################################################################################################################################################
NUT_reads_WT = t(sapply(1:nrow(NUT_annot), function(x) interval_reads(start = NUT_annot$start[x], end = NUT_annot$end[x], chrom = NUT_annot$chr[x], strand = NUT_annot$strand[x], plus_list = WT_clean[[1]], minus_list = WT_clean[[2]], plot_flag = 0)))[,1]
NUT_reads_WT_pilot = t(sapply(1:nrow(NUT_annot), function(x) interval_reads(start = NUT_annot$start[x], end = NUT_annot$end[x], chrom = NUT_annot$chr[x], strand = NUT_annot$strand[x], plus_list = WT_pilot_clean[[1]], minus_list = WT_pilot_clean[[2]], plot_flag = 0)))[,1]
NUT_reads_CCR4 = t(sapply(1:nrow(NUT_annot), function(x) interval_reads(start = NUT_annot$start[x], end = NUT_annot$end[x], chrom = NUT_annot$chr[x], strand = NUT_annot$strand[x], plus_list = CCR4_clean[[1]], minus_list = CCR4_clean[[2]], plot_flag = 0)))[,1]
NUT_reads_DHH1 = t(sapply(1:nrow(NUT_annot), function(x) interval_reads(start = NUT_annot$start[x], end = NUT_annot$end[x], chrom = NUT_annot$chr[x], strand = NUT_annot$strand[x], plus_list = DHH1_clean[[1]], minus_list = DHH1_clean[[2]], plot_flag = 0)))[,1]
NUT_reads_LSM1 = t(sapply(1:nrow(NUT_annot), function(x) interval_reads(start = NUT_annot$start[x], end = NUT_annot$end[x], chrom = NUT_annot$chr[x], strand = NUT_annot$strand[x], plus_list = LSM1_clean[[1]], minus_list = LSM1_clean[[2]], plot_flag = 0)))[,1]
NUT_reads_RPB4 = t(sapply(1:nrow(NUT_annot), function(x) interval_reads(start = NUT_annot$start[x], end = NUT_annot$end[x], chrom = NUT_annot$chr[x], strand = NUT_annot$strand[x], plus_list = RPB4_clean[[1]], minus_list = RPB4_clean[[2]], plot_flag = 0)))[,1]
NUT_reads_XRN1 = t(sapply(1:nrow(NUT_annot), function(x) interval_reads(start = NUT_annot$start[x], end = NUT_annot$end[x], chrom = NUT_annot$chr[x], strand = NUT_annot$strand[x], plus_list = XRN1_clean[[1]], minus_list = XRN1_clean[[2]], plot_flag = 0)))[,1]
NUT_reads_XRN1_pilot = t(sapply(1:nrow(NUT_annot), function(x) interval_reads(start = NUT_annot$start[x], end = NUT_annot$end[x], chrom = NUT_annot$chr[x], strand = NUT_annot$strand[x], plus_list = XRN1_pilot_clean[[1]], minus_list = XRN1_pilot_clean[[2]], plot_flag = 0)))[,1]
NUT_reads = cbind(NUT_reads_WT, NUT_reads_WT_pilot, NUT_reads_CCR4, NUT_reads_DHH1, NUT_reads_LSM1, NUT_reads_RPB4, NUT_reads_SFP1, NUT_reads_XRN1, NUT_reads_XRN1_pilot)
colnames(NUT_reads) = colnames(our_sense_reads)[-c(1:5)]
###################################################################################################################################################################################################################################################
SUT_reads_WT = t(sapply(1:nrow(SUT_annot), function(x) interval_reads(start = SUT_annot$start[x], end = SUT_annot$end[x], chrom = SUT_annot$chr[x], strand = SUT_annot$strand[x], plus_list = WT_clean[[1]], minus_list = WT_clean[[2]], plot_flag = 0)))[,1]
SUT_reads_WT_pilot = t(sapply(1:nrow(SUT_annot), function(x) interval_reads(start = SUT_annot$start[x], end = SUT_annot$end[x], chrom = SUT_annot$chr[x], strand = SUT_annot$strand[x], plus_list = WT_pilot_clean[[1]], minus_list = WT_pilot_clean[[2]], plot_flag = 0)))[,1]
SUT_reads_CCR4 = t(sapply(1:nrow(SUT_annot), function(x) interval_reads(start = SUT_annot$start[x], end = SUT_annot$end[x], chrom = SUT_annot$chr[x], strand = SUT_annot$strand[x], plus_list = CCR4_clean[[1]], minus_list = CCR4_clean[[2]], plot_flag = 0)))[,1]
SUT_reads_DHH1 = t(sapply(1:nrow(SUT_annot), function(x) interval_reads(start = SUT_annot$start[x], end = SUT_annot$end[x], chrom = SUT_annot$chr[x], strand = SUT_annot$strand[x], plus_list = DHH1_clean[[1]], minus_list = DHH1_clean[[2]], plot_flag = 0)))[,1]
SUT_reads_LSM1 = t(sapply(1:nrow(SUT_annot), function(x) interval_reads(start = SUT_annot$start[x], end = SUT_annot$end[x], chrom = SUT_annot$chr[x], strand = SUT_annot$strand[x], plus_list = LSM1_clean[[1]], minus_list = LSM1_clean[[2]], plot_flag = 0)))[,1]
SUT_reads_RPB4 = t(sapply(1:nrow(SUT_annot), function(x) interval_reads(start = SUT_annot$start[x], end = SUT_annot$end[x], chrom = SUT_annot$chr[x], strand = SUT_annot$strand[x], plus_list = RPB4_clean[[1]], minus_list = RPB4_clean[[2]], plot_flag = 0)))[,1]
SUT_reads_XRN1 = t(sapply(1:nrow(SUT_annot), function(x) interval_reads(start = SUT_annot$start[x], end = SUT_annot$end[x], chrom = SUT_annot$chr[x], strand = SUT_annot$strand[x], plus_list = XRN1_clean[[1]], minus_list = XRN1_clean[[2]], plot_flag = 0)))[,1]
SUT_reads_XRN1_pilot = t(sapply(1:nrow(SUT_annot), function(x) interval_reads(start = SUT_annot$start[x], end = SUT_annot$end[x], chrom = SUT_annot$chr[x], strand = SUT_annot$strand[x], plus_list = XRN1_pilot_clean[[1]], minus_list = XRN1_pilot_clean[[2]], plot_flag = 0)))[,1]
SUT_reads = cbind(SUT_reads_WT, SUT_reads_WT_pilot, SUT_reads_CCR4, SUT_reads_DHH1, SUT_reads_LSM1, SUT_reads_RPB4, SUT_reads_SFP1, SUT_reads_XRN1, SUT_reads_XRN1_pilot)
colnames(SUT_reads) = colnames(our_sense_reads)[-c(1:5)]
###################################################################################################################################################################################################################################################
XUT_reads_WT = t(sapply(1:nrow(XUT_annot), function(x) interval_reads(start = XUT_annot$start[x], end = XUT_annot$end[x], chrom = XUT_annot$chr[x], strand = XUT_annot$strand[x], plus_list = WT_clean[[1]], minus_list = WT_clean[[2]], plot_flag = 0)))[,1]
XUT_reads_WT_pilot = t(sapply(1:nrow(XUT_annot), function(x) interval_reads(start = XUT_annot$start[x], end = XUT_annot$end[x], chrom = XUT_annot$chr[x], strand = XUT_annot$strand[x], plus_list = WT_pilot_clean[[1]], minus_list = WT_pilot_clean[[2]], plot_flag = 0)))[,1]
XUT_reads_CCR4 = t(sapply(1:nrow(XUT_annot), function(x) interval_reads(start = XUT_annot$start[x], end = XUT_annot$end[x], chrom = XUT_annot$chr[x], strand = XUT_annot$strand[x], plus_list = CCR4_clean[[1]], minus_list = CCR4_clean[[2]], plot_flag = 0)))[,1]
XUT_reads_DHH1 = t(sapply(1:nrow(XUT_annot), function(x) interval_reads(start = XUT_annot$start[x], end = XUT_annot$end[x], chrom = XUT_annot$chr[x], strand = XUT_annot$strand[x], plus_list = DHH1_clean[[1]], minus_list = DHH1_clean[[2]], plot_flag = 0)))[,1]
XUT_reads_LSM1 = t(sapply(1:nrow(XUT_annot), function(x) interval_reads(start = XUT_annot$start[x], end = XUT_annot$end[x], chrom = XUT_annot$chr[x], strand = XUT_annot$strand[x], plus_list = LSM1_clean[[1]], minus_list = LSM1_clean[[2]], plot_flag = 0)))[,1]
XUT_reads_RPB4 = t(sapply(1:nrow(XUT_annot), function(x) interval_reads(start = XUT_annot$start[x], end = XUT_annot$end[x], chrom = XUT_annot$chr[x], strand = XUT_annot$strand[x], plus_list = RPB4_clean[[1]], minus_list = RPB4_clean[[2]], plot_flag = 0)))[,1]
XUT_reads_XRN1 = t(sapply(1:nrow(XUT_annot), function(x) interval_reads(start = XUT_annot$start[x], end = XUT_annot$end[x], chrom = XUT_annot$chr[x], strand = XUT_annot$strand[x], plus_list = XRN1_clean[[1]], minus_list = XRN1_clean[[2]], plot_flag = 0)))[,1]
XUT_reads_XRN1_pilot = t(sapply(1:nrow(XUT_annot), function(x) interval_reads(start = XUT_annot$start[x], end = XUT_annot$end[x], chrom = XUT_annot$chr[x], strand = XUT_annot$strand[x], plus_list = XRN1_pilot_clean[[1]], minus_list = XRN1_pilot_clean[[2]], plot_flag = 0)))[,1]
XUT_reads = cbind(XUT_reads_WT, XUT_reads_WT_pilot, XUT_reads_CCR4, XUT_reads_DHH1, XUT_reads_LSM1, XUT_reads_RPB4, XUT_reads_SFP1, XUT_reads_XRN1, XUT_reads_XRN1_pilot)
colnames(XUT_reads) = colnames(our_sense_reads)[-c(1:5)]
###################################################################################################################################################################################################################################################
CUT_reads_WT_NC = t(sapply(1:nrow(CUT_annot_c), function(x) interval_reads(start = CUT_annot_c$start[x], end = CUT_annot_c$end[x], chrom = CUT_annot_c$chr[x], strand = CUT_annot_c$strand[x], plus_list = WT_NC_plus, minus_list = WT_NC_minus, plot_flag = 0)))[,1]
CUT_reads_WT_mRNA = t(sapply(1:nrow(CUT_annot_c), function(x) interval_reads(start = CUT_annot_c$start[x], end = CUT_annot_c$end[x], chrom = CUT_annot_c$chr[x], strand = CUT_annot_c$strand[x], plus_list = WT_mRNA_plus, minus_list = WT_mRNA_minus, plot_flag = 0)))[,1]
CUT_reads_RCO1 = t(sapply(1:nrow(CUT_annot_c), function(x) interval_reads(start = CUT_annot_c$start[x], end = CUT_annot_c$end[x], chrom = CUT_annot_c$chr[x], strand = CUT_annot_c$strand[x], plus_list = RCO1_plus, minus_list = RCO1_minus, plot_flag = 0)))[,1]
CUT_reads_DST1 = t(sapply(1:nrow(CUT_annot_c), function(x) interval_reads(start = CUT_annot_c$start[x], end = CUT_annot_c$end[x], chrom = CUT_annot_c$chr[x], strand = CUT_annot_c$strand[x], plus_list = DST1_plus, minus_list = DST1_minus, plot_flag = 0)))[,1]
CUT_reads_EAF3 = t(sapply(1:nrow(CUT_annot_c), function(x) interval_reads(start = CUT_annot_c$start[x], end = CUT_annot_c$end[x], chrom = CUT_annot_c$chr[x], strand = CUT_annot_c$strand[x], plus_list = EAF3_plus, minus_list = EAF3_minus, plot_flag = 0)))[,1]
CUT_reads_SET1 = t(sapply(1:nrow(CUT_annot_c), function(x) interval_reads(start = CUT_annot_c$start[x], end = CUT_annot_c$end[x], chrom = CUT_annot_c$chr[x], strand = CUT_annot_c$strand[x], plus_list = SET1_plus, minus_list = SET1_minus, plot_flag = 0)))[,1]
CUT_reads_SET2 = t(sapply(1:nrow(CUT_annot_c), function(x) interval_reads(start = CUT_annot_c$start[x], end = CUT_annot_c$end[x], chrom = CUT_annot_c$chr[x], strand = CUT_annot_c$strand[x], plus_list = SET2_plus, minus_list = SET2_minus, plot_flag = 0)))[,1]
CUT_reads_c = cbind(CUT_reads_WT_NC, CUT_reads_WT_mRNA, CUT_reads_RCO1, CUT_reads_DST1, CUT_reads_EAF3, CUT_reads_SET1, CUT_reads_SET2)
colnames(CUT_reads_c) = colnames(their_sense_reads)[-c(1:6)]
###################################################################################################################################################################################################################################################
NUT_reads_WT_NC = t(sapply(1:nrow(NUT_annot_c), function(x) interval_reads(start = NUT_annot_c$start[x], end = NUT_annot_c$end[x], chrom = NUT_annot_c$chr[x], strand = NUT_annot_c$strand[x], plus_list = WT_NC_plus, minus_list = WT_NC_minus, plot_flag = 0)))[,1]
NUT_reads_WT_mRNA = t(sapply(1:nrow(NUT_annot_c), function(x) interval_reads(start = NUT_annot_c$start[x], end = NUT_annot_c$end[x], chrom = NUT_annot_c$chr[x], strand = NUT_annot_c$strand[x], plus_list = WT_mRNA_plus, minus_list = WT_mRNA_minus, plot_flag = 0)))[,1]
NUT_reads_RCO1 = t(sapply(1:nrow(NUT_annot_c), function(x) interval_reads(start = NUT_annot_c$start[x], end = NUT_annot_c$end[x], chrom = NUT_annot_c$chr[x], strand = NUT_annot_c$strand[x], plus_list = RCO1_plus, minus_list = RCO1_minus, plot_flag = 0)))[,1]
NUT_reads_DST1 = t(sapply(1:nrow(NUT_annot_c), function(x) interval_reads(start = NUT_annot_c$start[x], end = NUT_annot_c$end[x], chrom = NUT_annot_c$chr[x], strand = NUT_annot_c$strand[x], plus_list = DST1_plus, minus_list = DST1_minus, plot_flag = 0)))[,1]
NUT_reads_EAF3 = t(sapply(1:nrow(NUT_annot_c), function(x) interval_reads(start = NUT_annot_c$start[x], end = NUT_annot_c$end[x], chrom = NUT_annot_c$chr[x], strand = NUT_annot_c$strand[x], plus_list = EAF3_plus, minus_list = EAF3_minus, plot_flag = 0)))[,1]
NUT_reads_SET1 = t(sapply(1:nrow(NUT_annot_c), function(x) interval_reads(start = NUT_annot_c$start[x], end = NUT_annot_c$end[x], chrom = NUT_annot_c$chr[x], strand = NUT_annot_c$strand[x], plus_list = SET1_plus, minus_list = SET1_minus, plot_flag = 0)))[,1]
NUT_reads_SET2 = t(sapply(1:nrow(NUT_annot_c), function(x) interval_reads(start = NUT_annot_c$start[x], end = NUT_annot_c$end[x], chrom = NUT_annot_c$chr[x], strand = NUT_annot_c$strand[x], plus_list = SET2_plus, minus_list = SET2_minus, plot_flag = 0)))[,1]
NUT_reads_c = cbind(NUT_reads_WT_NC, NUT_reads_WT_mRNA, NUT_reads_RCO1, NUT_reads_DST1, NUT_reads_EAF3, NUT_reads_SET1, NUT_reads_SET2)
colnames(NUT_reads_c) = colnames(their_sense_reads)[-c(1:6)]
###################################################################################################################################################################################################################################################
SUT_reads_WT_NC = t(sapply(1:nrow(SUT_annot_c), function(x) interval_reads(start = SUT_annot_c$start[x], end = SUT_annot_c$end[x], chrom = SUT_annot_c$chr[x], strand = SUT_annot_c$strand[x], plus_list = WT_NC_plus, minus_list = WT_NC_minus, plot_flag = 0)))[,1]
SUT_reads_WT_mRNA = t(sapply(1:nrow(SUT_annot_c), function(x) interval_reads(start = SUT_annot_c$start[x], end = SUT_annot_c$end[x], chrom = SUT_annot_c$chr[x], strand = SUT_annot_c$strand[x], plus_list = WT_mRNA_plus, minus_list = WT_mRNA_minus, plot_flag = 0)))[,1]
SUT_reads_RCO1 = t(sapply(1:nrow(SUT_annot_c), function(x) interval_reads(start = SUT_annot_c$start[x], end = SUT_annot_c$end[x], chrom = SUT_annot_c$chr[x], strand = SUT_annot_c$strand[x], plus_list = RCO1_plus, minus_list = RCO1_minus, plot_flag = 0)))[,1]
SUT_reads_DST1 = t(sapply(1:nrow(SUT_annot_c), function(x) interval_reads(start = SUT_annot_c$start[x], end = SUT_annot_c$end[x], chrom = SUT_annot_c$chr[x], strand = SUT_annot_c$strand[x], plus_list = DST1_plus, minus_list = DST1_minus, plot_flag = 0)))[,1]
SUT_reads_EAF3 = t(sapply(1:nrow(SUT_annot_c), function(x) interval_reads(start = SUT_annot_c$start[x], end = SUT_annot_c$end[x], chrom = SUT_annot_c$chr[x], strand = SUT_annot_c$strand[x], plus_list = EAF3_plus, minus_list = EAF3_minus, plot_flag = 0)))[,1]
SUT_reads_SET1 = t(sapply(1:nrow(SUT_annot_c), function(x) interval_reads(start = SUT_annot_c$start[x], end = SUT_annot_c$end[x], chrom = SUT_annot_c$chr[x], strand = SUT_annot_c$strand[x], plus_list = SET1_plus, minus_list = SET1_minus, plot_flag = 0)))[,1]
SUT_reads_SET2 = t(sapply(1:nrow(SUT_annot_c), function(x) interval_reads(start = SUT_annot_c$start[x], end = SUT_annot_c$end[x], chrom = SUT_annot_c$chr[x], strand = SUT_annot_c$strand[x], plus_list = SET2_plus, minus_list = SET2_minus, plot_flag = 0)))[,1]
SUT_reads_c = cbind(SUT_reads_WT_NC, SUT_reads_WT_mRNA, SUT_reads_RCO1, SUT_reads_DST1, SUT_reads_EAF3, SUT_reads_SET1, SUT_reads_SET2)
colnames(SUT_reads_c) = colnames(their_sense_reads)[-c(1:6)]
###################################################################################################################################################################################################################################################
XUT_reads_WT_NC = t(sapply(1:nrow(XUT_annot_c), function(x) interval_reads(start = XUT_annot_c$start[x], end = XUT_annot_c$end[x], chrom = XUT_annot_c$chr[x], strand = XUT_annot_c$strand[x], plus_list = WT_NC_plus, minus_list = WT_NC_minus, plot_flag = 0)))[,1]
XUT_reads_WT_mRNA = t(sapply(1:nrow(XUT_annot_c), function(x) interval_reads(start = XUT_annot_c$start[x], end = XUT_annot_c$end[x], chrom = XUT_annot_c$chr[x], strand = XUT_annot_c$strand[x], plus_list = WT_mRNA_plus, minus_list = WT_mRNA_minus, plot_flag = 0)))[,1]
XUT_reads_RCO1 = t(sapply(1:nrow(XUT_annot_c), function(x) interval_reads(start = XUT_annot_c$start[x], end = XUT_annot_c$end[x], chrom = XUT_annot_c$chr[x], strand = XUT_annot_c$strand[x], plus_list = RCO1_plus, minus_list = RCO1_minus, plot_flag = 0)))[,1]
XUT_reads_DST1 = t(sapply(1:nrow(XUT_annot_c), function(x) interval_reads(start = XUT_annot_c$start[x], end = XUT_annot_c$end[x], chrom = XUT_annot_c$chr[x], strand = XUT_annot_c$strand[x], plus_list = DST1_plus, minus_list = DST1_minus, plot_flag = 0)))[,1]
XUT_reads_EAF3 = t(sapply(1:nrow(XUT_annot_c), function(x) interval_reads(start = XUT_annot_c$start[x], end = XUT_annot_c$end[x], chrom = XUT_annot_c$chr[x], strand = XUT_annot_c$strand[x], plus_list = EAF3_plus, minus_list = EAF3_minus, plot_flag = 0)))[,1]
XUT_reads_SET1 = t(sapply(1:nrow(XUT_annot_c), function(x) interval_reads(start = XUT_annot_c$start[x], end = XUT_annot_c$end[x], chrom = XUT_annot_c$chr[x], strand = XUT_annot_c$strand[x], plus_list = SET1_plus, minus_list = SET1_minus, plot_flag = 0)))[,1]
XUT_reads_SET2 = t(sapply(1:nrow(XUT_annot_c), function(x) interval_reads(start = XUT_annot_c$start[x], end = XUT_annot_c$end[x], chrom = XUT_annot_c$chr[x], strand = XUT_annot_c$strand[x], plus_list = SET2_plus, minus_list = SET2_minus, plot_flag = 0)))[,1]
XUT_reads_c = cbind(XUT_reads_WT_NC, XUT_reads_WT_mRNA, XUT_reads_RCO1, XUT_reads_DST1, XUT_reads_EAF3, XUT_reads_SET1, XUT_reads_SET2)
colnames(XUT_reads_c) = colnames(their_sense_reads)[-c(1:6)]

# load prior Xrn1 sensitivity ratings
xrn1_synth_ratings = read.csv('~/Documents/R_Projects/NETSEQ/data/xrn1_synth_ratings.csv')

