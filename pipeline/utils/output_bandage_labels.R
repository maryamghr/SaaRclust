log <- file(snakemake@log[[1]], open='wt')
sink(file=log, type='message')
sink(file=log, type='output')

library(data.table)
library(seqinr)
library(Rsamtools)
library(reshape2)
#setwd('aligns_k15_w1_f0.1_z500/SaaRclust_results_HG00733/Clusters/split/hifiasm/')


clust_pair <- snakemake@wildcards[['clust_pair']]
clust.to.chrom <- fread(snakemake@input[["clust_to_chrom_file"]])
ss.clust <- fread(snakemake@input[["ss_clust"]], col.names=c("SSname_lib", "SSclust"))
ss.phase <- fread(snakemake@input[["ss_phase"]], col.names=c("SSlib", "SSclust", "SShaplo"))
unitigs.fasta <- read.fasta(snakemake@input[["unitigs_fasta"]])
unitig.ss.map <- fread(snakemake@input[["valid_maps"]])
unitigs.phase <- fread(snakemake@input[["phased_unitigs"]])
unitigs.bam <- scanBam(BamFile(snakemake@input[["unitigs_bam"]]), param=ScanBamParam(what=c("rname", "qname"), tag="HP", flag=scanBamFlag(isSupplementaryAlignment=FALSE)))
ccs.bam <- scanBam(BamFile(snakemake@input[["ccs_bam"]]), param=ScanBamParam(what=c("rname", "qname"), tag="HP", flag=scanBamFlag(isSupplementaryAlignment=FALSE)))
hifiasm.gfa <- fread(snakemake@input[["hifiasm_gfa"]], fill=T)

################################
##### R Consule:
#clust_pair <- 'V25_V32'
#clust.to.chrom <- fread('../../../clust_partners.txt')
#ss.clust <- fread(paste0('../ss_reads/cluster', clust_pair, '.data'), header=F, col.names=c("SSname_lib", "SSclust"))
#ss.phase <- fread(paste0('phased_strand_states/haplo_strand_states_', clust_pair, '_r.data'), col.names=c("SSlib", "SSclust", "SShaplo"))
#unitigs.fasta <- read.fasta(paste0(clust_pair,'.r_utg.fa'))
#unitig.ss.map <- fread(paste0('bubble_ss_map/valid_',clust_pair,'_r_utg_maximal_uniqe_exact_match.data'))
#unitigs.phase <- fread(paste0('phased_bubbles/iteration1_',clust_pair,'_r.data'))
#unitigs.bam <- scanBam(BamFile(paste0('bubble_ref_aln/',clust_pair,'.r.haplotagged.bam')), param=ScanBamParam(what=c("rname", "qname"), tag="HP", flag=scanBamFlag(isSupplementaryAlignment=FALSE)))
#ccs.bam <- scanBam(BamFile(paste0('../HG00733.',clust_pair,'.haplotagged.bam')), param=ScanBamParam(what=c("rname", "qname"), tag="HP", flag=scanBamFlag(isSupplementaryAlignment=FALSE)))
#hifiasm.gfa <- fread('V25_V32.r_utg.gfa', fill=T)
################################

valid.chr <- paste0('chr', c(1:22, "X"))
clust1 <- strsplit(clust_pair, '_')[[1]][1]
chr.name <- clust.to.chrom[clust.forward==clust1, original.chrom] #'chr22'
chr.name <- strsplit(chr.name, '_')[[1]][1]

colnames(unitigs.phase)[1] <- 'unitig_name'

unitigs.gr <- data.table(ref_name=unitigs.bam[[1]]$rname, unitig_name=unitigs.bam[[1]]$qname, unitig_gr_haplo=as.character(unlist(unitigs.bam[[1]]$tag)))
ccs.gr <- data.table(ref_name=ccs.bam[[1]]$rname, ccs_name=ccs.bam[[1]]$qname, ccs_gr_haplo=as.character(unlist(ccs.bam[[1]]$tag)))

ss.clust[, SSname:=strsplit(SSname_lib, '_')[[1]][1], by=SSname_lib]
ss.clust[, SSlib:=strsplit(SSname_lib, '_')[[1]][2], by=SSname_lib]
ss.clust[, SSname_lib:=NULL]	
ss.phase <- merge(ss.clust, ss.phase, by=c("SSlib", "SSclust"))

unitigs.gr[unitig_gr_haplo==1, unitig_gr_haplo:="H1"]
unitigs.gr[unitig_gr_haplo==2, unitig_gr_haplo:="H2"]
ccs.gr[ccs_gr_haplo==1, ccs_gr_haplo:="H1"]
ccs.gr[ccs_gr_haplo==2, ccs_gr_haplo:="H2"]

# filter out not haplotagged ccs reads
#ccs.gr <- ccs.gr[!is.na(ccs_gr_haplo)]

hifiasm.gfa <- hifiasm.gfa[V1=="A", .(V2, V5)]
colnames(hifiasm.gfa) <- c("unitig_name", "ccs_name")
#unitigs.gr[, `:=`(ref_haplo=paste0(ref_name, "_", unitig_gr_haplo), ref_name=NULL, unitig_gr_haplo=NULL)]

unitig.lens <- sapply(unitigs.fasta, length)
unitig.lens <- data.table(unitig_name=names(unitig.lens), len=unitig.lens)

unitig.ss.cov <- merge(unitig.ss.map, ss.phase[, .(SSname, SShaplo)])
unitig.ss.cov[, unitig_haplo_ss_cov:=.N, by=.(unitig_name, SShaplo)]
unitig.ss.cov <- unitig.ss.cov[, head(.SD, 1), by=.(unitig_name, SShaplo, unitig_haplo_ss_cov)]
unitig.ss.cov <- data.table::dcast(unitig.ss.cov, unitig_name~SShaplo, value.var="unitig_haplo_ss_cov", fill=0)
colnames(unitig.ss.cov) <- c("unitig_name", "ss_cov_h1", "ss_cov_h2")
unitig.ss.cov[, ss_haplo_cov:=paste0("(", ss_cov_h1, ",", ss_cov_h2, ")"), by=unitig_name]

unitig.ccs.cov <- merge(hifiasm.gfa, ccs.gr[, .(ccs_name, ccs_gr_haplo)])
unitig.ccs.cov[, unitig_haplo_ccs_cov:=.N, by=.(unitig_name, ccs_gr_haplo)]
unitig.ccs.cov <- unitig.ccs.cov[, head(.SD, 1), by=.(unitig_name, ccs_gr_haplo, unitig_haplo_ccs_cov)]
unitig.ccs.cov <- data.table::dcast(unitig.ccs.cov, unitig_name~ccs_gr_haplo, value.var="unitig_haplo_ccs_cov", fill=0)
colnames(unitig.ccs.cov) <- c("unitig_name", "ccs_cov_untagged", "ccs_cov_h1", "ccs_cov_h2")
unitig.ccs.cov[, ccs_haplo_cov:=paste0("(", ccs_cov_h1, ",", ccs_cov_h2, ",", ccs_cov_untagged, ")"), by=unitig_name]

unitigs <- merge(unitig.lens, unitig.ss.cov[, .(unitig_name, ss_haplo_cov)], all=TRUE, by='unitig_name')
unitigs <- merge(unitigs, unitig.ccs.cov[, .(unitig_name, ccs_haplo_cov)], all=TRUE, by='unitig_name')
unitigs[is.na(ss_haplo_cov), ss_haplo_cov:="(0,0)"]
unitigs[is.na(ccs_haplo_cov), ss_haplo_cov:="(0,0)"]

unitigs <- merge(unitigs, unitigs.phase[, .(unitig_name, haplotype)], all=TRUE, by='unitig_name')
unitigs <- merge(unitigs, unitigs.gr, by="unitig_name")

#norm.factor <- median(unitigs[, unitig_ss_cov/len])
#unitigs[, CN:=(unitig_ss_cov/len)/norm.factor]
#unitigs[, CN:=round(CN, digits=1)]

# set the Colors for unitig nodes in Bandage
unitigs[, Color:="#FFA07A"]
unitigs[haplotype=="H2", Color:="#87CEEB"]
#unitigs[is.na(haplotype), Color:="gray"]
unitigs[haplotype=="None", Color:="gray"]
unitigs[ref_name!=chr.name, Color:="green"]

fwrite(unitigs[, .(unitig_name, unitig_gr_haplo, Color)], file=snakemake@output[["gr_haplo"]], row.names=F, sep='\t', quote=FALSE)
fwrite(unitigs[, .(unitig_name, ss_haplo_cov, Color)], file=snakemake@output[["ss_haplo_cov"]], row.names=F, sep='\t', quote=FALSE)
#fwrite(unitigs[, .(unitig_name, unitig_gr_haplo, Color)], file=paste0('Bandage_sessions/',clust_pair,'_r_gr_ref_haplo.csv'), row.names=F, sep='\t', quote=FALSE)
#fwrite(unitigs[, .(unitig_name, ss_haplo_cov, Color)], file=paste0('Bandage_sessions/',clust_pair,'_ss_haplo_cov.csv'), row.names=F, sep='\t', quote=FALSE)
#fwrite(unitigs[, .(unitig_name, ccs_haplo_cov, Color)], file=paste0('Bandage_sessions/',clust_pair,'_ccs_haplo_cov.csv'), row.names=F, sep='\t', quote=FALSE)
