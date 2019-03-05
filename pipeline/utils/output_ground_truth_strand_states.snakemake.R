log <- file(snakemake@log[[1]], open='wt')
sink(file=log, type='message')
sink(file=log, type='output')

library(data.table)

ss_count_watson <- fread(snakemake@input[["ss_count_watson"]])
colnames(ss_count_watson) <- c("W", "chrom")

ss_count_crick <- fread(snakemake@input[["ss_count_crick"]])
colnames(ss_count_crick) <- c("C", "chrom")

ss_count <- merge(ss_count_watson, ss_count_crick, by="chrom", all=T)

ss_count_watson_h1 <- fread(snakemake@input[["ss_haplotagged_count_watson_h1"]])
colnames(ss_count_watson_h1) <- c("W_h1", "chrom")

ss_count <- merge(ss_count, ss_count_watson_h1, by="chrom", all=T)

ss_count_watson_h2 <- fread(snakemake@input[["ss_haplotagged_count_watson_h2"]])
colnames(ss_count_watson_h2) <- c("W_h2", "chrom")

ss_count <- merge(ss_count, ss_count_watson_h2, by="chrom", all=T)

ss_count_crick_h1 <- fread(snakemake@input[["ss_haplotagged_count_crick_h1"]])
colnames(ss_count_crick_h1) <- c("C_h1", "chrom")

ss_count <- merge(ss_count, ss_count_crick_h1, by="chrom", all=T)

ss_count_crick_h2 <- fread(snakemake@input[["ss_haplotagged_count_crick_h2"]])
colnames(ss_count_crick_h2) <- c("C_h2", "chrom")

ss_count <- merge(ss_count, ss_count_crick_h2, by="chrom", all=T)
ss_count[is.na(ss_count)] <- 0

### calling W and C (mostly covered) haplotypes
ss_count[, `:=`(W_h=which.max(c(W_h1, W_h2)), C_h=which.max(c(C_h1, C_h2))), by=chrom]

### binomial tests for wc cells
# note that the last two tests are very strong and enough to distinguish WC cells
ss_count[, wc.p.value:=binom.test(x=c(W, C), p=0.5, alternative="two.sided")$p.value, by=chrom]
ss_count[, w.haplo.p.value:=binom.test(x=c(min(W_h1, W_h2), max(W_h1, W_h2)), p=as.numeric(snakemake@params[["background_rate"]]), alternative="greater")$p.value, by=chrom]
ss_count[, c.haplo.p.value:=binom.test(x=c(min(C_h1, C_h2), max(C_h1, C_h2)), p=as.numeric(snakemake@params[["background_rate"]]), alternative="greater")$p.value, by=chrom]

ss_count_w = ss_count[W_h!=C_h & w.haplo.p.value>0.01 & c.haplo.p.value>0.01, .(original.cluster=chrom, haplotype=W_h)]
ss_count_c = ss_count[W_h!=C_h & w.haplo.p.value>0.01 & c.haplo.p.value>0.01, .(original.cluster=chrom, haplotype=C_h)]
ss_count_w[, original.cluster:=paste0(original.cluster, "_16")]
ss_count_c[, original.cluster:=paste0(original.cluster, "_0")]

strand_states <- rbind(ss_count_w, ss_count_c)
setkey(strand_states, original.cluster)

cluster_names <- fread(snakemake@input[["clust_to_chrom"]])
colnames(cluster_names) <- c("original.cluster", "inferred.cluster", "inferred.cluster.pair")
strand_states <- merge(strand_states, cluster_names, by="original.cluster")

# reorder columns
strand_states <- strand_states[, .(inferred.cluster, haplotype, original.cluster)]

fwrite(strand_states, file=snakemake@output[[1]], sep="\t")
