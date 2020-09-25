library(data.table)

clusters <- c('V25', 'V32')

phase.V25_V32 <- fread(paste0('aligns_k15_w1_f0.1_z500/SaaRclust_results_HG00733/Clusters/split/hifiasm/phased_strand_states/haplo_strand_states_', clusters[1], '_', clusters[2], '_r.data'))

gr.phase.files <- list.files('ground_truth_strand_states', pattern='haplo_strand_states_k15_w1_f0.1_z500_HG00733.data$', full.names=TRUE)

gr.phase <- data.table()
for (f in gr.phase.files){
 filename <- basename(f)
 lib.name <- strsplit(filename, '_')[[1]][1]
 dt <- fread(f)
 dt[, lib:=lib.name]
 gr.phase <- rbind(gr.phase, dt)
}

merged.phase <- merge(gr.phase, phase.V25_V32, by=c('lib', 'cluster'))
merged.phase[haplotype.x==haplotype.y]
