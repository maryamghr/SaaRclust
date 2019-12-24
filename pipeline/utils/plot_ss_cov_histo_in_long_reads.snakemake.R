log <- file(snakemake@log[[1]], open='wt')
sink(file=log, type='message')
sink(file=log, type='output')

library(data.table)
library(ggplot2)
library(gridExtra)



#' This function reads the coverage of ss reads/libs in long reads and return the coverage quantiles
#' 
#' @param coverage_files a set of files containing the coverage of ss reads/libs in long reads. Each file has two columns: ss coverage and long read name, respectively
#' @param probs An \code{numeric} vecotr of probabilities with values in range [0,1] for which the quantiles are to be computed
#' @author Maryam Ghareghani
#' @export


getCoverageQuantiles <- function(coverage_files, probs)
{
  cov <- lapply(coverage_files, function(x) fread(x, col.names=c('coverage', 'long_read_name')))
  cov.all <- Reduce(rbind, cov)
  return(quantile(cov.all[, coverage], probs=probs))
}

cov <- lapply(snakemake@input[["ss_cov"]], function(x) fread(x, col.names=c('ss_coverage', 'long_read_name')))
cov.all <- Reduce(rbind, cov)
print(cov.all)

lib.cov <- lapply(snakemake@input[["ss_lib_cov"]], function(x) fread(x, col.names=c('ss_lib_coverage', 'long_read_name')))
lib.cov.all <- Reduce(rbind, lib.cov)
print(lib.cov.all)

q_ss_cov <- quantile(cov.all[, ss_coverage], probs=c(0.5, 0.75, 0.85, 0.95))
print('quantile of coverage of ss reads')
print(q_ss_cov)

q_ss_lib_cov <- quantile(lib.cov.all[, ss_lib_coverage], probs=c(0.5, 0.75, 0.85, 0.95))
print('quantile of coverage of ss libs')
print(q_ss_lib_cov)

cov.all.removetail.05 <- cov.all[ss_coverage < q_ss_cov[4]]
lib.cov.all.removetail.05 <- lib.cov.all[ss_lib_coverage < q_ss_lib_cov[4]]


cov.hist.removetail.05 <- ggplot(cov.all.removetail.05, aes(ss_coverage))+geom_histogram(binwidth=1)+ggtitle('coverage of ss reads in long reads removing 5% quantile tail')
lib.cov.hist <- ggplot(lib.cov.all, aes(ss_lib_coverage))+geom_histogram(binwidth=1)+ggtitle('coverage of ss libs in long reads')
lib.cov.hist.removetail.05 <- ggplot(lib.cov.all.removetail.05, aes(ss_lib_coverage))+geom_histogram(binwidth=1)+ggtitle('coverage of ss libs in long reads removing 5% quantile tail')

plt <- grid.arrange(grobs=list(cov.hist.removetail.05, lib.cov.hist, lib.cov.hist.removetail.05), ncol=3)

ggsave(filename=snakemake@output[[1]], plt)
