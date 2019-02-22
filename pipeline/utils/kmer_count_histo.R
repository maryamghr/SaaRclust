library(data.table)
library(ggplot2)
library(gridExtra)

args = commandArgs(trailingOnly=TRUE)
print(args)

kmer.count <- fread(args[1])
trim = as.numeric(args[2])
binwidth = as.numeric(args[3])

q = quantile(kmer.count[, km], probs=c(trim, 1-trim))

plt1 = ggplot(kmer.count[km > q[1] & km < q[2]]) + geom_histogram(aes(x=km, y=..density.., fill=hetstatus), alpha=0.5, position = 'identity', binwidth=binwidth) + ggtitle("kmer count density")
plt2 = ggplot(kmer.count[km > q[1] & km < q[2]]) + geom_histogram(aes(x=km, fill=hetstatus), alpha=0.5, position = 'identity', binwidth=binwidth) + ggtitle("kmer count frequency")

ggsave(args[4], grid.arrange(plt1, plt2, nrow=2))
