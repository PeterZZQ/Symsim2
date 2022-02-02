rm(list = ls())
gc()
library("devtools")
library("ape")
# setwd("~/Symsim2")
load_all("SymSim2")
########################################################################
#
# generate region to gene matrix
#
########################################################################

# A gene is regulated by k consecutive regions
# regions 1 to nregion are considered sequentially located on the genome
ngenes <- 1000
nregions <- 3000 
seed <- 1
# the probability that a gene is regulated by respectively 0, 1, 2 regions
p0 <- 0.01
prob_REperGene <- c(p0, (1-p0)*(1/10), (1-p0)*(1/10),(1-p0)*(1/5), (1-p0)*(1/5),(1-p0)*(1/5),(1-p0)*(1/10),(1-p0)*(1/10))
cumsum_prob <- cumsum(prob_REperGene)

region2gene <- matrix(0, nregions, ngenes)
set.seed(seed)
rand_vec <- runif(ngenes)

for (igene in 1:ngenes){
  if (rand_vec[igene] >= cumsum_prob[1] & rand_vec[igene] < cumsum_prob[2]) {
    region2gene[round(runif(1,min = 1, max = nregions)),igene] <- 1 
  } else if (rand_vec[igene] >= cumsum_prob[2]){
    startpos <- round(runif(1,min = 1, max = nregions-1))
    region2gene[startpos: (startpos+1),igene] <- c(1,1)
  }
}


ncells_total <- 10000
min_popsize <- 100

########################################################################
#
# Define trajectory structure
#
########################################################################
# linear
# phyla <- read.tree(text="(t1:1);")
# bifur
# phyla <- read.tree(text="(t1:1,t2:0.1);")
# trifur
phyla <- read.tree(text="(t1:1,t2:1.5,t3:0.5);")
plot(phyla)

########################################################################
#
# Simulate true scATAC-Seq and scRNA-Seq
#
########################################################################
# simulate the true count, the parameter setting the same as symsim
true_counts_res <- SimulateTrueCounts(ncells_total=ncells_total,min_popsize=min_popsize,i_minpop=2,ngenes=dim(region2gene)[2], 
                                      nregions=dim(region2gene)[1],region2gene=region2gene,atac_effect=0.8,
                                      evf_center=1,evf_type="continuous",nevf=12,
                                      n_de_evf=8,n_de_evf_atac = 3, impulse=F,vary='s',Sigma=0.30,
                                      phyla=phyla,geffect_mean=0,gene_effects_sd=1,gene_effect_prob=0.3,
                                      bimod=0,param_realdata="zeisel.imputed",scale_s=0.8,
                                      prop_hge=0.015, mean_hge=5, randseed=seed, gene_module_prop=0)
atacseq_data <- true_counts_res[[2]]
rnaseq_data <- true_counts_res[[1]]

# plot the tsne of true scRNA-Seq count 
tsne_rnaseq_true <- PlotTsne(meta=true_counts_res[[4]], data=log2(rnaseq_data+1),
                             evf_type="continous", n_pc=20, label='pop', saving = F, plotname="true rna-seq")

tsne_atacseq_true <- PlotTsne(meta=true_counts_res[[4]], data=log2(atacseq_data+1),
                             evf_type="continous", n_pc=20, label='pop', saving = F, plotname="true atac-seq")

########################################################################
#
# Simulate technical noise
#
########################################################################
# generate observed scRNA-Seq count
data(gene_len_pool)
gene_len <- sample(gene_len_pool, ngenes, replace = FALSE)
observed_rnaseq <- True2ObservedCounts(true_counts=true_counts_res[[1]], meta_cell=true_counts_res[[4]], 
                                       protocol="UMI", alpha_mean=0.2, alpha_sd=0.05, 
                                       gene_len=gene_len, depth_mean=5e5, depth_sd=3e3)
# cell_meta includes batch column
observed_rnaseq_loBE <- DivideBatches(observed_rnaseq, nbatch = 2, batch_effect_size = 1)
# Plot the tsne of observed scRNA-Seq count
tsne_rnaseq_noisy <- PlotTsne(meta=true_counts_res[[4]], data=log2(observed_rnaseq[[1]]+1), 
                              evf_type="continuous", n_pc=20, label='pop', saving = F, plotname="noisy rna-seq")

# generated observed scATAC-Seq count
atacseq_data <- round(atacseq_data)
atacseq_noisy <- atacseq_data
for (icell in 1:ncells_total){
  for (iregion in 1:nregions){
    if (atacseq_data[iregion, icell]>0){
      atacseq_noisy[iregion, icell] <- rbinom(n=1, size = atacseq_data[iregion, icell], prob = 0.3)}
    if (atacseq_noisy[iregion, icell] > 0){
      atacseq_noisy[iregion, icell] <- atacseq_noisy[iregion, icell]+rnorm(1, mean = 0, sd=atacseq_noisy[iregion, icell]/10)
    }
  }
}
atacseq_noisy[atacseq_noisy<0.1] <- 0
prop_1 <- sum(atacseq_noisy>0.1)/(dim(atacseq_noisy)[1]*dim(atacseq_noisy)[2])
target_prop_1 <- 0.1
if (prop_1 > target_prop_1) { # need to set larger threshold to have more non-zero values become 0s
  n2set0 <- ceiling((prop_1 - target_prop_1)*dim(atacseq_data)[1]*dim(atacseq_data)[2])
  threshold <- sort(atacseq_noisy[atacseq_noisy>0.1])[n2set0]
  atacseq_noisy[atacseq_noisy<threshold] <- 0
} else {
  print(sprintf("The proportion of 1s is %4.3f", prop_1))
}
# Plot the tsne of observed scATAC-Seq count
tsne_atacseq_noisy <- PlotTsne(meta=true_counts_res[[4]], data=atacseq_noisy, 
                               evf_type="continuous", n_pc=20, label='pop', saving = F, plotname="noisy atac-seq")

cellset_batch1 <- which(observed_rnaseq_loBE[[2]]$batch==1)
cellset_batch2 <- which(observed_rnaseq_loBE[[2]]$batch==2)

########################################################################
#
# Save the result
#
########################################################################
datapath <- sprintf("./seed_%d", seed)
system(sprintf("mkdir %s", datapath))
write.table(region2gene, sprintf("%s/region2gene.txt", datapath),
            quote=F, row.names = F, col.names = F, sep = "\t")
write.table(observed_rnaseq[[1]][, cellset_batch1], sprintf("%s/GxC1.txt", datapath), quote=F, row.names = F, col.names = F, sep = "\t")
write.table(observed_rnaseq[[1]][, cellset_batch2], sprintf("%s/GxC2.txt", datapath), quote=F, row.names = F, col.names = F, sep = "\t")
write.table(atacseq_noisy[, cellset_batch1], sprintf("%s/RxC1.txt", datapath), quote=F, row.names = F, col.names = F, sep = "\t")
write.table(atacseq_noisy[, cellset_batch2], sprintf("%s/RxC2.txt", datapath), quote=F, row.names = F, col.names = F, sep = "\t")
write.table(true_counts_res[[4]][cellset_batch1,1:2], sprintf("%s/cell_label1.txt", datapath), quote=F, row.names = F, col.names = T, sep = "\t")
write.table(true_counts_res[[4]][cellset_batch2,1:2], sprintf("%s/cell_label2.txt", datapath), quote=F, row.names = F, col.names = T, sep = "\t")
write.table(true_counts_res[[4]][cellset_batch1,3], sprintf("%s/pseudotime1.txt", datapath), 
            quote=F, row.names = F, col.names = F, sep = "\t")
write.table(true_counts_res[[4]][cellset_batch2,3], sprintf("%s/pseudotime2.txt", datapath), 
            quote=F, row.names = F, col.names = F, sep = "\t")

# save the plots
pdf(file=sprintf("%s/tsne.pdf", datapath))
print(tsne_rnaseq_true[[2]])
print(tsne_atacseq_true[[2]])
print(tsne_rnaseq_noisy[[2]])
print(tsne_atacseq_noisy[[2]])
dev.off()


