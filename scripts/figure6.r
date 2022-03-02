message("Generating Fig. 6 and S7")

outdir = paste0("figures/figure6/")
dir.create(outdir)

supdir = paste0("figures//figureS7/")
dir.create(supdir)

##########

breast_id = "breast"
breast_cl = scdb_mc(breast_id)
breast_2d = scdb_mc2d(breast_id)

breast_mc = which(breast_cl@colors != "white")
breast_cells = names(breast_cl@mc)[ breast_cl@mc %in% breast_mc]
breast_lfp = log2(breast_cl@mc_fp)
breast_umis = read_large_umis(breast_id, cells = breast_cells)
breast_n = sweep(breast_umis,2,colSums(breast_umis),"/") * 1000

color_scheme = breast_cl@color_key
color2name = as.vector(color_scheme$group); names(color2name) = color_scheme$color
name2color = as.vector(color_scheme$color); names(name2color) = color_scheme$group
breast_names = color2name[ breast_cl@colors[ breast_cl@mc]]; names(breast_names) = names(breast_cl@mc)

anno_colors = as.matrix(read.delim("annotations/human_lung_annotations.txt", stringsAsFactor=F, row.names=1))
anno_names = anno_colors[,2]
lin_ord = names(anno_names)
t_ord = setdiff(c(lin_ord[1:8], "Tht_1", "Tht_2", lin_ord[10], "IFN"), "T_CD4_memory")

breast_stats = scdb_mat(breast_id)@cell_metadata
tcr_cells = intersect(breast_cells, rownames(breast_tcr))
cell2patient = breast_stats[tcr_cells, "patient_id"]; names(cell2patient) = tcr_cells
cell2clone_id = breast_tcr[tcr_cells, "cdr3_nt"]; names(cell2clone_id) = tcr_cells
cell2tp = breast_stats[tcr_cells, "timepoint"]; names(cell2tp) = tcr_cells

###########

comb = with(breast_stats, paste0(patient_id, "@", expansion, "@", BC_type, ":", timepoint)); names(comb) = rownames(breast_stats)

pops = t_ord
good_cells = names(breast_names)[ breast_names %in% pops]
sample_dist = table(comb[ good_cells], breast_names[ good_cells])
sample_dist = sample_dist[ rowSums(sample_dist) > 50,]
samp2pat = vecsplit(rownames(sample_dist), ":", 1)
good_samps = names(which(table(samp2pat) == 2))

dist_n = sample_dist  / rowSums(sample_dist)
x = dist_n[ paste0(good_samps, ":Pre"),]
y = dist_n[ paste0(good_samps, ":On"),]

good_pops = c("Tht_1", "Tht_2")
dist_melt = melt(x)
colnames(dist_melt) = c("patient", "pop", "pre")
dist_melt$post = melt(y)$value
dist_melt$pre_x=as.numeric(factor(dist_melt$pop, levels = good_pops)) * 3 #runif(nrow(dist_melt), -0.1, 0.1)
#dist_melt$pre_x = dist_melt$prx + 1
dist_melt$post_x = dist_melt$pre_x + 1
dist_melt$pre_x_jitter = dist_melt$pre_x + runif(nrow(dist_melt), -0.1, 0.1)
dist_melt$post_x_jitter = dist_melt$pre_x_jitter + 1
dist_melt = dist_melt[ dist_melt$pop %in% good_pops,]
xlim = with(dist_melt, quantile(c(pre_x, post_x), c(0,1), na.rm=T))
ylim = with(dist_melt, quantile(c(pre, post), c(0,1)))

box_cols = rep(NA, xlim[2])
reg=0.02
good_pops = c("Tht_1", "Tht_2")
dist_melt$x = as.numeric(factor(dist_melt$pop, levels = good_pops)) + runif(nrow(dist_melt), -0.2, 0.2)
dist_melt$fc = with(dist_melt, log2((reg + post) / (reg + pre)))
fc_median = with(dist_melt, tapply(fc, pop, median))

ylim = max(abs(dist_melt[ !is.na(dist_melt$x), "fc"])) * c(-1,1)
pdf(paste0(outdir, "/Fig6a.pdf"), useDingbats=F, height=10, width=7)
par(mar = c(5,3,1,1))
with(dist_melt, plot(x, fc, pch = 21, bg = name2color[ rep(as.vector(pop))], axes=F, xlab="", ylab="", cex=1.5, ylim = ylim))
abline(h=0)
axis(2, las=2); axis(1, las=1, at = seq_along(good_pops), labels = good_pops)
dev.off()

##############

id = "ot2_apd-1"
sin_cl = scdb_mc(id); sin_mat = scdb_mat(id)

cells = names(sin_cl@mc)
umis = read_large_umis(id, cells = cells)
fp = sin_cl@mc_fp
lfp = log2(sin_cl@mc_fp)

sin_stats = sin_mat@cell_metadata[names(sin_cl@mc),]
sin_stats[ sin_stats$sorting.scheme == "Cd11c+", "sorting.scheme"] = "CD11c+"
sin_comb = with(sin_stats, paste0(sorting.scheme, "@", PIC, "@", tissue, "@", treatment, "@", timepoint, ":",  as.numeric(factor(date)), "-", replicate))
names(sin_comb) = rownames(sin_stats)
sin_group = vecsplit(sin_comb, ":", 1)
umis_n = sweep(umis,2,colSums(umis),"/") * 1000

#################

old_id = "ot2"
old_cl = scdb_mc(old_id); old_mat = scdb_mat(old_id)

old_cells = names(old_cl@mc)
old_umis = read_large_umis(old_id, cells = old_cells)
old_stats = old_mat@cell_metadata[old_cells,]
old_stats[ old_stats$sorting.scheme == "Cd11c+", "sorting.scheme"] = "CD11c+"
old_comb = with(old_stats, paste0(sorting.scheme, "@", PIC, "@", tissue, "@", timepoint))
names(old_comb) = rownames(old_stats)

both_umis = cbind(old_umis, umis)
both_n = sweep(both_umis,2,colSums(both_umis),"/") * 1000
all_stats = scdb_mat("mouse_singlets")@cell_metadata[ colnames(both_umis),]

t_cells = rownames(all_stats)[ all_stats$sorting.scheme %in% c("TCRb+", "TCRb+CD45.1+") & all_stats$timepoint %in% c("d10", "d17") &
        all_stats$treatment %in% c("none", "aPD-1", "Ova-Vaccination") &
        all_stats$tissue != "tumor"]
all_comb = with(all_stats, paste0(ifelse(treatment == "Ova-Vaccination", "old", "new"), ":",
	sorting.scheme, "@", PIC, "@", tissue, "@", timepoint, "@", treatment, ":", as.numeric(factor(date)), "-", replicate))
names(all_comb) = rownames(all_stats)
good_combs = names(which(table(all_comb[t_cells]) > 10))
t_cells = intersect(t_cells, names(all_comb)[ all_comb %in% good_combs])

bad_samps = c("TCRb+CD45.1+@Singlets@dLN@Ova-Vaccination@d10:4-2", grep("d17@Ova-Vaccination", good_combs, v=T))
bad_genes = grep("^Rpl|^Rps|^Gm[0-9]|^AC[0-9]|^Snor|^His[0-9]|^Hsp[0-9]|^Trim30|^mt-|^Ig", rownames(both_umis), v=T)
m = t(apply(both_n[ setdiff(rownames(both_n), bad_genes), t_cells], 1, tapply, all_comb[t_cells], mean))
m = m[, setdiff(colnames(m), bad_samps)]

samp_names = vecsplit(colnames(m), ":", 2)
samp_date = factor(vecsplit(colnames(m), ":", 1), levels = c("old", "new"))
samp2gate = factor(vecsplit(samp_names, "@", 1), levels = c("TCRb+", "TCRb+CD45.1+"))
samp2tissue = factor(vecsplit(samp_names, "@", 3), levels = c("cLN", "dLN"))
samp2tp = factor(vecsplit(samp_names, "@", 4), levels = c("d10", "d17"))
samp2treat = factor(vecsplit(samp_names, "@", 5), levels = c("Ova-Vaccination", "none", "aPD-1"))

samp2comb = interaction(interaction(samp2tissue, samp2treat), samp2gate)
samp2comb = factor(samp2comb, levels = names(which(table(samp2comb) > 0)))
names(samp2comb) = colnames(m)
m_med = t(apply(m, 1, tapply, samp2comb, median))
reg = 0.02
IM = log2((reg + m_med) / (reg + apply(m, 1, median)))
diff_genes = names(which(apply(abs(IM),1,max) > 1))
IM2 = IM[ diff_genes,]

library(tglkmeans)
k = round(nrow(IM2) / 20)
data = as.data.frame(IM2)
data$id = rownames(data)
data = data[,c(ncol(data), 1:(ncol(data) - 1))]
km <- TGL_kmeans_tidy(data, k=k, metric='euclid', verbose=TRUE, seed = 18)
centers = as.matrix(km$centers[,-1]); rownames(centers) = seq_len(k)
centers = centers[order(max.col(centers)),]
km_clusts = as.numeric(factor(km$cluster$clust, levels = rownames(centers))); names(km_clusts) = rownames(IM2)
rownames(centers) = seq_len(nrow(centers))

meta_mat = rbind(samp2gate, samp2tissue, samp2treat) + c(0,2,4)
colnames(meta_mat) = samp2comb[ colnames(meta_mat)]
cols = c("gray80", "gray20", "#E5E4D1", "#E5DE61", "#efc000", "#efc000", "#CD534C")
rownames(meta_mat) = c("Gate", "Tissue", "Treatment")
unique_meta = t(unique(t(meta_mat)))

pdf(paste0(supdir, "/FigS7a.pdf"), useDingbats=F, height=8, width=5)
par(mar = c(1,3,1,1), fig = c(0,1,0.2,1))
image.2(centers, b=T, annotate="rows"); box()
par(fig = c(0,1,0,0.2), new=T)
image.2(unique_meta[, colnames(IM2)], col = cols, annotate="rows"); box()
dev.off()

ks = c(4,11)
dir.create(paste0(supdir, "/FigS7a_individual_clusters"))
for (k in rownames(centers)) {
	genes = names(which(km_clusts == k))
	message(k)
	pdf(paste0(supdir, "/FigS7a_individual_clusters/", k, ".pdf"), useDingbats=F, height=8, width=5)
	par(mar = c(1,10,1,1), fig = c(0,1,0.2,1))
	image.2(IM2[genes,], b=T, annotate="rows", hct=km_clusts[genes]); box()
	par(fig = c(0,1,0,0.2), new=T)
	image.2(unique_meta[, colnames(IM2)], col = cols, annotate="rows"); box()
	dev.off()
}

#############

genes = c("Btla", "Cd200", "Cxcr5", "Hif1a", "Lag3", "Nfatc1", "Pou2f2",
	"Slamf6", "Tnfsf4", "Zbtb7b", "Vmp1", "Trp53inp1")
good_combs = c("dLN.none.TCRb+", "dLN.aPD-1.TCRb+", "dLN.none.TCRb+CD45.1+", "dLN.aPD-1.TCRb+CD45.1+")
m_rep = m[, samp2comb %in% good_combs]
sub2comb = factor(samp2comb[ colnames(m_rep)], levels = names(which(table(samp2comb[colnames(m_rep)]) > 0)))
means = t(apply(m_rep, 1, tapply, sub2comb[ colnames(m_rep)], mean))
#means = means[, colSums(is.na(means)) == 0]

cols = c(NA, "#EFC000", "#CD534C")
x = as.numeric(sub2comb)
l = levels(sub2comb)

pdf(paste0(supdir, "/FigS7c.pdf"), useDingbats=F, height=9, width=12)
par(mfrow = c(3,4), mar = c(3,3,1,1))
for (gene in genes) {
	plot(x + runif(length(x), -0.1, 0.1), m_rep[gene,], pch = 21, bg = cols[ as.numeric(samp2treat)], cex=2, xlim = c(0.5,4.5),
		ylim = c(0,max(m_rep[gene,] * 1.1)), main = gene, axes=F, xlab="", ylab="")
	segments(seq_along(l) -0.2, means[gene,l], seq_along(l) + 0.2, col = "red", lwd=3)
	axis(2, las=2); axis(1, at=c(1.5, 3.5), labels = levels(samp2gate))
}
dev.off()

###########