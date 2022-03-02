##################
#
# Figure 4
#
###################

message("Generating Fig. 4")

# in vitro
outdir = paste0("figures/figure4/")
dir.create(outdir)

id = "in_vitro"
sin_cl = scdb_mc(id); sin_mat = scdb_mat(id)

cells = names(sin_cl@mc)
umis = read_large_umis(id, cells = cells)
fp = sin_cl@mc_fp
lfp = log2(sin_cl@mc_fp)

sin_stats = sin_mat@cell_metadata[names(sin_cl@mc),]
sin_stats[ sin_stats$sorting.scheme == "Trbc+", "sorting.scheme"] = "TCRb+"
sin_stats[ sin_stats$treatment == "OVA + LPS T only", "treatment"] = "T mono-culture"
sin_stats[ sin_stats$date == 20200423, "replicate"] = 1

sin_comb = with(sin_stats, paste0(sorting.scheme, "@", PIC, "@", tissue, "@", treatment, "@", as.numeric(factor(date)), "-", replicate))

names(sin_comb) = rownames(sin_stats)
sin_umis = as.matrix(sin_mat@mat[, names(sin_cl@mc)])
sin_n = sweep(sin_umis,2,colSums(sin_umis),"/") * 1000

color_scheme = sin_cl@color_key
color2name = as.vector(color_scheme$group); names(color2name) = color_scheme$color
name2color = as.vector(color_scheme$color); names(name2color) = color_scheme$group
sin_names = color2name[ sin_cl@colors[ sin_cl@mc]]; names(sin_names) = names(sin_cl@mc)

lin_ord = c("T", "Activated_T", "Tht-like", "DC", "B") #in vitro

umicount = colSums(umis)

bad_mc = colnames(lfp)[ lfp["Trbc2",] > -2 & lfp["Cst3",] > 1]
clust_ord = setdiff(order(factor(color2name[ sin_cl@colors], levels = lin_ord)), bad_mc)

t_pops = c("T", "Activated_T", "Tht-like")
t_clusts = intersect(clust_ord, which(color2name[ sin_cl@colors] %in% t_pops))
t_cells = setdiff(names(sin_cl@mc)[ sin_cl@mc %in% t_clusts], names(which(umicount > 1e4)))
bad_genes = grep("^Rpl|^Rps|^Gm[0-9]|^AC[0-9]|^SNOR|^His[0-9]|^Hsp[0-9]", rownames(umis), v=T)

############

all_genes = read.table(paste0(outdir, "/Fig4b.txt"), stringsAsFactor=F)[,1]
ds = .downsamp(sin_umis, 500)
C = cor(t(log(1 + ds[all_genes, intersect(t_cells, colnames(ds))])))
diag(C) = NA
hc = hclust(as.dist(1-C), "ward.D2"); small_ct = 3 - cutree(hc,2)

grad = colorRampPalette(c("turquoise3", "white", "violetred"))(1000)
pdf(paste0(outdir, "/Fig4b.pdf"), height=14, width=14, useDingbats=F)
par(mar = c(5,5,0,0))
image.2(C, balance=T, hct=small_ct, vct = small_ct, annotate="both", col = grad); box()
dev.off()

#########

nms = names(sort(small_ct))
pdf(paste0(outdir, "/Fig4c.pdf"), height=14, width=14, useDingbats=F)
par(mar = c(5,5,0,0))
cell_ord = plot_sc_heatmap(id, id, nms, annotate=F, cells = t_cells, good_clusts = intersect(clust_ord, t_clusts)); box()
dev.off()

pdf(paste0(outdir, "/Fig4c_bottom.pdf"), height=3, width=14, useDingbats=F)
par(mar = c(0,5,0,0))
image(matrix(seq_along(cell_ord)), col = sin_cl@colors[ sin_cl@mc[ cell_ord]], axes=F); box()
dev.off()

########

sub_cells = intersect(t_cells, rownames(sin_stats)[sin_stats$sorting.scheme == "TCRb+"])
tfh_dist = table(sin_comb[sub_cells], sin_names[sub_cells] == "Tht-like")
val = tfh_dist[,2] / rowSums(tfh_dist)
condition = factor(vecsplit(names(val), "@", 4), levels = c("T mono-culture", "B16+T co-culture", "B16+T+DC co-culture"))

pdf(paste0(outdir, "/Fig4d.pdf"), height=10, width=7, useDingbats=F)
par(mar = c(10,3,1,1))
x = as.numeric(condition) + runif(length(condition), -0.2, 0.2)
plot(x, val, xlim = quantile(x, c(0,1)) + c(-0.5,0.5), axes=F, xlab = "", ylab="", type="n")
grid(); axis(2, las=2); axis(1, at = seq_along(levels(condition)), labels = levels(condition), las=2)
points(x, val, pch=21, bg = name2color["Tht-like"], cex=3, lwd=1)
dev.off()

sub_cells = intersect(t_cells, rownames(sin_stats)[sin_stats$sorting.scheme == "TCRb+CXCR5+PD1+"])
tfh_dist = table(sin_comb[sub_cells], sin_names[sub_cells] == "Tht-like")
val = tfh_dist[,2] / rowSums(tfh_dist)
condition = factor(vecsplit(names(val), "@", 4), levels = c("T mono-culture", "B16+T co-culture", "B16+T+DC co-culture"))

pdf(paste0(outdir, "/Fig4e.pdf"), height=10, width=7, useDingbats=F)
par(mar = c(10,3,1,1))
x = as.numeric(condition) + runif(length(condition), -0.2, 0.2)
plot(x, val, xlim = quantile(x, c(0,1)) + c(-0.5,0.5), axes=F, xlab = "", ylab="", type="n")
grid(); axis(2, las=2); axis(1, at = seq_along(levels(condition)), labels = levels(condition), las=2)
points(x, val, pch=21, bg = name2color["Tht-like"], cex=3, lwd=1)
dev.off()

#############
# compare to other datasets
# in vitro activation

all_mat = scdb_mat("pic-seq_all")
all_stats = all_mat@cell_metadata

act_cells = rownames(all_stats)[ all_stats$tissue == "in vitro" & all_stats$treatment == "OVA + LPS"]
naive_cells = rownames(all_stats)[ all_stats$tissue == "in vitro" & all_stats$treatment == "OVA + LPS T only"]
old_cells = union(act_cells, naive_cells)

act_umis = read_large_umis("pic-seq_all", cells = old_cells)
act_umis = act_umis[ rownames(umis),]
act_umis = act_umis[, colSums(act_umis) >= 500]
act_cells = colnames(act_umis)
act_stats = all_stats[act_cells,]
act_t_cells = rownames(act_stats)[ act_stats$sorting.scheme == "Trbc+"]

comb_umis = cbind(umis[,t_cells], act_umis[ rownames(umis), act_t_cells])
comb_names = c(sin_names[t_cells], ifelse(act_stats[ act_t_cells, "treatment"] == "OVA + LPS", "LPS_activated", "LPS_naive")); 
names(comb_names) = colnames(comb_umis)
comb_n = sweep(comb_umis, 2, colSums(comb_umis), "/") * 1000

good_genes = setdiff(rownames(comb_n), bad_genes)
m = t(apply(comb_n[good_genes,], 1, tapply, comb_names, mean)) * min(table(comb_names))
m = m[,c("T", "Activated_T", "Tht-like", "LPS_naive", "LPS_activated")]

#########

z_tvl = log2((10 + m[,"Tht-like"]) / (10 + m[,"LPS_activated"]))
z_tvn = log2((10 + m[,"Tht-like"]) / (10 + m[,"T"]))
z_lvn = log2((10 + m[,"LPS_activated"]) / (10 + m[,"LPS_naive"]))
z_a = apply(cbind(z_tvn, z_lvn),1,max)

genes = setdiff(union(scr_chi_square_diff_genes(comb_umis, g1 = names(which(comb_names == "LPS_activated")),
	g2 = names(which(comb_names == "Tht-like")), fdr = T, pval = 1e-3),
	scr_chi_square_diff_genes(comb_umis, g1 = names(which(comb_names == "LPS_activated" | comb_names == "Tht-like")),
	g2 = names(which(comb_names == "T" | comb_names == "LPS_naive")), fdr = T, pval = 1e-3)), bad_genes)

disp_genes = genes #intersect(genes, names(which(abs(z_tvl) > 1 | abs(z_a) > 0)))
disp_genes = c(names(head(sort(z_tvl[intersect(genes, names(which(z_a > 1)))],T),20)), 
	names(head(sort(z_tvl[intersect(genes, names(which(z_a > 2)))],F),20)), 
	names(head(sort(z_a[disp_genes],T),30)), 
	names(head(sort(z_a[disp_genes],F),15)))

must_haves = c("Btla", "Ctla4", "Tigit", "Pdcd1", "Maf")
disp_genes = union(disp_genes, must_haves)
#refined_genes = names(which(z[disp_genes] > 0))
df = data.frame(x=z_tvl, y=z_a,
        color = names(z_a) %in% disp_genes, text = ifelse(names(z_a) %in% disp_genes, names(z_a), ""))
df = df[ setdiff(rownames(df), bad_genes),]
df$rx = round(df$x, 3); df$ry = round(df$y, 3)
df = df[!duplicated(df[,c("color", "rx", "ry")]),]
df = df[order(df$color),]
xlim = max(abs(df$x)) * c(-1,1); ylim = max(abs(df$y)) * c(-1,1)
ggsave(paste0(outdir, "/Fig4f.pdf"),
        ggplot(df, aes(x,y, color = factor(color), label=text)) + geom_point() +
	geom_hline(yintercept=0) + geom_vline(xintercept=0) + xlim(xlim) + ylim(ylim) + 
        geom_text_repel(color="black", segment.color = "gray20", segment.alpha = 0.5, nudge_x=0.05, nudge_y=0.1) + scale_color_manual(values=c("gray70", "tomato2")) +
	theme(legend.position = "none"), useDingbats=F)

activation_genes = intersect(genes, names(which(z_tvl > 1 & z_tvn > 1)))
write.table(activation_genes, row.names=F, quote=F, col.names=F, file = paste0(outdir, "/activation_genes.txt"))
write.table(z_tvn, sep = "\t", quote=F, col.names=NA, file = paste0(outdir, "/tht_vs_naive.txt"))

