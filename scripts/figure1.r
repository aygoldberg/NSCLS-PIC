message("Generating Fig. 1 and S1")
set.seed(1111)

id = "human_singlets"
sin_2d = scdb_mc2d(id); sin_cl = scdb_mc(id); sin_mat = scdb_mat(id)

cells = names(sin_cl@mc)
cell_stats = sin_mat@cell_metadata[cells,]
fp = sin_cl@mc_fp
lfp = log2(sin_cl@mc_fp)

outdir = paste0("figures/figure1")
supdir = paste0("figures/figureS1")
dir.create(outdir)
dir.create(supdir)

sin_stats = sin_mat@cell_metadata[names(sin_cl@mc),]
sin_umis = read_large_umis(id, cells = cells)
sin_n = sweep(sin_umis,2,colSums(sin_umis),"/") * 1000

color_scheme = sin_cl@color_key
color2name = as.vector(color_scheme$group); names(color2name) = color_scheme$color
name2color = as.vector(color_scheme$color); names(name2color) = color_scheme$group
sin_names = color2name[ sin_cl@colors[ sin_cl@mc]]; names(sin_names) = names(sin_cl@mc)

sin_comb = with(sin_stats, paste0(sorting.scheme, ":", PIC, "@", organ, "@", site, "@", patient))
names(sin_comb) = rownames(sin_stats)

anno_colors = as.matrix(read.delim("annotations/human_lung_annotations.txt", stringsAsFactor=F, row.names=1))
anno_names = anno_colors[,2]
lin_ord = names(anno_names)
umicount = colSums(sin_umis)

bad_clusts = colnames(lfp)[ lfp["TRBC2",] > -1 & lfp["CST3",] > 2]

t_pops = names(which(anno_colors[,3] == "T"))
dc_pops = names(which(anno_colors[,3] == "Myeloid"))

t_clusts = setdiff(which(color2name[ sin_cl@colors] %in% t_pops), bad_clusts)
t_cells = setdiff(names(sin_cl@mc)[ sin_cl@mc %in% t_clusts], names(which(umicount > 1e4)))

dc_clusts =  setdiff(which(color2name[ sin_cl@colors] %in% dc_pops), bad_clusts)
dc_cells = setdiff(names(sin_cl@mc)[ sin_cl@mc %in% dc_clusts], names(which(umicount > 1e4)))

cells = union(t_cells, dc_cells)

##############

lateral_genes = as.matrix(read.delim("annotations/lateral_genes.txt", stringsAsFactor=F, row.names=1))[,1]
human_lateral = names(which(lateral_genes == "human_lateral"))
bad_genes = grep("^RPL|^RPS|^MT[A-Z]|^HLA-|^AC[0-9]|^SNOR|^HIST[0-9]|^HSP[0-9]|^MIR[0-9]", rownames(sin_umis), v=T)
bad_markers = c("Metazoa_SRP", "MTND1P23", "MALAT1",bad_genes, human_lateral)

lr_features =  choose_lr_features(id, t_cells, dc_cells, bad_genes, must_haves = names(scdb_gset(id)@gene_set))
mle_features = choose_mle_features(id, id, t_cells, dc_cells, union(bad_markers, bad_genes), existing_list= names(scdb_gset(id)@gene_set))

mle_features = union(mle_features, "CXCL13")

############

mc = names(table(sin_cl@mc[cells]))
mc = mc[ order(factor(color2name[ sin_cl@colors[ as.numeric(mc)]], levels = lin_ord))]
patient = as.vector(sin_stats$patient); names(patient) = rownames(sin_stats)
sample_dist = table(sin_cl@mc[cells], patient[cells])
sample_dist = t(apply(sample_dist[mc,],1,sort,T))
dist_n = sample_dist / rowSums(sample_dist)

cols = colorRampPalette(c("chocolate4", "white"))(length(table(patient)))
pdf(paste0(supdir, "/FigS1e.pdf"), useDingbats=F) #, height=1000, width=1500)
par(mar = c(0.5,3,0.5,0.5), fig=c(0,1,0.1,1))
barplot(t(dist_n), xaxs = "i", names.arg = rep("", length(mc)), col = cols, las=2)
par(fig	= c(0,1,0,0.1), new=T)
image(matrix(seq_along(mc)), axes=F, col = sin_cl@colors[ as.numeric(mc)]); box(lwd=2)
dev.off()

#############

pdf(paste0(supdir, "/blue_red_colorbar.pdf"), useDingbats=F)
par(mar=rep(0,4))
image(matrix(seq_len(1000)), col=colorRampPalette(c("blue", "white", "red"))(1000), axes=F); box()
dev.off()

bad_genes = grep("^RPL|^RPS|^MT[A-Z]|^HLA-|^AC[0-9]|^SNOR|^HIST[0-9]|^HSP[0-9]|^MIR[0-9]", rownames(sin_umis), v=T)

t_nms = rev(read.table(paste0(outdir, "/Fig1d.txt"), stringsAsFactor=F)[[1]])
#t_nms = setdiff(choose_genes_from_clust(paste0(id, "_a"), id, nms_per_clust=5, nms_thresh=1.4, ord = "max.col", must_haves = c(t_nms, "CXCR4", "CXCR6")), bad_genes)
t_ord = intersect(mc, t_clusts)
t_cl = scdb_mc(paste0(id, "_a"))
t_fp = log2(t_cl@mc_fp); colnames(t_fp) = names(table(sin_cl@mc[t_cells]))

#vct = factor(color2name[ sin_cl@colors[ as.numeric(t_ord)]], levels = intersect(lin_ord, names(table(sin_names[t_cells])))); names(vct) = t_ord
IM = t_fp[t_nms, as.character(t_ord)]
IM = IM[ order(max.col(IM)),]
vct = factor(color2name[ sin_cl@colors[ as.numeric(colnames(IM))]], levels = lin_ord); names(vct) = colnames(IM)
pdf(paste0(outdir, "/Fig1d.pdf"), useDingbats=F, height=10, width=7)
par(mar = c(0.5,5,0.5,0.5), fig = c(0,1,0.1,1))
image.2(IM, balance=T, vct=vct, annotate="rows"); box()
par(fig = c(0,1,0,0.1), new=T)
image(matrix(seq_len(ncol(IM))), axes=F, col = sin_cl@colors[ as.numeric(colnames(IM))]); box()
dev.off()

dc_nms = rev(read.table(paste0(outdir, "/Fig1e.txt"), stringsAsFactor=F)[[1]])
#dc_nms = setdiff(choose_genes_from_clust(paste0(id, "_b"), id, nms_per_clust=5, nms_thresh=1.4, ord = "max.col", must_haves = dc_nms), bad_genes)
dc_ord = intersect(mc, dc_clusts)
dc_cl = scdb_mc(paste0(id, "_b"))
dc_fp = log2(dc_cl@mc_fp); colnames(dc_fp) = names(table(sin_cl@mc[dc_cells]))

IM = dc_fp[intersect(dc_nms, rownames(dc_fp)), as.character(dc_ord)]
IM = IM[ order(max.col(IM)),]
vct = factor(color2name[ sin_cl@colors[ as.numeric(colnames(IM))]], levels = lin_ord); names(vct) = colnames(IM)
pdf(paste0(outdir, "/Fig1e.pdf"), useDingbats=F, height=10, width=7)
par(mar = c(0.5,5,0.5,0.5), fig = c(0,1,0.1,1))
image.2(IM, balance=T, vct=vct, annotate="rows"); box()
par(fig = c(0,1,0,0.1), new=T)
image(matrix(seq_len(ncol(IM))), axes=F, col = sin_cl@colors[ as.numeric(colnames(IM))]); box()
dev.off()

########

bad_clusts = which( color2name[ sin_cl@colors] %in% c("Mast_cells", "B_Plasma"))
clust_ord = setdiff(order(factor(color2name[ sin_cl@colors], levels = lin_ord)), bad_clusts)
good_mc = intersect(clust_ord, c(t_clusts, dc_clusts))
confu = mcell_mc_confusion_mat(id, id, 1000, ignore_mismatch=F)
r_confu = rowSums(confu)
c_confu = colSums(confu)
norm = r_confu %*% t(c_confu)
confu_n = confu/norm

confu_n = confu_n[clust_ord, clust_ord]
confu_nodiag = confu_n
diag(confu_nodiag) = 0
confu_n = pmin(confu_n, max(confu_nodiag))
confu_n = pmin(confu_n, quantile(confu_n, 1-3/nrow(confu_n)))
shades = colorRampPalette(c("white",
        "pink", "orange", "brown3", "gray20"))(1000)

cls = factor(color2name[sin_cl@colors[ as.numeric(clust_ord)]], levels = lin_ord); names(cls) = clust_ord
pdf(paste0(supdir, "/FigS1f.pdf"), useDingbats=F)
par(mar = rep(0.5,4), fig = c(0.05,1,0.05,1))
image.2(confu_n, annotate="none", col=shades, hct=cls, vct=cls); box()
par(fig = c(0,0.05,0.05,1), new=T)
image(t(matrix(seq_along(clust_ord))), axes=F, col = sin_cl@colors[ as.numeric(clust_ord)]); box()
par(fig = c(0.05,1,0,0.05), new=T)
image(matrix(seq_along(clust_ord)), axes=F, col = sin_cl@colors[ as.numeric(clust_ord)]); box()
dev.off()

pdf(paste0(supdir, "/confusion_colorbar.pdf"), useDingbats=F)
par(mar=rep(0,4))
image(matrix(seq_along(shades)), col=shades, axes=F); box()
dev.off()

###########

pdf(paste0(outdir, "/Fig1c.pdf"), useDingbats=F)
plot(sin_2d@sc_x[cells], sin_2d@sc_y[cells], pch = 21, bg = sin_cl@colors[ sin_cl@mc[cells]], cex = 1,
	axes=F, xlab="", ylab="")
dev.off()

pdf(paste0(outdir, "/legend.pdf"), useDingbats=F)
plot(1,1,type="n", axes=F, xlab="", ylab="")
legend("topleft", anno_names[lin_ord], pch = 21, cex=1, pt.bg=name2color[ lin_ord])
dev.off()

site = as.vector(sin_stats$site); names(site) = rownames(sin_stats)
gating = as.vector(sin_stats$sorting.scheme); names(gating) = rownames(sin_stats)
patient = as.vector(sin_stats$patient); names(patient) = rownames(sin_stats)

comb2 = interaction(gating, site); names(comb2) = names(gating)
pdf(paste0(outdir, "/Fig1f.pdf"), useDingbats=F, height=10, width=10)
par(mar=rep(0.5,4), mfrow = c(2,2))
for (c in names(table(comb2))) {
	sub_cells = intersect(cells, names(comb2)[ comb2 == c])
	plot(sin_2d@sc_x[cells], sin_2d@sc_y[cells], pch = 20, col="gray80", axes=F, xlab="", ylab="", cex=0.7)
	points(sin_2d@sc_x[sub_cells], sin_2d@sc_y[sub_cells], pch = 20, cex=1.5,
	        col = sin_cl@colors[ sin_cl@mc[sub_cells]])
}
dev.off()

##############

comb = with(sin_stats, paste0(site, "@", patient)); names(comb) = rownames(sin_stats)
t_dist = table(comb[t_cells], sin_names[ t_cells])
t_n = t_dist / rowSums(t_dist)

dc_dist = table(comb[dc_cells], sin_names[ dc_cells])
dc_dist = dc_dist[,setdiff(colnames(dc_dist), c("NK_CX3CR1", "NK"))]
dc_dist = dc_dist[ rowSums(dc_dist) > 0,]
dc_n = dc_dist / rowSums(dc_dist)

shared_samps = intersect(rownames(t_dist), rownames(dc_dist))
dist_n = cbind(t_n[shared_samps,], dc_n[shared_samps,])

C = cor(dist_n, m="spearman"); diag(C) = NA
hc = hclust(as.dist(1-C), "ward.D2");   ord = rownames(C)[hc$order]

reg = 0.2

dist_melt = melt(dist_n)
dist_melt$patient = vecsplit(as.character(dist_melt$Var1), "@", 2)
dist_melt$site = vecsplit(as.character(dist_melt$Var1), "@", 1)
dist_melt$comb = interaction(dist_melt$patient, factor(dist_melt$Var2, levels = ord))
merged = merge(dist_melt[dist_melt$site == "normal",], dist_melt[dist_melt$site == "tumor",], by.x = "comb", by.y = "comb")
merged$fc = log2((reg + merged$value.y) / (reg + merged$value.x))
merged$y = as.numeric(factor(merged$Var2.x, levels = ord)) + runif(nrow(merged), -0.1,0.1)
zlim = max(abs(merged$fc))

grad = colorRampPalette(c("turquoise3", "white", "violetred"))(1000)
rownames(C) = anno_names[ rownames(C)]; colnames(C) = rownames(C)
ord = anno_names[ord]
pdf(paste0(outdir, "/Fig1g-h.pdf"), useDingbats=F, height=8, width=12)
par(mar = c(10,10,1,1), fig = c(0,0.714,0,1))
image.2(C[ord,ord], balance=T, annotate="both", col=grad); box()
par(mar=c(10,1,1,10), fig = c(0.714,1,0,1), new=T)
with(merged, plot(fc, y, xlim = c(-zlim,zlim), ylim = round(quantile(y, c(0,1))) + c(-0.5,0.5), yaxs="i", pch = 21, bg = name2color[ as.character(Var2.x)],
	axes=F, xlab="", ylab="", cex=1))
abline(v=0, lwd=2)
axis(1); axis(4, at = seq_along(ord), labels = ord, las=2)
dev.off()

pdf(paste0(outdir, "/cor_cb.pdf"), useDingbats=F)
par(mar = rep(0.5,4))
image(matrix(seq_along(grad)), axes=F, col = grad); box()
dev.off()

