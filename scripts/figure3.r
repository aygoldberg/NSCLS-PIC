message("Generating Fig. 3, S3 and S4")

# generate Fig3c-d

outdir = paste0("figures/figure3/")
dir.create(outdir)

supdir = paste0("figures//figureS3-4/")
dir.create(supdir)

############

ds_f = ds[setdiff(rownames(sin_cl@e_gc), bad_genes), good_pics]
us = umis[ rownames(ds_f), colnames(ds_f)]
good_genes = setdiff(rownames(us), bad_genes)

exp_us = generate_expected_pics_from_mle(id, mle_res[good_pics, c("a_mc", "b_mc")], mle_res[good_pics, "alpha"],
	colSums(us[,good_pics]), bad_genes = bad_genes)
exp_n = generate_expected_pics_from_mle(id, mle_res[good_pics, c("a_mc", "b_mc")], mle_res[good_pics, "alpha"],
	colSums(ds_f[,good_pics]), bad_genes = bad_genes)
t_n = generate_expected_pics_from_mle(id, mle_res[good_pics, c("a_mc", "b_mc")], rep(1, length(good_pics)),
        colSums(ds_f[,good_pics]) * alpha[good_pics], bad_genes = bad_genes)
dc_n = generate_expected_pics_from_mle(id, mle_res[good_pics, c("a_mc", "b_mc")], rep(0, length(good_pics)),
        colSums(ds_f[,good_pics]) * (1 - alpha[good_pics]), bad_genes = bad_genes)
genes = rownames(exp_us)

############

pop = "Tht"
sub_cells = intersect(rownames(cell_stats)[ cell_stats$site == "tumor"], names(which(parser_t == pop)))

diff_genes = read.table(paste0(outdir, "/Fig3c.txt"), stringsAsFactors=F, h=T)[[1]]

reg=20
y = rowSums(us[genes,sub_cells]); x = rowSums(exp_us[genes,sub_cells])
sum_t = rowSums(t_n[genes, sub_cells]); sum_dc = rowSums(dc_n[genes, sub_cells])
z = log2((reg + y) / (reg + x))
disp_genes = union(names(which(abs(z) > 1)), diff_genes)
col_z = log2((reg + sum_dc) / (reg + sum_t))
grad = colorRampPalette(c("limegreen", "gray40", "firebrick3"))(101)
val = col_z; zlim = max(abs(val[ disp_genes]))
val_n = round((val + zlim) / (2 * zlim) * 100) + 1

df = data.frame(gene = names(x), x=log2(reg + x), y=log2(reg + y), 
        fill = val_n, 
        text = ifelse(names(x) %in% disp_genes, good_genes, ""))
df = df[ setdiff(rownames(df), bad_genes),]
df$rx = round(df$x, 3); df$ry = round(df$y, 3)
df = df[!duplicated(df[,c("text", "rx", "ry")]),]
#df = df[order(df$color),]
lim = quantile(c(df$x, df$y), c(0,1))

pdf(paste0(outdir, "/Fig3c.pdf"), useDingbats=F)
with(df, plot(x, y, pch = 20, col = "gray", cex = 2,
        xlim = lim, ylim = lim, axes = F, xlab = "", ylab = ""))
axis(1, las=2); axis(2, las=2);
abline(coef = c(1,1), lty = 2); abline(coef = c(-1,1), lty = 2); abline(coef = c(0,1))
with(df[disp_genes,], points(x, y, cex = 2, pch = 21, bg = grad[fill]))
#with(df[disp_genes,], text(x, y, gene, adj = 1))
dev.off()

##############

pic_comb = parser_t[good_pics] 
pic_comb[ !(pic_comb %in% c("Tht", "T_CD4_reg"))] = "Other"
pic_comb = factor(pic_comb, levels = c("Other", intersect(lin_ord, names(table(pic_comb))))); names(pic_comb) = good_pics

comb = with(cell_stats, paste0(site, "@", patient)); names(comb) = rownames(cell_stats)
comb2 = interaction(comb[good_pics], pic_comb, sep = "@")
comb2 = factor(comb2, levels = names(which(table(comb2) > 0))); names(comb2) = good_pics
cells2 = intersect(good_pics, colnames(ds_f));
cells2 = cells2[ comb2[cells2] %in% names(which(table(comb2[cells2]) >= 5))]
comb2 = factor(comb2, levels = names(which(table(comb2[cells2]) > 0))); names(comb2) = good_pics

disp_genes = c("CXCL13", "CTLA4", "MAF", "CCL4", "IGF2R", "DAP", "RAN", "PDCD1")
db_int = comb2
real_m = t(apply(ds_f[disp_genes	,cells2], 1, tapply, db_int[cells2], mean));
exp_m =  t(apply(exp_n[disp_genes	,cells2], 1, tapply, db_int[cells2], mean));
t_m =  t(apply(t_n[disp_genes	,cells2], 1, tapply, db_int[cells2], mean));
dc_m =  t(apply(dc_n[disp_genes	,cells2], 1, tapply, db_int[cells2], mean));

############

reg = 0.001
grad = colorRampPalette(c("limegreen", "gray40", "firebrick3"))(101)
high_patients = c(1,6,7,9)
disp_genes = c("CXCL13", "CTLA4", "MAF", "CCL4", "IGF2R", "DAP", "RAN", "PDCD1")
pdf(paste0(outdir, "/Fig3d.pdf"), useDingbats=F, height=14, width=7)
par(mar = c(2,2,1,1), mfrow = c(4,2))
for (gene in disp_genes) {
	obs_vals = real_m[gene,]; 
	vals_melt = melt(obs_vals)
	vals_melt$exp_value = melt(exp_m[gene,])$value
	vals_melt$t_value = melt(t_m[gene,])$value
	vals_melt$dc_value = melt(dc_m[gene,])$value
	vals_melt$site = interaction(vecsplit(rownames(vals_melt), "@", 1),
		factor(vecsplit(rownames(vals_melt), "@", 3), levels = levels(pic_comb)), sep = "@")
	vals_melt$site = factor(vals_melt$site, levels = names(which(table(vals_melt$site) > 0)))
	vals_melt$patient = vecsplit(rownames(vals_melt), "@", 2)
	vals_melt$x = as.numeric(factor(vals_melt$site)) * 2
	vals_melt$exp_x = vals_melt$x + runif(nrow(vals_melt), -0.3,0.3)
	vals_melt$obs_x = vals_melt$x + 1 + runif(nrow(vals_melt), -0.3,0.3)
	vals_melt$contrib = with(vals_melt, log2((reg + dc_value) / (reg + t_value)))
	zlim = max(abs(vals_melt$contrib))
	vals_melt$col_v=round(100 * with(vals_melt, (contrib + zlim) / (2 * zlim))) + 1
	vals_melt$col = grad[ vals_melt$col_v]
	exp_med = with(vals_melt, tapply(exp_value, site, median))
	obs_med = with(vals_melt, tapply(value, site, median))
	ylim = c(0,max(c(vals_melt$exp_value, vals_melt$value)))
	with(vals_melt, plot(1,1, type="n", xlim = round(quantile(x, c(0,1))) + c(-0.5,1.5), ylim = ylim, axes=F, main = gene, cex.main=1, xlab = "", ylab = ""))
	with(vals_melt, segments(exp_x, exp_value, obs_x, value))
	with(vals_melt, points(exp_x, exp_value, cex=2, pch = 21, bg = col))
	with(vals_melt, points(obs_x, value, cex = 2, pch = 21, bg = "gray"))
	segments(seq_along(exp_med) * 2 - 0.4, exp_med, seq_along(exp_med) * 2 + 0.4, lwd = 1, col = "blue")
        segments(seq_along(obs_med) * 2 + 0.6, obs_med, seq_along(obs_med) * 2 + 1.4, lwd = 1, col = "blue")
	axis(2, las=2); #axis(1, at = seq_along(levels(vals_melt$site)) * 2 + 0.5, labels = names(table(vals_melt$site)), las=2)
}
dev.off()

##################

mel_id = "melanoma"
breast_id = "breast"

mel_cl = scdb_mc(mel_id)
breast_cl = scdb_mc(breast_id)

mel_mc = which(mel_cl@colors !="white")
breast_mc = which(breast_cl@colors != "white")

mel_cells = names(mel_cl@mc)[ mel_cl@mc %in% mel_mc]
breast_cells = names(breast_cl@mc)[ breast_cl@mc %in% breast_mc]

mel_umis = read_large_umis(mel_id, cells = mel_cells)
breast_umis = read_large_umis(breast_id, cells = breast_cells)

mel_n = sweep(mel_umis,2,colSums(mel_umis),"/") * 1000
breast_n = sweep(breast_umis,2,colSums(breast_umis),"/") * 1000

color_scheme = unique(rbind(sin_cl@color_key, mel_cl@color_key, breast_cl@color_key))
color_scheme = color_scheme[!duplicated( color_scheme$color),]
color2name = as.vector(color_scheme$group); names(color2name) = color_scheme$color
name2color = as.vector(color_scheme$color); names(name2color) = color_scheme$group
sin_names = color2name[ sin_cl@colors[ sin_cl@mc]]; names(sin_names) = names(sin_cl@mc)
mel_names = color2name[ mel_cl@colors[ mel_cl@mc]]; names(mel_names) = names(mel_cl@mc)
breast_names = color2name[ breast_cl@colors[ breast_cl@mc]]; names(breast_names) = names(breast_cl@mc)

lung_tht2 = colnames(lfp)[ lfp["NMB",] > 1]
t_ord = c(lin_ord[1:9], "Tht_2", lin_ord[10], "IFN")
shared = intersect(rownames(mel_umis), intersect(rownames(breast_umis), rownames(sin_umis)))
all_n = cbind(sin_n[shared,], mel_n[shared,], breast_n[shared,])
all_names = c(sin_names, mel_names, breast_names)[ colnames(all_n)]
all_names[ names(which(sin_cl@mc == lung_tht2))] = "Tht_2"
t_cells_all = names(all_names)[ all_names %in% t_ord]
cell_source = rep(c("lung", "melanoma", "breast"), sapply(list(sin_n, mel_n, breast_n), ncol))
names(cell_source) = colnames(all_n)
comb = paste0(cell_source, "@", all_names); names(comb) = names(cell_source)

m_all = t(apply(all_n[,t_cells_all], 1, tapply, comb[ t_cells_all], mean))

comb2 = paste0(cell_source, "@", ifelse(all_names %in% c("Tht_1", "Tht_2"), "Tht", "other")); names(comb2) = names(cell_source)
m_tht = t(apply(all_n[,t_cells_all], 1, tapply, comb2[ t_cells_all], mean))

###########

#combined_names = all_names

bulk2source = factor(vecsplit(colnames(m_all), "@", 1), levels = c("lung", "breast", "melanoma"))
bulk2pop = factor(vecsplit(colnames(m_all), "@", 2), levels = t_ord)

tht2source = factor(vecsplit(colnames(m_tht), "@", 1), levels = c("lung", "breast", "melanoma"))
tht2pop = vecsplit(colnames(m_tht), "@", 2)

reg = 0.02
sources = names(table(bulk2source))
Z_between = log2((reg + m_tht[,paste0(sources, "@Tht")]) / (reg + m_tht[,paste0(sources, "@other")]))
Z_within = log2((reg + m_all[,paste0(sources, "@Tht_2")]) / (reg + m_all[,paste0(sources, "@Tht")]))
colnames(Z_between) = sources; colnames(Z_within) = sources

tht_genes = names(which(rowSums(Z_between > 1) == 3))
tht2_genes = names(which(rowSums(Z_within > 0.7) == 3))
tht1_genes = names(which(rowSums(Z_within < -0.7) == 3))

#########

all_cells = c(t_cells, good_pics)
t_ord = intersect(lin_ord, names(table(parser_t)))
all_comb = c(paste0(sin_names[t_cells], "@A"), paste0(parser_t[good_pics], "@B")); 
all_comb = factor(all_comb, levels = paste0(rep(t_ord, each=3), "@", rep(c("A", "B", "C"), length(t_ord))))
names(all_comb) = all_cells
good_cells = intersect(colnames(ds), all_cells)

cols = rep("white", length(levels(all_comb))); cols[ seq(2, length(cols), by=3)] = name2color[t_ord]
borders = rep("black", length(levels(all_comb))); borders[ seq(1, length(borders), by=3)] = name2color[t_ord]
pdf(paste0(outdir, "/Fig3b.pdf"), useDingbats=F)
par(mar = c(10,3,1,1))
boxplot(colSums(log(1 + ds[ tht_genes, good_cells])) ~ all_comb[ good_cells], las=2, col=cols, border=borders, axes=F)
axis(2, las=2); axis(1, at = c(seq(1.5, length(levels(all_comb[good_cells])), by=3)), labels = t_ord, las=2)
dev.off()

###########

t_nms = union(rev(read.table(paste0("figures/figure1/Fig1d.txt"), stringsAsFactor=F)[[1]]), c("IFIT1", "IRF7", "STAT1"))

breast_lfp = log2(breast_cl@mc_fp)
IM = breast_lfp[t_nms, ]
vct = factor(color2name[ breast_cl@colors[ as.numeric(colnames(IM))]], levels = t_ord); names(vct) = colnames(IM)
IM = IM[,names(sort(vct))]

pdf(paste0(supdir, "/FigS3a_right.pdf"), useDingbats=F, height=10, width=7)
par(mar = c(0.5,5,0.5,0.5), fig = c(0,1,0.1,1))
image.2(IM, balance=T, vct=vct[ colnames(IM)], annotate="rows"); box()
par(fig = c(0,1,0,0.1), new=T)
image(matrix(seq_len(ncol(IM))), axes=F, col = breast_cl@colors[ as.numeric(colnames(IM))]); box()
dev.off()

mel_lfp = log2(mel_cl@mc_fp)
IM = mel_lfp[t_nms, ]
vct = factor(color2name[ mel_cl@colors[ as.numeric(colnames(IM))]], levels = t_ord); names(vct) = colnames(IM)
IM = IM[,names(sort(vct))]

pdf(paste0(supdir, "/Fig3a_left.pdf"), useDingbats=F, height=10, width=7)
par(mar = c(0.5,5,0.5,0.5), fig = c(0,1,0.1,1))
image.2(IM, balance=T, vct=vct[ colnames(IM)], annotate="rows"); box()
par(fig = c(0,1,0,0.1), new=T)
image(matrix(seq_len(ncol(IM))), axes=F, col = mel_cl@colors[ as.numeric(colnames(IM))]); box()
dev.off()

#########

tht_genes = c("SLAMF6", "IL21", "MAF", "ICOS", "BATF", "PDCD1", "SH2D1A", "TIAM1",
	"BCL6", "CXCR5", "CXCL13", "IFNG", "NR3C1", "BTLA", "TIGIT", "CD200", "ZBED2",
	"CD40LG", "CD4", "PRDM1", "IGFL2", "CTLA4", "CXCR3", "TNFRSF4",
	"BHLHE40", "TNFRSF9", "TNFRSF18", "RGS1", "CD2", "RBPJ", "CD7", "SLA","TNFAIP3", "NFKBIA", "PKM")
other_genes = c("CCR7", "IL7R", "TCF7", "FOXO1", "SOX4", "TXNIP", "FOXP3", "IL2RA", "IL2RB", "LAG3", "GZMA", "GZMB", "HAVCR2", "ITGAE")
bad_genes = grep("^RPL|^RPS|^MT[A-Z]|^HLA-|^AC[0-9]|^SNOR|^HIST[0-9]|^HSP[0-9]|^MIR[0-9]", rownames(umis), v=T)

important_genes = union(tht_genes, other_genes)
Z_lung = log2((0.02 + m_tht[,"lung@Tht"]) / (0.02 + m_all[,setdiff(grep("lung", colnames(m_all), v=T), c("lung@Tht_1", "lung@Tht_2"))]))
Z_mel = log2((0.02 + m_tht[,"melanoma@Tht"]) / (0.02 + m_all[,setdiff(grep("melanoma", colnames(m_all), v=T), c("melanoma@Tht_1", "melanoma@Tht_2"))]))
Z_breast = log2((0.02 + m_tht[,"breast@Tht"]) / (0.02 + m_all[,setdiff(grep("breast", colnames(m_all), v=T), c("breast@Tht_1", "breast@Tht_2"))]))

good_genes = setdiff(shared, bad_genes)
genes = union(important_genes, names(which(apply(abs(cbind(Z_breast[good_genes,], Z_mel[good_genes,], Z_lung[good_genes,])),1,function(x) sort(x,T)[2]) > 1.5)))


medians = t(apply(m_all, 1, tapply, bulk2source, median))
IM = log((0.02 + m_all[genes,]) / (0.02 + medians[genes, as.vector(bulk2source)]))
IM = IM[, order(bulk2pop, bulk2source)]
IM = IM[ order(max.col(IM)),]

k = nrow(IM) / 20
data = as.data.frame(IM)
data$id = rownames(data)
data = data[,c(ncol(data), 1:(ncol(data) - 1))]
km <- TGL_kmeans_tidy(data, k=k, metric='euclid', verbose=TRUE, seed = 18)
#gc = km$cluster$clust; names(gc) = km$cluster$id
#cls = cumsum(table(gc)) / length(gc)
centers = as.matrix(km$centers[,-1]); rownames(centers) = seq_len(k)
centers = centers[ order(max.col(centers)),]
km_clusts = as.numeric(factor(km$cluster$clust, levels = rownames(centers))); names(km_clusts) = rownames(IM)

bad_clusts = c()
disp_genes = setdiff(genes, names(km_clusts)[ km_clusts %in% bad_clusts])
sources = names(table(bulk2source))
pdf(paste0(supdir, "/FigS2b.pdf"), height=10, width=12)
par(mar = c(1,3,1,1))
for (i in seq_along(sources)) {
	s = sources[i]
	par(fig = c((i - 1) / 3, i/3, 0.1, 1), new = (i > 1))
	pops = grep(paste0(s,"@"), colnames(IM), v=T)
	image.2(IM[disp_genes,pops], balance=T, hct = km_clusts[disp_genes], zlim = c(-1,1) * max(abs(IM)), annotate="none"); box()
	par(fig = c((i - 1) / 3, i/3, 0,0.1), new = T)
	image(matrix(seq_along(pops)), axes = F, col = name2color[ vecsplit(pops, "@", 2)]); box()
}
dev.off()

########

indir = "output/published_data/breast/"

t_ord = setdiff(c(lin_ord[1:9], "Tht_2", lin_ord[10], "IFN"), "T_CD4_memory")
breast_tcr = rbind(read.delim(paste0(indir, "/1879-BIOKEY_barcodes_vdj_combined_cohort1.csv"), stringsAsFactor=F, sep = ","),
	read.delim(paste0(indir, "/1880-BIOKEY_barcodes_vdj_combined_cohort2.csv"), stringsAsFactor=F, sep = ","))
breast_tcr = breast_tcr[ breast_tcr$barcode %in% names(which(table(breast_tcr$barcode) == 1)),]
rownames(breast_tcr) = paste0(breast_tcr$barcode, "-1")

breast_stats = scdb_mat(breast_id)@cell_metadata
tcr_cells = intersect(breast_cells, rownames(breast_tcr))
cell2patient = breast_stats[tcr_cells, "patient_id"]; names(cell2patient) = tcr_cells
cell2clone_id = breast_tcr[tcr_cells, "cdr3_nt"]; names(cell2clone_id) = tcr_cells
cell2tp = breast_stats[tcr_cells, "timepoint"]; names(cell2tp) = tcr_cells

#########

cell2comb = paste0(cell2patient, "@", cell2tp, "@", cell2clone_id); names(cell2comb) = tcr_cells
clones = names(which(table(cell2comb) > 1))

pre_cells = intersect(tcr_cells, rownames(breast_stats)[ breast_stats$timepoint == "Pre"])
comb = paste0(breast_names[pre_cells], "@", cell2patient[ pre_cells]); names(comb) = pre_cells
X = table(comb[ pre_cells], cell2comb[pre_cells] %in% clones)
X = X[ rowSums(X) > 2,]
Xn = X / rowSums(X)
samp2pop = factor(vecsplit(rownames(X), "@", 1), levels = t_ord)
samp2patient = vecsplit(rownames(X), "@", 2)
good_patients = names(which(table(samp2patient) >= 7))
X_small = Xn[ samp2patient %in% good_patients,]

x = as.numeric(samp2pop[rownames(Xn)]) + runif(nrow(Xn), -0.1, 0.1); y = Xn[,2]
meds = tapply(y, as.numeric(samp2pop[ names(y)]), median)
pdf(paste0(outdir, "/Fig3g.pdf"), useDingbats=F) #, height=1000, width=1000)
plot(x,y, axes=F, xlab = "", ylab = "", pch = 21, cex = 2, bg = name2color[ as.vector(samp2pop)], lwd=1)
segments(seq_along(meds)-0.2, meds, seq_along(meds) + 0.2, lwd=4, col="red")
axis(2, las=2); axis(1, at = seq_along(levels(samp2pop)), labels = levels(samp2pop), las=2)
dev.off()

clone_cells = pre_cells[ cell2comb[pre_cells] %in% clones]
k=1e4
cells = rep(NA,k); neighbors = cells; s_neighbors = cells
for (i in seq_len(k)) {
	cells[i] = sample(clone_cells, 1)
	patient = cell2patient[cells[i]]
	clone = cell2comb[cells[i]]
	neighbors[i] = sample(setdiff(names(which(cell2comb == clone)), cells[i]), 1)
	s_neighbors[i] = sample(setdiff(names(which(cell2patient == patient)), cells[i]), 1)	
}

fac = factor(breast_names, levels = t_ord); names(fac) = names(breast_names)
true_share = table(fac[cells], fac[neighbors])
true_share = true_share + t(true_share)
shuf_share = table(fac[cells], fac[s_neighbors])
shuf_share = shuf_share + t(shuf_share)
IM = log2((1 + true_share) / (1 + shuf_share))
IM[is.infinite(IM)] = NA

pdf(paste0(outdir, "/Fig3h.pdf"), useDingbats=F)
par(mar = c(5,5,0.5,0.5))
image.2(IM, balance=T, annotate="both"); box()
dev.off()

############

mel_tcr = read.delim("annotations/melanoma_tcr-seq.txt", stringsAsFactor=F)
tcr_cells = intersect(mel_cells, rownames(mel_tcr))
cell2patient = mel_tcr[tcr_cells, "Patient"]; names(cell2patient) = tcr_cells
cell2clone_id = mel_tcr[tcr_cells, "clone_id"]; names(cell2clone_id) = tcr_cells

cell2comb = paste0(cell2patient, "@", cell2clone_id); names(cell2comb) = tcr_cells
clones = names(which(table(cell2comb) > 1))

comb = paste0(mel_names[tcr_cells], "@", cell2patient); names(comb) = tcr_cells
X = table(comb[ tcr_cells], cell2comb[tcr_cells] %in% clones)
X = X[ rowSums(X) > 2,]
Xn = X / rowSums(X)
samp2pop = factor(vecsplit(rownames(X), "@", 1), levels = t_ord)
samp2patient = vecsplit(rownames(X), "@", 2)
#good_patients = names(which(table(samp2patient) == max(table(samp2patient))))
good_patients = names(which(table(samp2patient) >= 8))
X_small = Xn[ samp2patient %in% good_patients,]

x = as.numeric(samp2pop[rownames(Xn)]) + runif(nrow(Xn), -0.1, 0.1); y = Xn[,2]
meds = tapply(y, as.numeric(samp2pop[ names(y)]), median)
pdf(paste0(supdir, "/Fig4e.pdf"), useDingbats=F) #, height=1000, width=1000)
plot(x,y, axes=F, xlab = "", ylab = "", pch = 21, cex = 2, bg = name2color[ as.vector(samp2pop)], lwd=1)
segments(seq_along(meds)-0.2, meds, seq_along(meds) + 0.2, lwd=4, col="red")
axis(2, las=2); axis(1, at = seq_along(levels(samp2pop)), labels = levels(samp2pop), las=2)
dev.off()

clone_cells = tcr_cells[ cell2comb[tcr_cells] %in% clones]
k=1e4
cells = rep(NA,k); neighbors = cells; s_neighbors = cells
for (i in seq_len(k)) {
        cells[i] = sample(clone_cells, 1)
        patient = cell2patient[cells[i]]
        clone = cell2comb[cells[i]]
        neighbors[i] = sample(setdiff(names(which(cell2comb == clone)), cells[i]), 1)
        s_neighbors[i] = sample(setdiff(names(which(cell2patient == patient)), cells[i]), 1)
}

fac = factor(mel_names, levels = t_ord); names(fac) = names(mel_names)
true_share = table(fac[cells], fac[neighbors])
true_share = true_share + t(true_share)
shuf_share = table(fac[cells], fac[s_neighbors])
shuf_share = shuf_share + t(shuf_share)
IM = log2((1 + true_share) / (1 + shuf_share))
IM[is.infinite(IM)] = NA

pdf(paste0(supdir, "/FigS4f.pdf"), useDingbats=F)
par(mar = c(5,5,0.5,0.5))
image.2(IM, balance=T, annotate="both"); box()
dev.off()

############

tumor_t = intersect(t_cells, rownames(sin_stats)[ sin_stats$site == "tumor"])
disp_genes = c("CXCL13", "IL21", "CCDC50", "ZBED2", "PTPN13", "BHLHE40",
	"LAG3", "IFNG", "CCL3", "HAVCR2", "PDCD1", "TIGIT",
	"BTLA", "IGFL2", "NMB", "NR3C1", "MAF", "CTLA4")

pops = intersect(lin_ord, names(table(sin_names[ tumor_t])))
fac = factor(sin_names[tumor_t], levels = pops); names(fac) = tumor_t
m = t(apply(sin_n[disp_genes, tumor_t],1,tapply, fac, sum))
n = tapply(colSums(sin_n[,tumor_t]),fac,sum)
Y = sweep(t(apply(m,1,binconf,n)), 2, rep(tapply(colSums(sin_n[,tumor_t]), fac, mean), 3), "*")
colnames(Y) = paste0(rep(pops,3), rep(c("", "-", "+"), each=length(pops)))
pdf(paste0(outdir, "/Fig3a.pdf"), useDingbats=F, height=7, width=14)
par(mfrow=c(3,6), mar = c(3,3,1,1))
for (gene in disp_genes) {
	X = barplot(Y[gene, pops], col = name2color[pops], ylim = c(0,max(Y[gene,]) * 1.05), main = gene, las=2, names.arg = rep("", length(pops)))
	obs_coords = X
	ci.l = Y[gene, paste0(pops, "-")]; ci.u = Y[gene, paste0(pops, "+")]
	segments(obs_coords, ci.l, y1 = ci.u);
	segments(obs_coords-0.2, ci.l, x1 = obs_coords + 0.2); segments(obs_coords-0.2, ci.u, x1 = obs_coords + 0.2);
}
dev.off()

other_genes = setdiff(tht_genes, disp_genes)

pops = intersect(lin_ord, names(table(sin_names[ tumor_t])))
fac = factor(sin_names[tumor_t], levels = pops); names(fac) = tumor_t
m = t(apply(sin_n[other_genes, tumor_t],1,tapply, fac, sum))
n = tapply(colSums(sin_n[,tumor_t]),fac,sum)
Y = sweep(t(apply(m,1,binconf,n)), 2, rep(tapply(colSums(sin_n[,tumor_t]), fac, mean), 3), "*")
colnames(Y) = paste0(rep(pops,3), rep(c("", "-", "+"), each=length(pops)))
pdf(paste0(supdir, "/FigS3c.pdf"), useDingbats=F, height=7, width=16)
par(mfrow=c(3,7), mar = c(3,3,1,1))
for (gene in other_genes) {
        #png(paste0(outdir, "/lung_genes/", gene, ".png"), height=1000, width=1000)
        X = barplot(Y[gene, pops], col = name2color[pops], ylim = c(0,max(Y[gene,]) * 1.05), main = gene, las=2, names.arg = rep("", length(pops)))
        obs_coords = X
        ci.l = Y[gene, paste0(pops, "-")]; ci.u = Y[gene, paste0(pops, "+")]
        segments(obs_coords, ci.l, y1 = ci.u);
	segments(obs_coords-0.2, ci.l, x1 = obs_coords + 0.2); segments(obs_coords-0.2, ci.u, x1 = obs_coords + 0.2);
}
dev.off()

##############

hybrid_cells = c(dc_cells, "WMC1609889")
thts = names(breast_names)[ breast_names %in% c("Tht", "Tht_2")]
hybrid_umis = cbind(breast_umis[shared, thts], sin_umis[shared, hybrid_cells])

hybrid_clusts = c(paste0("T", breast_cl@mc[ thts]),
	paste0("DC", sin_cl@mc[ dc_cells]))
hybrid_levels = c(paste0("T", sort(unique(breast_cl@mc[ thts]))), paste0("DC", sort(unique(sin_cl@mc[ dc_cells]))))
clusts = as.numeric(factor(hybrid_clusts, levels = hybrid_levels)); names(clusts) = c(thts, dc_cells)

cell_metadata = data.frame(group = rep(c("T", "DC"), c(length(thts), length(hybrid_cells)))); rownames(cell_metadata) = colnames(hybrid_umis)
scdb_add_mat("hybrid", tgScMat(mat=hybrid_umis, cell_metadata=cell_metadata))
scdb_add_mc("hybrid", tgMCCov(clusts, "WMC1609889", scdb_mat("hybrid")))

hybrid_id = "hybrid"
thts = names(breast_names)[ breast_names %in% c("Tht", "Tht_2")]
lr_features =  choose_lr_features(hybrid_id, thts, dc_cells, bad_genes, must_haves = names(scdb_gset(id)@gene_set))
mle_features = choose_mle_features(hybrid_id, hybrid_id, thts, dc_cells, union(bad_markers, bad_genes), existing_list= names(scdb_gset(id)@gene_set))

mle_features = union(mle_features, "CXCL13")

numis=800
tht_pics = names(parser_t)[ parser_t %in% c("Tht")]
hybrid_res = run_pic_seq(hybrid_id, hybrid_id, ds[, tht_pics], thts, dc_cells, intersect(lr_features, shared), intersect(mle_features, shared), paste0(supdir, "/alpha_estimation_hybrid.png"),
        numis, reg = 1e-4)
hybrid_res$well = rownames(hybrid_res)

#hybrid_t_mc = mle_res[good_pics, "a_mc"]; names(t_mc) = good_pics
hybrid_t = color2name[ breast_cl@colors[  as.numeric(gsub("T", "", hybrid_levels[hybrid_res$a_mc]))]]; names(hybrid_t) = tht_pics

#############

tht_pics = names(parser_t)[ parser_t %in% c("Tht")]
db_clusts = factor(hybrid_t)
db_vec = as.vector(db_clusts); names(db_vec) = names(hybrid_t)


nms = rev(read.table(paste0(supdir, "FigS3d.txt"), stringsAsFactor=F, h=T)[[1]])
pdf(paste0(supdir, "/FigS3d.pdf"), useDingbats=F, height=10, width=28)
par(mar = c(1,5,1,1), fig = c(0,0.25,0.1,1))
message("breast:")
breast_thts = names(breast_names)[ breast_names %in% c("Tht", "Tht_2")]
sin_vec = breast_names[ breast_thts]; names(sin_vec) = breast_thts
sin_ord = plot_sc_heatmap(breast_id, breast_id, nms, clusts = sin_vec, good_clusts = names(table(sin_vec)), cells = breast_thts,
        annotate=T, annotate_what="rows", normalize=T, lty=1, draw_cls=T); box(lwd=1)
par(fig=c(0,0.25,0,0.1), new=T)
image(matrix(seq_along(sin_ord)), axes = F, col = name2color[ breast_names[ sin_ord]]); box(lwd=1)
par(fig = c(0.25,0.5,0.1,1), new=T)
mel_thts = names(mel_names)[ mel_names %in% c("Tht", "Tht_2")]
sin_vec = all_names[mel_thts]; names(sin_vec) = mel_thts
message("mleanoma:")
sin_ord = plot_sc_heatmap(mel_id, mel_id, nms, clusts = sin_vec, good_clusts = names(table(sin_vec)), cells = mel_thts,
        annotate=T, annotate_what="rows", normalize=T, lty=1, draw_cls=T); box(lwd=1)
par(fig=c(0.25,0.5,0,0.1), new=T)
image(matrix(seq_along(sin_ord)), axes = F, col = name2color[ all_names[ sin_ord]]); box(lwd=1)
par(fig = c(0.5,0.75,0.1,1), new=T)
lung_thts = names(which(sin_names == "Tht"))
sin_vec = all_names[lung_thts]; names(sin_vec) = lung_thts
message("lung:")
sin_ord = plot_sc_heatmap(id, id, nms, clusts = sin_vec, good_clusts = names(table(sin_vec)), cells = lung_thts,
        annotate=T, annotate_what="rows", normalize=T, lty=1, draw_cls=T); box(lwd=1)
par(fig=c(0.5,0.75,0,0.1), new=T)
image(matrix(seq_along(sin_ord)), axes = F, col = name2color[ all_names[ sin_ord]]); box(lwd=1)
par(fig = c(0.75,1,0.1,1), new=T)
message("PICs:")
db_ord = plot_sc_heatmap(id_d, id_d, nms, clusts = db_vec, cells = tht_pics, good_clusts = names(table(db_clusts)),
        annotate=T, annotate_what="rows", normalize=T, lty=1, draw_cls=T); box(lwd=1)
par(fig = c(0.75,1,0,0.1), new=T)
image(matrix(seq_along(db_ord)), axes = F, col = name2color[ hybrid_t[ db_ord]]); box(lwd=1)
dev.off()

db_tht_names = rep("Other", length(parser_t)); names(db_tht_names) = names(parser_t)
db_tht_names[ names(hybrid_t)] = hybrid_t

sin_tht_names = all_names[ t_cells]
sin_tht_names[ sin_names[t_cells] != "Tht"] = "Other"

comb = with(cell_stats, paste0(site, "@", patient, "@", organ, ":", gate)); names(comb) = rownames(cell_stats)
t_dist = rbind(table(comb[good_pics],  db_tht_names[good_pics]),
        table(comb[t_cells],  sin_tht_names[t_cells]))
t_dist = t_dist[ rowSums(t_dist) >= 5,]
t_dist = t_dist[ -grep(":APC|normal@", rownames(t_dist)),]
t_dist = t_dist[ order(vecsplit(rownames(t_dist), ":", 1), factor(vecsplit(rownames(t_dist), ":", 2), levels = c("T", "APC", "PIC"))),]
desc2rep = vecsplit(rownames(t_dist), ":", 1)
t_dist = t_dist[ sapply(names(which(table(desc2rep) == 2)), grep, rownames(t_dist), v=T), ]
t_n = t_dist / rowSums(t_dist)
desc2rep = desc2rep[ rownames(t_dist)]
t_melt = melt(t_n[grep(":T", rownames(t_n)),-1])
t_melt$pic_value = melt(t_n[grep(":PIC", rownames(t_n)),-1])$value
t_melt$x = as.numeric(t_melt$Var2) * 2

pdf(paste0(supdir, "/FigS3e.pdf"), useDingbat=F)
plot(1,1, type="n", xlim = c(1,6), ylim = c(0,max(c(t_melt$value, t_melt$pic_value))), axes=F)
with(t_melt, segments(x, value, x + 1, pic_value))
with(t_melt, points(x, value, pch=21, cex=2, bg = name2color[ as.vector(Var2)]))
with(t_melt, points(x + 1, pic_value, pch=21, cex=2, bg = name2color[ as.vector(Var2)]))
axis(2, las=2); axis(1, at = c(2,3, 4,5), labels = rep("", 4))
pvals = with(t_melt, sapply(names(table(Var2)), function(x) wilcox.test(value[Var2 == x], pic_value[Var2 == x], paired=T)$p.value))
dev.off()
