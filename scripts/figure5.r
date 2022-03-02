message("Generating Fig. 5 and S5")

outdir = "figures/figure5/"
dir.create(outdir)

supdir = "figures/figureS5/"
dir.create(supdir)

id = "ot2"
sin_2d = scdb_mc2d(id); sin_cl = scdb_mc(id); sin_mat = scdb_mat(id)

cells = names(sin_cl@mc)
umis = read_large_umis(id, cells = cells)
fp = sin_cl@mc_fp
lfp = log2(sin_cl@mc_fp)

sin_stats = sin_mat@cell_metadata[cells,]
sin_stats[ sin_stats$sorting.scheme == "Cd11c+", "sorting.scheme"] = "CD11c+"
sin_comb = with(sin_stats, paste0(sorting.scheme, "@", PIC, "@", tissue, "@", timepoint, "@",  as.numeric(factor(date)), "-", replicate))
names(sin_comb) = rownames(sin_stats)
umis_n = sweep(umis,2,colSums(umis),"/") * 1000

color_scheme = sin_cl@color_key
color2name = as.vector(color_scheme$group); names(color2name) = color_scheme$color
name2color = as.vector(color_scheme$color); names(name2color) = color_scheme$group
sin_names = color2name[ sin_cl@colors[ sin_cl@mc]]; names(sin_names) = names(sin_cl@mc)

lin_ord = c("T_naive", "CD8_T", "T_ot2", "T_ot2_act", "T_Vim", "CD8_T_Ccl5")

umicount = colSums(umis)

bad_mc = c() 
clust_ord = setdiff(order(factor(color2name[ sin_cl@colors], levels = lin_ord)), bad_mc)

t_pops = lin_ord
t_clusts = intersect(clust_ord, which(color2name[ sin_cl@colors] %in% t_pops))
t_cells = setdiff(names(sin_cl@mc)[ sin_cl@mc %in% t_clusts], names(which(umicount > 1e4)))

############
#cluster LN and TME

activation_genes = read.table("figures/figure4/activation_genes.txt", stringsAsFactors=F)[[1]]

sin_comb = with(sin_stats, paste0(sorting.scheme, "@", PIC, "@", tissue, "@", treatment, "@", timepoint, ":",
        as.numeric(factor(date)), "-", replicate))
names(sin_comb) = rownames(sin_stats)

rel_cells = intersect(t_cells, rownames(sin_stats)[ sin_stats$sorting.scheme %in% c("TCRb+", "TCRb+CD45.1+")])
bad_genes = grep("^Rpl|^Rps|^Gm[0-9]|^AC[0-9]|^Snor|^His[0-9]|^Hsp[0-9]|^Trim30|^mt-|^Ig", rownames(umis), v=T)
bad_samps = union("TCRb+CD45.1+@Singlets@dLN@Ova-Vaccination@d10:4-2", names(which(table(sin_comb[ t_cells]) < 20)))
m = t(apply(umis_n[ setdiff(rownames(umis_n), bad_genes), rel_cells], 1, tapply, sin_comb[rel_cells], mean))
m = m[, setdiff(colnames(m), bad_samps)]

samp_names = vecsplit(colnames(m), ":", 1)
comb_names = vecsplit(sin_comb, ":", 1)
samp2sort = factor(vecsplit(samp_names, "@", 1))
samp2tissue = factor(vecsplit(samp_names, "@", 3))
samp2tp = factor(vecsplit(samp_names, "@", 5))
samp2treat = factor(vecsplit(samp_names, "@", 4))

m_med = t(apply(m, 1, tapply, samp_names[ colnames(m)], median))

good_genes = setdiff(rownames(m), c(bad_genes))
samp_ord = c("TCRb+CD45.1+@Singlets@cLN@Ova-Vaccination@d10", "TCRb+CD45.1+@Singlets@dLN@Ova-Vaccination@d10",
	"TCRb+@Singlets@tumor@Ova-Vaccination@d10", "TCRb+@Singlets@tumor@Ova-Vaccination@d17", "TCRb+CD45.1+@Singlets@tumor@Ova-Vaccination@d10",
	"TCRb+CD45.1+@Singlets@tumor@Ova-Vaccination@d17")
IM = log2((0.05 + m_med[,samp_ord]) / (0.05 + apply(m_med[,samp_ord], 1, median)))
nms = unique(as.vector(apply(IM[good_genes,] ,2,function(x){ head(names(sort(x, decreasing = T)),20)})))
nms = union(intersect(nms, names(which(apply(IM,1, function(x) sort(x,T)[2]) > 1))), activation_genes)
IM2 = IM[nms,]
k = nrow(IM2) / 8
data = as.data.frame(IM2)
data$id = rownames(data)
data = data[,c(ncol(data), 1:(ncol(data) - 1))]
km <- TGL_kmeans_tidy(data, k=k, metric='euclid', verbose=TRUE, seed = 18)
centers = as.matrix(km$centers[,-1]); rownames(centers) = seq_len(k)
centers = centers[ order(max.col(centers)),]
km_clusts = as.numeric(factor(km$cluster$clust, levels = rownames(centers))); names(km_clusts) = rownames(IM2)

meta_mat = rbind(samp2sort, samp2tissue, samp2tp) + c(0,2,5)
cols = c("gray80", "gray20", "#E5E4D1", "#E5DE61", "#CD9D57", "powderblue", "pink")
rownames(meta_mat) = c("Gate", "Tissue", "Timepoint")
unique_meta = t(unique(t(meta_mat)))
colnames(unique_meta) = vecsplit(colnames(unique_meta), ":", 1)

pdf(paste0(outdir, "/Fig5g.pdf"), height=15, width=10)
par(fig = c(0,1,0.1,1), mar = c(1,10,1,5))
image.2(IM2, balance=T, hct = km_clusts, annotate="both"); box()
par(fig = c(0,1,0,0.1), new=T)
image.2(unique_meta[, colnames(IM2)], col = cols, annotate="rows"); box()
dev.off()

##########


score = colSums(log(1 + 7 * umis_n[activation_genes, t_cells]))
sin_group = with(sin_stats, paste0(sorting.scheme, "@", tissue, "@", timepoint))
sin_rep = with(sin_stats, paste0(treatment, "-", as.numeric(factor(date)), "-", replicate))
names(sin_group) = rownames(sin_stats); names(sin_rep) = rownames(sin_stats)
sub_cells = intersect(t_cells, rownames(sin_stats)[ sin_stats$sorting.scheme %in% c("TCRb+", "TCRb+CD45.1+") & 
	sin_stats$tissue %in% c("cLN", "dLN")])
medians = tapply(score[sub_cells], paste0(sin_group[sub_cells], ":", sin_rep[sub_cells]), median)
medians2group = vecsplit(names(medians), ":", 1)

pairs_x = c("TCRb+@dLN@d10", "TCRb+CD45.1+@dLN@d10", "TCRb+CD45.1+@cLN@d10", "TCRb+CD45.1+@dLN@d10", "TCRb+CD45.1+@dLN@d10")
pairs_y = c("TCRb+@cLN@d10", "TCRb+CD45.1+@cLN@d10", "TCRb+@cLN@d10", "TCRb+@dLN@d10", "TCRb+CD45.1+@dLN@d17")
pvals = rep(NA, length(pairs_x)); names(pvals) = paste0(pairs_x, "-", pairs_y)
for (i in seq_along(pvals)) {
	pvals[i] = wilcox.test(medians[ medians2group == pairs_x[i]], medians[ medians2group	== pairs_y[i]])$p.value
}

medians2rep = gsub("Ova-Vaccination-", "", vecsplit(names(medians), ":", 2))
groups = names(table(medians2group))
levels = seq_len(length(groups)*3/2)
levels[ setdiff(levels, seq(3,length(levels),by=3))] = groups
medians2group = factor(medians2group, levels = levels)
x = as.numeric(medians2group) + runif(length(medians2group), -0.1, 0.1)
pdf(paste0(outdir, "/Fig5b.pdf"), useDingbats=F)
par(mar = c(5,3,1,1))
boxplot(score[sub_cells] ~ factor(sin_group[ sub_cells], levels = levels), las=2, outline=T, col = c("powderblue", "pink", "white"), axes=F)
points(x, medians, cex=2, pch=21, bg = c("blue", "red")[ as.numeric(factor(vecsplit(as.vector(medians2group[names(medians)]), "@", 3)))])
#text(x, medians, medians2rep[ names(medians)])
axis(2, las=2); 
#axis(1, at = setdiff(seq_along(levels), seq(3,length(levels),by=3)), labels = names(which(table(medians2group) > 0)), las=2)
axis(1, at = seq(1.5, length(levels), by=3), labels = rep(c("cLN", "tdLN"), 2))
dev.off()

###############

rel_cells = intersect(t_cells, rownames(sin_stats)[ sin_stats$sorting.scheme == "TCRb+CD45.1+"])# & sin_stats$treatment == "Ova-Vaccination"])
names(sin_group) = names(sin_comb)
m_group = t(apply(umis_n[ setdiff(rownames(umis_n), bad_genes), rel_cells], 1, tapply, sin_group[rel_cells], mean))

samp1 = "TCRb+CD45.1+@dLN@d10"; samp2 = "TCRb+CD45.1+@cLN@d10"
diff_genes = intersect(good_genes, scr_chi_square_diff_genes(umis, g1 = names(which(sin_group == samp1)), g2 = names(which(sin_group == samp2)), pval=1e-3, fdr=T))

##############

virus_id = "virus"
virus_mat = scdb_mat(virus_id)
virus_cl = scdb_mc(virus_id)
virus_stats = virus_mat@cell_metadata[ names(virus_cl@mc),]
virus_comb = with(virus_stats, paste0(sorting.scheme, "@", treatment, ":", replicate)); names(virus_comb) = rownames(virus_stats)
virus_umis = as.matrix(virus_mat@mat)
virus_n = sweep(virus_umis, 2, colSums(virus_umis), "/") * 1000
virus_m = t(apply(virus_n[ rownames(m), names(virus_comb)], 1, tapply, virus_comb, mean))
v_samp_names = vecsplit(colnames(virus_m), ":", 1)
virus_med = t(apply(virus_m, 1, tapply, v_samp_names[ colnames(virus_m)], median))

virus_z = log2((virus_med[,"SMARTA+Tfh+@rVSV-infected"] + 0.02) / (virus_med[,"SMARTA+@not-infected"] + 0.02))
#z = log2((m_group[,"] + 0.02) / (m_group[,1] + 0.02))
cancer_z = log2((m_med[,"TCRb+CD45.1+@Singlets@dLN@Ova-Vaccination@d10"] + 0.02) / (m_med[,"TCRb+CD45.1+@Singlets@cLN@Ova-Vaccination@d10"] + 0.02))
iv_z = as.matrix(read.delim("figures/figure4/tht_vs_naive.txt", stringsAsFactor=F, row.names=1))[,1]

shared = intersect(names(iv_z), names(cancer_z))
IM = cbind(cancer = cancer_z[shared], iv = iv_z[ shared], virus = virus_z[shared])
pos_genes = unlist(apply(IM, 2, function(x) names(head(sort(x,T),50))))
neg_genes = unlist(apply(IM, 2, function(x) names(head(sort(x,F),20))))
IM2 = IM[union(pos_genes, neg_genes),]
k = nrow(IM2) / 10
data = as.data.frame(IM2)
data$id = rownames(data)
data = data[,c(ncol(data), 1:(ncol(data) - 1))]
km <- TGL_kmeans_tidy(data, k=k, metric='euclid', verbose=TRUE, seed = 18)
centers = as.matrix(km$centers[,-1]); rownames(centers) = seq_len(k)
centers = centers[ order(max.col(centers)),]
km_clusts = as.numeric(factor(km$cluster$clust, levels = rownames(centers))); names(km_clusts) = rownames(IM2)
rownames(centers) = seq_len(nrow(centers))

good_clusts = rownames(centers)
intensity = apply(abs(IM2), 1, max)
disp_genes = union(unlist(sapply( as.character(good_clusts), function(i) names(head(sort(intensity[ names(which(km_clusts == i))], T), 3)))), c("Il21", "Lag3"))

xlim = max(abs(IM[disp_genes,])) * c(-1,1)
cols = c("#1D9E02", "#CE13AA", "#89613A"); names(cols) = colnames(IM)
pdf(paste0(outdir, "/Fig5c-e.pdf"), useDingbats=F, height=8, width=12)
par(mfrow= c(1,3))
for (cond in colnames(IM)) {
	barplot(t(IM[ disp_genes, cond]), horiz=T, las=2, xlim = xlim, col = cols[cond], main=cond, axes=F)
	abline(v=0); axis(1); grid(col="black"); 
}
dev.off()

###############

sub_cells = intersect(rel_cells, names(sin_group)[ sin_group %in% "TCRb+CD45.1+@dLN@d10"])
act_cells = names(head(sort(score[sub_cells],T), 400)) #intersect(sub_cells, names(which(score > 20)))

treg_cells = read.table("annotations/ln_tregs.txt", stringsAsFactors=F)[[1]]
treg_umis = read_large_umis("pic-seq_all", cells = intersect(treg_cells, scdb_mat("pic-seq_all")@cells))

all_umis = cbind(umis[,act_cells], treg_umis[ rownames(umis),])
all_n = sweep(all_umis, 2, colSums(all_umis), "/") * 1000
all_m = t(apply(all_n, 1, tapply, colnames(all_n) %in% act_cells, mean))
tht_vs_treg_m = log2((0.02 + all_m[,2]) / (0.02 + all_m[,1]))

naive_cells = names(sin_group)[ sin_group == "TCRb+@dLN@d10"]
us_n = umis_n[,union(act_cells, naive_cells)]
compare_m = t(apply(us_n, 1, tapply, colnames(us_n) %in% act_cells, mean))
tht_vs_naive_m = log2((0.02 + compare_m[,2]) / (0.02 + compare_m[,1]))

###############
# Load human breast data for comparison

breast_id = "breast"
breast_cl = scdb_mc(breast_id)
breast_mc = which(breast_cl@colors != "white")
breast_cells = names(breast_cl@mc)[ breast_cl@mc %in% breast_mc]
breast_umis = read_large_umis(breast_id, cells = breast_cells)
breast_n = sweep(breast_umis,2,colSums(breast_umis),"/") * 1000
breast_color_scheme = breast_cl@color_key
breast_color2name = as.vector(breast_color_scheme$group); names(breast_color2name) = breast_color_scheme$color
breast_names = breast_color2name[ breast_cl@colors[ breast_cl@mc]]; names(breast_names) = names(breast_cl@mc)
breast_names[ breast_names %in% c("Tht_1", "Tht_2")] = "Tht"
breast_t = names(which(!is.na(breast_names)))
breast_m = t(apply(breast_n[,breast_t], 1, tapply, breast_names[ breast_t], mean))

###############

tht_vs_treg_h = log2((breast_m[,"Tht"] + 0.02) / (breast_m[,"T_CD4_reg"] + 0.02))
tht_vs_naive_h = log2((breast_m[,"Tht"] + 0.02) / (breast_m[,"T_naive"] + 0.02))

mouse2human = read.delim("annotations/mouse2human.txt", stringsAsFactor=F)
m2h = mouse2human[ mouse2human[,1] %in% names(tht_vs_naive_m),]
unique_genes = m2h[m2h[,2] %in% intersect(unique(m2h[,2]), names(tht_vs_treg_h)),]
unique_genes = unique_genes[ unique_genes[,1] %in% names(which(table(unique_genes[,1]) ==1 )),]
unique_genes = unique_genes[ unique_genes[,2] %in% names(which(table(unique_genes[,2]) ==1 )),]

x = tht_vs_treg_h[ unique_genes[,2]]; y =  tht_vs_treg_m[ unique_genes[,1]]

disp_human = names(which(tht_vs_treg_h > 1 | tht_vs_treg_h < -2))
disp_mouse = names(which(tht_vs_treg_m > 1 | tht_vs_treg_m < -2))
disp_genes = union(disp_mouse, unique_genes[ unique_genes[,2] %in% disp_human, 1])
df = data.frame(x=x, y=y, color = names(y) %in% disp_genes,
        text = ifelse(names(y) %in% disp_genes, names(y), ""))
df = df[ setdiff(rownames(df), bad_genes),]
IM = as.matrix(df[,1:2]); colnames(IM) = c("human", "mouse")
pos_genes = unlist(apply(IM, 2, function(x) names(head(sort(x,T),50))))
neg_genes = unlist(apply(IM, 2, function(x) names(head(sort(x,F),20))))
IM = IM[union(pos_genes, neg_genes),]
k = nrow(IM) / 10
data = as.data.frame(IM)
data$id = rownames(data)
data = data[,c(ncol(data), 1:(ncol(data) - 1))]
km <- TGL_kmeans_tidy(data, k=k, metric='euclid', verbose=TRUE, seed = 18)
centers = as.matrix(km$centers[,-1]); rownames(centers) = seq_len(k)
centers = centers[ order(max.col(centers)),]
km_clusts = as.numeric(factor(km$cluster$clust, levels = rownames(centers))); names(km_clusts) = rownames(IM)
rownames(centers) = seq_len(nrow(centers))

good_clusts = rownames(centers)
intensity = apply(abs(IM), 1, max)
disp_genes = union(unlist(sapply( as.character(good_clusts), function(i) names(head(sort(intensity[ names(which(km_clusts == i))], T), 4)))), c())
xlim = max(abs(IM[disp_genes,])) * c(-1,1)
cols = c("#41AFB5", "#E2522B"); names(cols) = colnames(IM)
pdf(paste0(outdir, "/Fig5f_left.pdf"), useDingbats=F, height=8, width=4)
barplot(t(IM[ disp_genes, ]), beside=T, horiz=T, las=2, xlim = xlim, col = cols, axes=F)
abline(v=0); axis(1); grid(col="black");
dev.off()

m2h = mouse2human[ mouse2human[,1] %in% names(tht_vs_naive_m),]
unique_genes = m2h[m2h[,2] %in% intersect(unique(m2h[,2]), names(tht_vs_naive_h)),]
unique_genes = unique_genes[ unique_genes[,1] %in% names(which(table(unique_genes[,1]) ==1 )),]
unique_genes = unique_genes[ unique_genes[,2] %in% names(which(table(unique_genes[,2]) ==1 )),]

x = tht_vs_naive_h[ unique_genes[,2]]; y =  tht_vs_naive_m[ unique_genes[,1]]

disp_human = names(which(tht_vs_naive_h > 1 | tht_vs_naive_h < -2))
disp_mouse = names(which(tht_vs_naive_m > 1 | tht_vs_naive_m < -2))
disp_genes = union(disp_mouse, unique_genes[ unique_genes[,2] %in% disp_human, 1])
df = data.frame(x=x, y=y, color = names(y) %in% disp_genes,
        text = ifelse(names(y) %in% disp_genes, names(y), ""))
df = df[ setdiff(rownames(df), bad_genes),]
IM = as.matrix(df[,1:2]); colnames(IM) = c("human", "mouse")
pos_genes = unlist(apply(IM, 2, function(x) names(head(sort(x,T),50))))
neg_genes = unlist(apply(IM, 2, function(x) names(head(sort(x,F),20))))
IM = IM[union(pos_genes, neg_genes),]
k = nrow(IM) / 10
data = as.data.frame(IM)
data$id = rownames(data)
data = data[,c(ncol(data), 1:(ncol(data) - 1))]
km <- TGL_kmeans_tidy(data, k=k, metric='euclid', verbose=TRUE, seed = 18)
centers = as.matrix(km$centers[,-1]); rownames(centers) = seq_len(k)
centers = centers[ order(max.col(centers)),]
km_clusts = as.numeric(factor(km$cluster$clust, levels = rownames(centers))); names(km_clusts) = rownames(IM)
rownames(centers) = seq_len(nrow(centers))

good_clusts = rownames(centers)
intensity = apply(abs(IM), 1, max)
disp_genes = union(unlist(sapply( as.character(good_clusts), function(i) names(head(sort(intensity[ names(which(km_clusts == i))], T), 4)))), c())

xlim = max(abs(IM[disp_genes,])) * c(-1,1)
cols = c("#41AFB5", "#E2522B"); names(cols) = colnames(IM)
pdf(paste0(outdir, "/Fig5f_right.pdf"), useDingbats=F, height=8, width=4)
barplot(t(IM[ disp_genes, ]), beside=T, horiz=T, las=2, xlim = xlim, col = cols, axes=F)
abline(v=0); axis(1); grid(col="black");
dev.off()


##############

good_cells = rownames(sin_stats)[ sin_stats$sorting.scheme == "TCRb+CD45.1+"]
t_modules = read.delim(paste0(supdir, "/t_modules.txt"), stringsAsFactor=F, row.names=1)
t_modules = t_modules[ t_modules$annotation != "",]
ct = factor(t_modules$annotation, levels = c("Naive", "Naive II", "Cell cycle", "Activation I", "Activation II",
	"Tumor I", "Tumor II")); names(ct) = rownames(t_modules)

ds = .downsamp(umis[, good_cells], 500)

bad_genes = grep("^Rpl|^Rps|^Gm[0-9]|^AC[0-9]|^Snor|^His[0-9]|^Hsp[0-9]|^Trim30|^mt-|^Ig", rownames(umis), v=T)
nms = setdiff(names(ct), bad_genes)

C = cor(t(log(1 + ds[nms,intersect(colnames(ds), good_cells)]))); diag(C) = NA
good_nms = nms #names(which(apply(C,1,min,na.rm=T) < -0.04 | apply(C,1,function(x) sort(x,T)[2]) > 0.12))
C2 = C[good_nms, good_nms]
m = t(apply(umis_n[ good_nms, good_cells], 1, tapply, sin_group[good_cells], mean)) * min(table(sin_group[good_cells]))
IM = log2((10 + m) / (10 + rowMedians(m)))

grad = colorRampPalette(c("turquoise3", "white", "violetred"))(1000)
pdf(paste0(supdir, "/FigS5c.pdf"), height=12, width=15)
par(mar = rep(0.5,4), fig = c(0,0.8,0,1))
image.2(C2, balance=T, hct=ct, vct=ct, annotate="rows", col=grad, las=1); box()
par(fig = c(0.8,1,0,1), new=T)
image.2(IM, balance=T, hct=ct, annotate=F); box()
dev.off()

foc = log(1 + 7 * umis_n)
modules = apply(foc[ good_nms, good_cells], 2, tapply, ct, sum)
cell_col = ifelse(vecsplit(sin_group[good_cells], "@", 3) == "d10", "powderblue", "pink")
score = modules["Activation II",good_cells]
cell_ord = order(sin_group[ good_cells], score + runif(length(good_cells), -1, 1) * sqrt(var(score)))
cls = cumsum(table(sin_group[good_cells]))

pdf(paste0(supdir, "/FigS5d.pdf"), height=16, width=12, useDingbats=F)
par(mar = c(3,6,0.5,0.5), mfrow = c(nrow(modules), 1))
invisible(sapply(rev(rownames(modules)), function(x) {
	plot(modules[x, cell_ord], cex=1, pch = 20, col = cell_col[ cell_ord], xaxs="i", yaxs="i", axes=F, xlab="", ylab = x); 
	axis(2, las=2); abline(v=cls); box(lwd=1)}))
dev.off()

#########


lateral_genes = as.matrix(read.delim("annotations/lateral_genes.txt", stringsAsFactor=F, row.names=1))[,1]
mouse_lateral = names(which(lateral_genes == "mouse_lateral"))

sub_cells = names(sin_group)[ sin_group %in% c("TCRb+CD45.1+@dLN@d10")]
cc_genes = intersect(rownames(umis_n), mouse_lateral)
cc_score = colSums(log(1 + 7 * umis_n[cc_genes, t_cells]))

X = table(cut(score[sub_cells], quantile(score[sub_cells]), include.lowest=T), cc_score[ sub_cells] > 32)
Xn = X / rowSums(X)
pdf(paste0(supdir, "/FigS5e.pdf"), useDingbats=F)
barplot(Xn[,2], las=2, col="gray40")
dev.off()

##########

