###################
#
# Figure 6
#
###################

message("Generating Fig. S6")

# PIC-seq

supdir = paste0("figures//figureS6/")
dir.create(supdir)

id = "ot2_ln"

sin_2d = scdb_mc2d(id); sin_cl = scdb_mc(id); sin_mat = scdb_mat(id)

cells = names(sin_cl@mc)
umis = read_large_umis(id, cells = cells)
fp = sin_cl@mc_fp
lfp = log2(sin_cl@mc_fp)

sin_stats = sin_mat@cell_metadata[names(sin_cl@mc),]
sin_stats[ sin_stats$sorting.scheme == "Cd11c+", "sorting.scheme"] = "CD11c+"
sin_comb = with(sin_stats, paste0(sorting.scheme, "@", PIC, "@", tissue, "@", treatment, "@", replicate))
names(sin_comb) = rownames(sin_stats)
umis_n = sweep(umis,2,colSums(umis),"/") * 1000
umicount = colSums(umis)

color_scheme = sin_cl@color_key
color2name = as.vector(color_scheme$group); names(color2name) = color_scheme$color
name2color = as.vector(color_scheme$color); names(name2color) = color_scheme$group
sin_names = color2name[ sin_cl@colors[ sin_cl@mc]]; names(sin_names) = names(sin_cl@mc)

annotations = as.matrix(read.delim("annotations/mouse_annotations.txt", stringsAsFactor=F, row.names=1))[,1]
lin_ord = names(annotations)
bad_mc = colnames(lfp)[lfp["Malat1",] > 4]

clust_ord = setdiff(order(factor(color2name[ sin_cl@colors], levels = lin_ord)), bad_mc)

t_pops = names(which(annotations == "T"))
t_clusts = intersect(clust_ord, which(color2name[ sin_cl@colors] %in% t_pops))
t_cells = setdiff(names(sin_cl@mc)[ sin_cl@mc %in% t_clusts], names(which(umicount > 1e4)))

dc_pops = names(which(annotations == "Myeloid"))
dc_clusts =  intersect(clust_ord, which(color2name[ sin_cl@colors] %in% dc_pops))
dc_cells = setdiff(names(sin_cl@mc)[ sin_cl@mc %in% dc_clusts], names(which(umicount > 1e4)))

cells = union(t_cells, dc_cells)

#####################

bad_genes = grep("^Rpl|^Rps|^Gm[0-9]|^AC[0-9]|^Snor|^His[0-9]|^Hsp[0-9]|^Trim30", rownames(umis), v=T)
bad_markers = unique(c("Malat1" ,bad_genes, mouse_lateral))

must_haves = c("Cd8b1", "Cd8a", "Cd4", "Pdcd1", "Tigit", "Ctla4", "Bhlhe40", "Il21", "Cxcr5", "Prdm1", "Bcl6", "Pou2af1", "Slamf6", "Tiam1",
        "Tpi1", "Pkm", "Tnfrsf4", "Tnfrsf18", "Pdcd4", "Smc4", "Gzma", "Nkg7", "Gzmb", "Tcf7", "Il7r", "Ccr7", "S1pr1", "Sell", "Cd200",
        "Cxcr3", "Icos", "Lgals1", "Top2a", "Pcna", "Tnfrsf4", "Ldha", "Aldoa", "Il6st", "Bcl2l11", "Neurl3", "Txnip", "Foxo1", "Itga4", "Bcl11b", "Il4ra",
        "Sox4", "Ccl5", "Klf6", "Ctsw", "Cd96", "Satb1", "Lef1")

lr_features =  setdiff(choose_lr_features(id, t_cells, dc_cells, bad_genes, cor_n=150,
        must_haves = names(scdb_gset(id)@gene_set)), bad_genes)
mle_features = choose_mle_features(id, id, t_cells, dc_cells, union(bad_markers, bad_genes), nms_thresh=1, shared_thresh=1,
        existing_list= union(names(scdb_gset(id)@gene_set), must_haves))

mle_features = setdiff(union(mle_features, c("Cd4", "Cd8b1", "Cd8a")), bad_genes)

##############

t_nms = rev(read.table(paste0(supdir, "/FigS6b.txt"), stringsAsFactor=F)[[1]])
t_ord = intersect(clust_ord, t_clusts)
t_cl = scdb_mc(paste0(id, "_a"))
t_fp = log2(t_cl@mc_fp); colnames(t_fp) = names(table(sin_cl@mc[t_cells]))

IM = t_fp[t_nms, as.character(t_ord)]
IM = IM[ order(max.col(IM)),]
vct = factor(color2name[ sin_cl@colors[ as.numeric(colnames(IM))]], levels = lin_ord); names(vct) = colnames(IM)
pdf(paste0(supdir, "/FigS6b.pdf"), useDingbats=F, height=10, width=7)
par(mar = c(0.5,5,0.5,0.5), fig = c(0,1,0.1,1))
image.2(IM, balance=T, vct=vct, annotate="rows"); box()
par(fig = c(0,1,0,0.1), new=T)
image(matrix(seq_len(ncol(IM))), axes=F, col = sin_cl@colors[ as.numeric(colnames(IM))]); box()
dev.off()

###################


numis = 800; k = 5000
bad_cells = c()
nice_cells = setdiff(cells, bad_cells)
res = simulate_doublets(id, setdiff(t_cells, bad_cells), setdiff(dc_cells, bad_cells), k, numis = rep(numis, k))

sim_umis = res$sim_umis; sim_info = res$info
sim_cells = names(which(colSums(sim_umis) == numis))
sim_umis = sim_umis[,sim_cells]; sim_info = sim_info[sim_cells,]

sim_alpha = sim_info$alpha.1; names(sim_alpha) = rownames(sim_info)

sim_cells = sample(rownames(sim_info), 4000)
sim_mle_res = assign_pics_to_singlets(id, id, sim_umis[,sim_cells], setdiff(t_cells, bad_cells), setdiff(dc_cells, bad_cells), sim_alpha[sim_cells],
        verbose=T, bad_genes = bad_genes, markers = mle_features, reg = 1e-4)

sim_cells = rownames(sim_mle_res)

t_confu = table(sin_cl@mc[ as.vector( sim_info[sim_cells , "sim.1"])], sim_mle_res[sim_cells, "a_mc"])
t_n = t_confu / rowSums(t_confu)

grad = colorRampPalette(c("white", "#FDC51D", "#CA531C", "#951851", "#36277A", "black"))(1000)
t_cls = factor(color2name[ sin_cl@colors[ as.numeric(rownames(t_n))]], levels = lin_ord)
pdf(paste0(supdir, "FigS6d.pdf"), useDingbats=F, height=10, width=10)
par(mar = rep(0.5,4), fig = c(0.05,1,0.05,1))
image.2(t_n, zlim = c(0,1), col = grad, annotate = "none", hct = t_cls, vct = t_cls); box()
par(fig = c(0.05,1,0,0.05), new = T)
image(matrix(seq_along(t_cls)), axes = F, col = name2color[ as.vector(sort(t_cls))]); box()
par(fig = c(0,0.05,0.05,1), new = T)
image(t(seq_along(t_cls)), axes = F, col = name2color[ as.vector(sort(t_cls))]); box()
dev.off()

alpha_fit = estimate_mixing(sim_umis, sim_info$alpha.1, lr_features, paste0(supdir, "/alpha_estimation_all.png"))
alpha_tag = predict(alpha_fit, newx = t(sim_umis[lr_features,]), s = "lambda.min")[,1]

i = which(alpha_fit$lambda == alpha_fit$lambda.min)
pdf(paste0(supdir, "/FigS6c.pdf"), useDingbats=F)
plot(sim_info[, "alpha.1"], alpha_tag, pch = 20, col = rgb(0,0,0,0.6), xlim = c(0,1), ylim = c(0,1),
        cex=2, cex.main=2, main = round(1 - alpha_fit$cvm[i] / var(sim_info$alpha.1),4), axes=F)
abline(coef = c(0,1)); axis(1); axis(2,las=2)
dev.off()

###################

numis=800
db_mat = scdb_mat("mouse_PIC")
db_stats = db_mat@cell_metadata[ db_mat@cells,]
db_cells = rownames(db_stats)[db_stats$timepoint == "d10"]
db_umis = read_large_umis("mouse_PIC", cells = db_cells)
ds = .downsamp(db_umis, numis)
bad_cells = c()
pics = intersect(colnames(ds), db_mat@cells)
mle_res = run_pic_seq(id, id, ds[,pics], setdiff(t_cells, bad_cells), setdiff(dc_cells, bad_cells),
        lr_features, mle_features, paste0(supdir, "/alpha_estimation.png"),
        numis, reg = 1e-4) #, comb=sin_group)
mle_res$well = rownames(mle_res)

#################

good_pics = setdiff(rownames(mle_res)[ mle_res$alpha > 0 & mle_res$alpha < 1], c(names(which((ds["Igkc", ] > 1 & ds["Jchain", ] > 1)))))
#good_pics = intersect(good_pics, rownames(cell_stats)[ cell_stats$timepoint == "d10"])

alpha = mle_res[ good_pics, "alpha"]; names(alpha) = good_pics
t_mc = mle_res[good_pics, "a_mc"]; names(t_mc) = good_pics
dc_mc = mle_res[good_pics, "b_mc"]; names(dc_mc) = good_pics
#parser_t = t_anno[ as.character(t_mc[ good_pics])]; names(parser_t) = good_pics
parser_dc = color2name[ sin_cl@colors[ dc_mc[ good_pics]]]; names(parser_dc) = good_pics
parser_t = color2name[ sin_cl@colors[ t_mc[ good_pics]]]; names(parser_t) = good_pics

##############

cell_stats = rbind(sin_stats, db_stats)
gate = as.vector(cell_stats$sorting.scheme); names(gate) = rownames(cell_stats)
gate[ gate == "Cd11c+"] = "CD11c+"

bad_cells = names(sin_cl@mc)[ sin_cl@mc %in% bad_mc]
sin_cells = cells

comb = with(cell_stats, paste0(PIC, ":", gate, "@", tissue, "@", treatment, "@", timepoint, "@", as.numeric(factor(date)), "-", replicate)); names(comb) = rownames(cell_stats)
tcrb_cells = intersect(t_cells, rownames(cell_stats)[ cell_stats$sorting.scheme == "TCRb+"])

t_dist = rbind(table(comb[good_pics],  parser_t[good_pics]),
        table(comb[tcrb_cells],  sin_names[tcrb_cells]))
t_dist = t_dist[ rowSums(t_dist) >= 10,]
#t_dist = t_dist[ grep("Ova-Vaccination", rownames(t_dist)), intersect(lin_ord, colnames(t_dist))]
samp2rep = paste0(vecsplit(rownames(t_dist), "@",2), "@", vecsplit(rownames(t_dist), "@",5))
good_reps = names(which(table(samp2rep) == 2))
t_dist = t_dist[ samp2rep %in% good_reps, intersect(lin_ord, colnames(t_dist))]
t_n = t_dist / rowSums(t_dist)

rownames(t_n) = gsub(":[a-zA-Z0-9,\\+]*@", ":@", rownames(t_n))
samp2gate = vecsplit(rownames(t_n), "@", 1)
samp2gate = factor(samp2gate, levels = c("Singlets:", "PIC:"))

samp2site = vecsplit(rownames(t_n), "@", 2)
samp2site[ samp2site == "LN"] = "cLN"
samp2site = factor(samp2site, levels = c("cLN", "dLN", "tumor"))

samp2tp = vecsplit(rownames(t_n), "@", 4)
samp2treat = factor(vecsplit(rownames(t_n), "@", 3), levels = c("Ova-Vaccination", "No vaccination"))

samp2rep = vecsplit(rownames(t_n), "@", 5)

samp2comb = interaction(interaction(samp2gate, samp2site), samp2tp) #, interaction(samp2tp, samp2treat))
samp2comb = factor(samp2comb, levels = names(which(table(samp2comb) > 0)))
#samp2comb = factor(samp2comb, levels = c("Singlets:.cLN", "PIC:.cLN", "Singlets:.dLN", "PIC:.dLN"))

resampled_t = c(); resampled_pics = c()
patient = vecsplit(comb, "@", 5)
site = vecsplit(comb, "@", 2)
for (p in names(table(patient))) {
        for (s in names(table(site))) {
                sub_cells = intersect(union(t_cells, good_pics), names(which(patient == p & site == s)))
                if (length(sub_cells) == 0) { next}
                l = min(table(factor(vecsplit(comb[ sub_cells], ":",1), levels = c("Singlets", "PIC"))))
                resampled_t = c(resampled_t, sample(intersect(t_cells, sub_cells), l))
                resampled_pics = c(resampled_pics, sample(intersect(good_pics, sub_cells), l))
        }
}

t_all = rbind(table(paste0(site[resampled_pics], "@PIC"),  parser_t[resampled_pics]),
        table(paste0(site[ resampled_t], "@T"),  sin_names[resampled_t]))
t_all = t_all[c("cLN@T", "cLN@PIC", "dLN@T", "dLN@PIC"), intersect(lin_ord, colnames(t_all))]
t_all_n = t_all / rowSums(t_all)


t_pvals = matrix(NA, nrow=ncol(t_all), ncol=2, dimnames=list(colnames(t_all), c("cLN", "dLN")))
for (pop in colnames(t_n)) {
        X = cbind(t_all[,pop], rowSums(t_all[, setdiff(colnames(t_all), pop)]))
	t_pvals[pop,] = sapply(colnames(t_pvals), function(x) fisher.test(X[grep(x, rownames(X)),])$p.value)
}
t_qvals = matrix(p.adjust(t_pvals, "fdr"), nrow = nrow(t_pvals), dimnames = dimnames(t_pvals))

pdf(paste0(supdir, "/FigS6e.pdf"), useDingbats=F)
t_all_n = t_all / rowSums(t_all)
#t_all_n = t_all_n[ grep("d10", rownames(t_all_n)),]
t_all_n = t_all_n[ rev(rownames(t_all_n)),]
X = barplot(t(t_all_n), las=2, col = name2color[ colnames(t_all_n)], space = c(2,1), axes = F, names.arg = rep("", nrow(t_all_n)), horiz=T)
names(X) = rownames(t_all_n)
t_all_sum = t(apply(t_all_n,1,cumsum))
t_all_melt = melt(t_all_sum)
t_all_melt$y = X[as.vector(t_all_melt$Var1)]
left_melt = t_all_melt[grep("@PIC", t_all_melt$Var1),]
right_melt = t_all_melt[grep("@T", t_all_melt$Var1),]
both_melt = cbind(left_melt, right_melt[,c(3,4)])
colnames(both_melt) = c("name", "pop", "left_val", "left_y", "right_val", "right_y")
with(both_melt, segments(left_val, left_y + 0.5, right_val, right_y - 0.5, lty=2))
both_melt$desc = vecsplit(as.vector(both_melt$name), ":", 1)
both_melt$mid_y = with(both_melt, rowMeans(cbind(left_y, right_y)))
with(both_melt, axis(2, at = mid_y, labels = desc, las = 2, cex.axis=2))
axis(1)
dev.off()

disp_pops = rownames(t_qvals); n=length(disp_pops)
pdf(paste0(supdir, "/FigS6f.pdf"), useDingbats=F, height=7, width=5 * n)
par(mfrow = c(1,n), mar=c(5,10,5,1))
for (pop in disp_pops) {
    	x = as.numeric(samp2comb)# + runif(nrow(t_n), -0.2,0.2)
        idx = which(!is.na(x)); y = t_n[idx,pop];
        names(x) = rownames(t_n); names(y) = rownames(t_n)
        labels = unique(vecsplit(names(x), ":", 2))
        plot(1,1,xlim = c(1,length(table(samp2comb))), ylim = c(0,max(y)), axes=F, xlab="", ylab="", pch = 21, type = "n", main=pop)
        segments(x[paste0("Singlets:", labels)], y[paste0("Singlets:", labels)], x[paste0("PIC:", labels)], y[paste0("PIC:", labels)], lwd=1)
        points(x[idx], y, pch = 21, bg = name2color[pop], cex=3, main=pop, cex.main = 1)
        #text(x, t_n[,pop], samp2rep)
        #grid()
        pvals = rep(NA,2)
        pvals[1] = wilcox.test(y[paste0("Singlets:", grep("cLN", labels, v=T))], y[paste0("PIC:", grep("cLN", labels, v=T))], paired = TRUE, alternative = "two.sided")$p.value
        pvals[2] = wilcox.test(y[paste0("Singlets:", grep("dLN", labels, v=T))], y[paste0("PIC:", grep("dLN", labels, v=T))], paired = TRUE, alternative = "two.sided")$p.value
        axis(2); axis(1, at = seq_along(levels(samp2comb)), labels = levels(samp2comb), las=2, cex.axis=1)
        axis(3, at = c(1.5,3.5), labels = paste0("P=",round(pvals, 4)), cex.axis=1)
}
dev.off()

#############

good_pics = names(which(ds["Igkc", good_pics] < 2 | ds["Ighm", good_pics] < 2))
ds_f = ds[setdiff(rownames(sin_cl@e_gc), bad_genes), good_pics]
us =db_umis[ rownames(ds_f), colnames(ds_f)]

exp_us = generate_expected_pics_from_mle(id, mle_res[good_pics, c("a_mc", "b_mc")], mle_res[good_pics, "alpha"],
        colSums(us[,good_pics]), bad_genes = bad_genes)
exp_n = generate_expected_pics_from_mle(id, mle_res[good_pics, c("a_mc", "b_mc")], mle_res[good_pics, "alpha"],
        colSums(ds_f[,good_pics]), bad_genes = bad_genes)
t_n = generate_expected_pics_from_mle(id, mle_res[good_pics, c("a_mc", "b_mc")], rep(1, length(good_pics)),
        colSums(ds_f[,good_pics]) * alpha[good_pics], bad_genes = bad_genes)
dc_n = generate_expected_pics_from_mle(id, mle_res[good_pics, c("a_mc", "b_mc")], rep(0, length(good_pics)),
        colSums(ds_f[,good_pics]) * (1 - alpha[good_pics]), bad_genes = bad_genes)
genes = rownames(exp_us)

#############

pic_comb = factor(parser_t[ good_pics], levels = lin_ord[1:6])
comb2 = interaction(paste0(vecsplit(comb[good_pics], "@", 2), "@", vecsplit(comb[good_pics], "@", 5)), pic_comb, sep = ":")
names(comb2) = good_pics

t_pops = lin_ord[1:6]
cells2 = names(parser_t)[ parser_t %in% t_pops]
cells2 = cells2[ comb2[cells2] %in% names(which(table(comb2[cells2]) >= 10))]
comb2 = factor(comb2, levels = names(which(table(comb2[cells2]) > 0)))

bad_genes = grep("^Rpl|^Rps|^Gm[0-9]|^AC[0-9]|^Snor|^His[0-9]|^Hsp[0-9]|^Trim30|^Ig[k,h]", rownames(ds_f), v=T)
good_genes = setdiff(rownames(ds_f), bad_genes)
db_int = comb2[cells2]
names(db_int) = cells2
real_m = t(apply(ds_f[good_genes,  cells2], 1, tapply, db_int[cells2], sum))
exp_m =  t(apply(exp_n[good_genes, cells2], 1, tapply, db_int[cells2], sum))
x2 = rowSums(((real_m - exp_m)^2 )/exp_m)
qvals = p.adjust(1 - pchisq(x2, df = ncol(real_m) - 1), "fdr")
z = log2(real_m / exp_m);
z[real_m==0] = NA

y = rowSums(us[names(qvals),cells2]);
genes = union(setdiff(names(which(qvals < 1e-3 & apply(us[ names(qvals), cells2], 1, function(x) sort(x,T)[3]) > 2 & y > 10)), c()),
	c(must_haves, "Btla", "Pdcd1", "Cd274", "Tigit", "Ctla4", "Icos", "Cd40lg", 
	"Icosl", "Il12b", "Ebi3", "Cd40", "Ccl22", "Ccl17", "Il21", "Cd86", "Cd28", "Cd80",
	"Top2a", "Npm1", "Ube2c", "S100a6", "S100a4", "Lgals1", "Maf", "Mif", "Pou2af1", "Srm", "Eomes",
	"S100a10", "Slamf6", "Irf4", "Batf", "Sh2d1a", "Tiam1", "Tnfsf11", "Cd44",
	"Hopx", "Ikzf2", "Il22", "Cxcr6", "Cxcl10", "Cxcl9", "Myc", "Hif1a", "Tox2", "Cd82", "Ly6a", "Vim", "Pkm", "Tmem176b",
	"Mtap"))

z_reg = log2((real_m + 50) / (exp_m + 50))
#z_reg = log2((real_m + 10) / (exp_m + 10))
IM = z_reg[intersect(rownames(z_reg), genes), ]#grep(paste0("^|\\.", pop, "\\.|$"), colnames(z))]
IM = IM[rowSums(!is.na(IM)) > 0,]

x_min = apply(IM,1,min); x_max = apply(IM,1,max)
x = ifelse(abs(x_min) > x_max, x_min, x_max)

genes2 = c("Xcl1", "Mif", "Npm1", "Eomes", "Tigit", "Tpi1")
db_int = comb2 #paste0(comb[cells2], "@", parser_dc[cells2]); names(db_int) = cells2
real_m = t(apply(ds_f[genes2,cells2], 1, tapply, db_int[cells2], mean));
exp_m =  t(apply(exp_n[genes2,cells2], 1, tapply, db_int[cells2], mean));
t_m =  t(apply(t_n[genes2,cells2], 1, tapply, db_int[cells2], mean));
dc_m =  t(apply(dc_n[genes2,cells2], 1, tapply, db_int[cells2], mean));

library(Hmisc)
m = t(apply(ds_f[genes2, cells2],1,tapply, db_int[cells2], sum))
n = tapply(colSums(ds_f[,cells2]),db_int[cells2],sum)
Y = sweep(t(apply(m,1,binconf,n)), 2, rep(tapply(colSums(ds_f[,cells2]), db_int[cells2], mean), 3), "*")

reg = 0.001
grad = colorRampPalette(c("limegreen", "gray40", "firebrick3"))(101)
good_samps = grep("LN.*ot2", colnames(real_m), v=T)
pdf(paste0(supdir, "/FigS6g.pdf"), height=8, width=12, useDingbats=F)
par(mfrow = c(2,3))
for (gene in genes2) {
	obs_vals = real_m[gene, good_samps];
        vals_melt = melt(obs_vals)
        vals_melt$exp_value = melt(exp_m[gene, good_samps])$value
        vals_melt$t_value = melt(t_m[gene, good_samps])$value
        vals_melt$dc_value = melt(dc_m[gene, good_samps])$value
        vals_melt$site = #interaction(vecsplit(rownames(vals_melt), "@", 1),
		vecsplit(rownames(vals_melt), ":", 2) #, sep = ".")
        #vals_melt$site = factor(vecsplit(rownames(vals_melt), ":", 2))
	vals_melt$site = factor(vals_melt$site, levels = names(which(table(vals_melt$site) > 0)))
        #vals_melt$site = factor(vals_melt$site, levels = intersect(levels(pic_comb), names(which(table(vals_melt$site) > 0))))
        vals_melt$patient = vecsplit(rownames(vals_melt), "\\.", 1)
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
        with(vals_melt, points(exp_x, exp_value, cex=3, pch = 21, bg = col))
        with(vals_melt, points(obs_x, value, cex = 3, pch = 21, bg = "gray"))
	segments(seq_along(exp_med) * 2 - 0.4, exp_med, seq_along(exp_med) * 2 + 0.4, lwd = 1, col = "blue")
	segments(seq_along(obs_med) * 2 + 0.6, obs_med, seq_along(obs_med) * 2 + 1.4, lwd = 1, col = "blue")
        axis(2, las=2); axis(1, at = seq_along(levels(vals_melt$site)) * 2 + 0.5, labels = names(table(vals_melt$site)), las=2)
}
dev.off()

pdf(paste0(supdir, "/t_dc_colorbar.pdf"), useDingbats=F)
image(matrix(seq_along(grad)), col=grad, axes=F); box()
dev.off()

#############
