message("Generating Fig. 2 and S2")

#######

outdir = paste0("figures/figure2")
dir.create(outdir)
supdir = paste0("figures/figureS2")
dir.create(supdir)

############

numis = 800; k = 5000
comb = with(cell_stats, paste0(organ, "@", site, "@", patient))
names(comb) = rownames(cell_stats)
bad_cells = c()
nice_cells = setdiff(cells, bad_cells)
res = simulate_doublets(id, t_cells, dc_cells, k, comb, numis = rep(numis, k))

sim_umis = res$sim_umis; sim_info = res$info
sim_cells = names(which(colSums(sim_umis) == numis))
sim_umis = sim_umis[,sim_cells]; sim_info = sim_info[sim_cells,]

sim_alpha = sim_info$alpha.1; names(sim_alpha) = rownames(sim_info)

sim_cells = sample(rownames(sim_info), 4000)
sim_mle_res = assign_pics_to_singlets(id, id, sim_umis[,sim_cells], t_cells, dc_cells, sim_alpha[sim_cells],
        verbose=T, bad_genes = bad_genes, markers = mle_features, reg = 1e-4)

sim_cells = rownames(sim_mle_res)
dc_confu = table(sin_cl@mc[ as.vector( sim_info[ sim_cells, "sim.2"])], sim_mle_res$b_mc)
dc_n = dc_confu / rowSums(dc_confu)

grad = colorRampPalette(c("white", "#FDC51D", "#CA531C", "#951851", "#36277A", "black"))(1000)
dc_cls = factor(color2name[ sin_cl@colors[ as.numeric(rownames(dc_n))]], levels = lin_ord)
pdf(paste0(supdir, "/FigS2b_right.pdf"), useDingbats=F, height=8, width=8)
par(mar = rep(0.5,4), lwd = 3, fig = c(0.1,1,0.1,1))
image.2(dc_n, zlim = c(0,1), col = grad, annotate = "none", hct = dc_cls,vct = dc_cls); box()
par(fig = c(0.1,1,0,0.1), new = T)
image(matrix(seq_along(dc_cls)), axes = F, col = name2color[ as.vector(sort(dc_cls))]);box()
par(fig = c(0,0.1,0.1,1), new = T)
image(t(seq_along(dc_cls)), axes = F, col = name2color[ as.vector(sort(dc_cls))]); box()
dev.off()

t_confu = table(sin_cl@mc[ as.vector( sim_info[sim_cells , "sim.1"])], sim_mle_res[sim_cells, "a_mc"])
t_n = t_confu / rowSums(t_confu)

t_cls = factor(color2name[ sin_cl@colors[ as.numeric(rownames(t_n))]], levels = lin_ord)
pdf(paste0(supdir, "/FigS2b_left.pdf"), useDingbats=F, height=8, width=8)
par(mar = rep(0.5,4), lwd = 3, fig = c(0.1,1,0.1,1))
image.2(t_n, zlim = c(0,1), col = grad, annotate = "none", hct = t_cls, vct = t_cls); box()
par(fig = c(0.1,1,0,0.1), new = T)
image(matrix(seq_along(t_cls)), axes = F, col = name2color[ as.vector(sort(t_cls))]); box()
par(fig = c(0,0.1,0.1,1), new = T)
image(t(seq_along(t_cls)), axes = F, col = name2color[ as.vector(sort(t_cls))]); box()
dev.off()

pdf(paste0(supdir, "/mle_confu_cb.pdf"), useDingbats=F)
par(mar = rep(0,4))
image(matrix(seq_along(grad)), axes=F, col = grad); box()
dev.off()

##############

nk_sim = intersect(sim_cells, rownames(res$info)[ sin_names[as.vector(res$info$sim.2)] %in% c("NK", "Cyto_NK")])
mye_sim = setdiff(sim_cells, nk_sim)
alpha_fit = estimate_mixing(sim_umis, sim_info$alpha.1, lr_features, paste0(supdir, "/alpha_estimation_all.png"))
alpha_tag = predict(alpha_fit, newx = t(sim_umis[lr_features,]), s = "lambda.min")[,1]

pdf(paste0(supdir, "/FigS2a_left.pdf"), useDingbats=F)
plot(sim_info[mye_sim, "alpha.1"], alpha_tag[mye_sim], pch = 20, col = rgb(0,0,0,0.6), xlim = c(0,1), ylim = c(0,1),
	cex=1, cex.main=2, main = round(cor(sim_info[mye_sim, "alpha.1"], alpha_tag[mye_sim])^2,3), axes=F, xlab="", ylab="")
axis(1); axis(2)
abline(coef = c(0,1))
dev.off()

pdf(paste0(supdir, "/FigS2a_right.pdf"), useDingbats=F)
plot(sim_info[nk_sim, "alpha.1"], alpha_tag[nk_sim], pch = 20, col = rgb(0,0,0,0.6), xlim = c(0,1), ylim = c(0,1),
	cex=1, cex.main=2, main = round(cor(sim_info[nk_sim, "alpha.1"], alpha_tag[nk_sim])^2,3), axes=F, xlab="", ylab="")
abline(coef = c(0,1))
axis(1); axis(2)
dev.off()

##############

id_d = "human_PIC"
db_mat = scdb_mat(id_d)
db_umis = read_large_umis(id_d)
db_cells = db_mat@cells
umis = cbind(sin_umis, db_umis[ rownames(sin_umis),])
numis=800
ds = .downsamp(umis, numis)
comb = with(cell_stats, paste0(organ, "@", site, "@", patient)); names(comb) = rownames(cell_stats)
X = table(comb[ cells], cells %in% t_cells)
good_combs = names(which(rowSums(X < 20) == 0))
good_cells = cells[ comb[cells] %in% good_combs]
mle_res = run_pic_seq(id, id, ds, intersect(good_cells, t_cells), intersect(good_cells, dc_cells), lr_features, mle_features, paste0(outdir, "/alpha_estimation.png"),
	numis, comb = comb, reg = 1e-4)
mle_res$well = rownames(mle_res)

cell_stats = scdb_mat("all_human")@cell_metadata[ union(cells, db_cells),]
mle_res$type = as.vector(cell_stats[rownames(mle_res), "PIC"])
mle_res[ rownames(mle_res) %in% t_cells, "type"] = "T"
mle_res[ rownames(mle_res) %in% dc_cells, "type"] = "DC"
mle_res$sin_alpha = with(mle_res, ifelse(type == "PIC", 1 * (alpha >= 0.5), 1 * (well %in% t_cells)))
mle_res$sin_ll = with(mle_res, pic_ll_to_pair(id, id, ds[,well], a_mc, b_mc, sin_alpha, reg = 1e-4, markers = mle_features))
mle_res$diff = with(mle_res, ll - sin_ll)

##############

gate = as.vector(cell_stats$sorting.scheme); names(gate) = rownames(cell_stats)
gate = c("T", "PIC", "APC")[ as.numeric(factor(gate))]; names(gate) = rownames(cell_stats)

nk_pic = rownames(mle_res)[ color2name[ sin_cl@colors[ mle_res$b_mc]] %in% c("NK", "NK_CX3CR1")]
good_pics = setdiff(rownames(mle_res)[ mle_res$type == "PIC" & mle_res$diff > 0 & mle_res$alpha > 0 & mle_res$alpha < 1], c(NA, nk_pic))

alpha = mle_res[ good_pics, "alpha"]; names(alpha) = good_pics
t_mc = mle_res[good_pics, "a_mc"]; names(t_mc) = good_pics
dc_mc = mle_res[good_pics, "b_mc"]; names(dc_mc) = good_pics
parser_t = color2name[ sin_cl@colors[ t_mc[ good_pics]]]; names(parser_t) = good_pics
parser_dc = color2name[ sin_cl@colors[ dc_mc[ good_pics]]]; names(parser_dc) = good_pics

sin_cells = setdiff(intersect(cells, names(gate)[ gate %in% c("APC", "T")]), 
	c(names(sin_names)[ sin_names %in% c("NK", "NK_CX3CR1", "Baso")], bad_cells))

##############

t_col = "limegreen"; dc_col =  "firebrick3"; db_col = "orange2"
cols = c(t_col, dc_col, db_col); names(cols) = c("T", "DC", "PIC")
pdf(paste0(supdir, "/FigS2c.pdf"), useDingbats=F)
plot(1,1,xlim = quantile(mle_res$diff, c(0,1), na.rm=T), ylim = c(0,1), type="n", axes=F, xlab="", ylab="")
abline(v=0, lwd=3); grid(col="black")
axis(1); axis(2)
with(mle_res, sapply(names(cols), function(x) lines(ecdf(diff[type == x]), col = cols[x], lwd=4)))
dev.off()

all_pics = db_mat@cells
patient_dist = table(cell_stats[ all_pics, "patient"], all_pics %in% good_pics)
patient_dist = patient_dist[ rowSums(patient_dist) > 0,]
rownames(patient_dist) = seq_len(nrow(patient_dist))
dist_n = patient_dist / rowSums(patient_dist)
ylim = c(0, max(dist_n[,2]) + 0.05)
pdf(paste0(supdir, "/FigS2d.pdf"), useDingbats=F)
X = barplot(dist_n[,2], ylim = ylim)
text(X, dist_n[,2] + 0.02, patient_dist[,2])
dev.off()

##############

cells = union(t_cells, dc_cells)
t_nms = rev(read.table(paste0(supdir, "/FigS2f-h_bottom.txt"), stringsAsFactor=F)[[1]])
dc_nms = rev(read.table(paste0(supdir, "/FigS2f-h_top.txt"), stringsAsFactor=F)[[1]])

a_lfp = log2(scdb_mc(paste0(id, "_a"))@mc_fp)
t_nms = t_nms[ order(max.col(a_lfp[t_nms, as.numeric(factor(intersect(clust_ord, t_clusts)))]))]

b_lfp = log2(scdb_mc(paste0(id, "_b"))@mc_fp)
dc_nms = dc_nms[ order(max.col(b_lfp[dc_nms, as.numeric(factor(intersect(clust_ord, dc_clusts)))]))]

nms = union(t_nms, dc_nms)
bad_cells = names(sin_cl@mc)[ sin_cl@mc %in% bad_clusts]

all_pics = good_pics; 

pdf(paste0(supdir, "/FigS2f-j.pdf"), height=20, width=20, useDingbats=F)
sin_vec = color2name[ sin_cl@colors[sin_cl@mc[sin_cells]]]; names(sin_vec) = sin_cells
par(mar = c(0.5,5,0.5,0.5), fig = c(0,0.5,0.2,0.575), lwd = 1)
sin_ord = plot_sc_heatmap(id, id, t_nms, clusts = sin_vec, good_clusts = intersect(lin_ord, names(table(sin_vec))), cells = sin_cells, 
	annotate=T, annotate_what="rows", normalize=T, lty=1, draw_cls=F); box(lwd=1)
par(fig = c(0,0.5,0.575,0.95), new = T)
sin_ord = plot_sc_heatmap(id, id, dc_nms, clusts = sin_vec, good_clusts = intersect(lin_ord, names(table(sin_vec))), cells = sin_cells, 
	annotate=T, annotate_what="rows", normalize=T, lty=1, draw_cls=F); box(lwd=1)

par(fig=c(0,0.5,0.95,1), new=T)
image(matrix(seq_along(sin_ord)), col = ifelse(cell_stats[sin_ord, "site"] == "tumor", "gray20", "gray80"), axes=F); box(lwd=1)

par(fig=c(0,0.5,0.1,0.15), new=T)
image(matrix(ifelse(sin_ord %in% t_cells, as.numeric(factor(sin_vec[sin_ord])), NA)),  axes = F, col = name2color[names(table(sin_vec))], zlim = c(1, length(table(sin_vec)))); box(lwd=1)
par(fig=c(0,0.5,0.15,0.2), new=T)
image(matrix(ifelse(sin_ord %in% dc_cells, as.numeric(factor(sin_vec[sin_ord])), NA)), axes = F, col = name2color[names(table(sin_vec))], zlim = c(1, length(table(sin_vec)))); box(lwd=1)

db_clusts = interaction(factor(parser_dc[good_pics], levels = lin_ord), factor(parser_t[good_pics], levels = lin_ord))
db_vec = as.vector(db_clusts); names(db_vec) = good_pics
par(fig=c(0.5,1,0.2,0.575), new=T)
db_ord = plot_sc_heatmap(id_d, id_d, t_nms, clusts = db_vec, cells = good_pics, good_clusts = names(table(db_clusts)),
        annotate=T, annotate_what="rows", normalize=T, lty=1, draw_cls=F); box(lwd=1)
par(fig = c(0.5,1,0.575,0.95), new = T)
db_ord = plot_sc_heatmap(id_d, id_d, dc_nms, clusts = db_vec, cells = good_pics, good_clusts = names(table(db_clusts)),
        annotate=T, annotate_what="rows", normalize=T, lty=1, draw_cls=F); box(lwd=1)
par(fig=c(0.5,1,0.95,1), new=T)
image(matrix(seq_along(db_ord)), col = ifelse(cell_stats[db_ord, "site"] == "tumor", "gray20", "gray80"), axes=F); box(lwd=1)

par(fig=c(0.5,1,0.1,0.15), new=T)
image(matrix(seq_along(db_ord)), axes = F, col = sin_cl@colors[ as.numeric(t_mc[db_ord])]); box(lwd=1)
par(fig=c(0.5,1,0.15,0.2), new=T)
image(matrix(dc_mc[ db_ord]), axes = F, col = sin_cl@colors, zlim = c(1, max(sin_cl@mc))); box(lwd=1)

t_col = "limegreen"; dc_col =  "firebrick3"
par(fig=c(0.5,1,0,0.1), new=T)
split_count = cbind(alpha[good_pics], 1 - alpha[good_pics]) * umicount[good_pics]
split_count = split_count / rowSums(split_count)
barplot(t(split_count[db_ord,]), col = c(t_col, dc_col), border = NA, xaxs = "i", yaxs="i", space = 0, names.arg = rep("", length(db_ord)), las = 2)
box(lwd=1)
dev.off()

##############

comb = with(cell_stats, paste0(site, "@", patient, "@", organ, ":", gate)); names(comb) = rownames(cell_stats)
t_dist = rbind(table(comb[good_pics],  parser_t[good_pics]),
	table(comb[intersect(t_cells, sin_cells)],  sin_names[intersect(t_cells, sin_cells)]))
t_dist = t_dist[ rowSums(t_dist) >= 5,]
t_dist = t_dist[ -grep(":APC", rownames(t_dist)),]
t_dist = t_dist[ order(vecsplit(rownames(t_dist), ":", 1), factor(vecsplit(rownames(t_dist), ":", 2), levels = c("T", "APC", "PIC"))),]
desc2rep = vecsplit(rownames(t_dist), ":", 1)
t_dist = t_dist[ sapply(names(which(table(desc2rep) == 2)), grep, rownames(t_dist), v=T), intersect(lin_ord, colnames(t_dist))]
t_n = t_dist / rowSums(t_dist)
desc2rep = desc2rep[ rownames(t_dist)]

desc2all = paste0(vecsplit(rownames(t_dist), "@", 1), "@", vecsplit(rownames(t_dist), ":", 2))
t_all = apply(t_dist, 2, tapply, factor(desc2all, levels = c("normal@T", "normal@PIC", "tumor@T", "tumor@PIC")), sum)
t_all_n = t_all / rowSums(t_all)

t_sum = t(apply(t_n,1,cumsum))
t_melt = melt(t_sum)
t_melt$x = X[as.vector(t_melt$Var1)]
t_melt$freq = melt(t_n)[,3]
right_melt = t_melt[grep(":PIC", t_melt$Var1),]
left_melt = t_melt[grep(":T", t_melt$Var1),]
both_melt = cbind(left_melt, right_melt[,c(3,4)])
colnames(both_melt) = c("name", "pop", "left_val", "left_x", "right_val", "right_x")
both_melt$desc = vecsplit(as.vector(both_melt$name), ":", 1)
both_melt$right_freq = t_melt[grep(":PIC", t_melt$Var1),"freq"]
both_melt$left_freq = t_melt[grep(":T", t_melt$Var1),"freq"]
both_melt$patient = vecsplit(both_melt$desc, "@", 2)
both_melt$site = vecsplit(both_melt$desc, "@", 1)
both_melt$pop = as.vector(both_melt$pop)

resampled_t = c(); resampled_pics = c()
patient = vecsplit(comb, "@", 2)
site = vecsplit(comb, "@", 1)
for (p in names(table(patient))) {
	for (s in names(table(site))) {
		sub_cells = intersect(union(t_cells, good_pics), names(which(patient == p & site == s)))
		if (length(sub_cells) == 0) { next}
		l = min(table(factor(vecsplit(comb[ sub_cells], ":",2), levels = c("T", "PIC"))))
		resampled_t = c(resampled_t, sample(intersect(t_cells, sub_cells), l))
		resampled_pics = c(resampled_pics, sample(intersect(good_pics, sub_cells), l))
	}
}

t_all = rbind(table(paste0(site[resampled_pics], "@PIC"),  parser_t[resampled_pics]),
	table(paste0(site[ resampled_t], "@T"),  sin_names[resampled_t]))
t_all = t_all[c("normal@T", "normal@PIC", "tumor@T", "tumor@PIC"), intersect(lin_ord, colnames(t_all))]
t_all_n = t_all / rowSums(t_all)

disp_pops = c("T_CD8_Exhausted", "Tht")
n=length(disp_pops)
pdf(paste0(outdir, "/Fig2a.pdf"), useDingbats=F, height=7, width=5 * n)
par(mfrow = c(1,n), mar=c(5,10,5,1))
for (pop in disp_pops) {
        pop_melt = both_melt[ both_melt$pop == pop,]
        pop_melt$left_coord = ifelse(pop_melt$site == "tumor", 2, 0)
        pop_melt$right_coord = 1 + pop_melt$left_coord
        with(pop_melt, plot(c(left_coord, right_coord), c(left_freq, right_freq), type="n", axes=F, xlab="", ylab="", main = unique(anno_names[pop])))
        with(pop_melt, segments(left_coord, left_freq, right_coord, right_freq, lwd=1))
        with(pop_melt, points(c(left_coord, right_coord), c(left_freq, right_freq),
                pch= 21, cex=3, bg=name2color[pop]))
        pvals = rep(NA, 2); names(pvals) = names(table(pop_melt$site))
        for (s in names(table(pop_melt$site))) { pvals[s] =  with(pop_melt[ pop_melt$site == s,], wilcox.test(left_freq, right_freq, paired=T)$p.value)}
        axis(3, at = c(0.5,2.5), labels = paste0(c("Normal tissue\n", "TME\n"), "(P=",round(pvals, 4), ")"), cex.axis=1)
        axis(2, las=2)
	axis(1,at=0:3,labels=rep(c(expression(paste("TCR", beta^"+")), "PIC"),2))
}
dev.off()

pdf(paste0(outdir, "/Fig2c.pdf"), useDingbats=F)
par(lwd=1)
t_all_n = t_all / rowSums(t_all)
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
#with(both_melt, axis(2, at = mid_y, labels = desc, las = 2, cex.axis=2))
axis(1)
dev.off()

#############

dc_dist = rbind(table(comb[good_pics],  color2name[ sin_cl@colors[dc_mc[good_pics]]]),
	table(comb[intersect(dc_cells, sin_cells)],  color2name[ sin_cl@colors[sin_cl@mc[intersect(dc_cells, sin_cells)]]]))
dc_dist = dc_dist[ rowSums(dc_dist) > 6,]
dc_dist = dc_dist[ -grep(":T", rownames(dc_dist)),]
dc_dist = dc_dist[ order(vecsplit(rownames(dc_dist), ":", 1), factor(vecsplit(rownames(dc_dist), ":", 2), levels = c("T", "APC", "PIC"))),]
desc2rep = vecsplit(rownames(dc_dist), ":", 1)
dc_dist = dc_dist[ sapply(names(which(table(desc2rep) == 2)), grep, rownames(dc_dist), v=T), intersect(lin_ord, colnames(dc_dist))]
#dc_dist = dc_dist[seq_len(nrow(dc_dist)) + c(1,-1),]
dc_n = dc_dist / rowSums(dc_dist)
desc2rep = desc2rep[ rownames(dc_dist)]

desc2all = paste0(vecsplit(rownames(dc_dist), "@", 1), "@", vecsplit(rownames(dc_dist), ":", 2))
dc_all = apply(dc_dist, 2, tapply, factor(desc2all, levels = c("normal@APC", "normal@PIC", "tumor@APC", "tumor@PIC")), sum)
dc_all_n = dc_all / rowSums(dc_all)

dc_sum = t(apply(dc_n,1,cumsum))
dc_melt = melt(dc_sum)
dc_melt$x = X[as.vector(dc_melt$Var1)]
right_melt = dc_melt[grep(":PIC", dc_melt$Var1),]
left_melt = dc_melt[grep(":APC", dc_melt$Var1),]
both_melt = cbind(left_melt, right_melt[,c(3,4)])
colnames(both_melt) = c("name", "pop", "left_val", "left_x", "right_val", "right_x")
both_melt$desc = vecsplit(as.vector(both_melt$name), ":", 1)
dc_melt$freq = melt(dc_n)[,3]
both_melt$right_freq = dc_melt[grep(":PIC", dc_melt$Var1),"freq"]
both_melt$left_freq = dc_melt[grep(":APC", dc_melt$Var1),"freq"]
both_melt$patient = vecsplit(both_melt$desc, "@", 2)
both_melt$site = vecsplit(both_melt$desc, "@", 1)
both_melt$pop = as.vector(both_melt$pop)

resampled_dc = c(); resampled_pics = c()
patient = vecsplit(comb, "@", 2)
site = vecsplit(comb, "@", 1)
for (p in names(table(patient))) {
        for (s in names(table(site))) {
                sub_cells = intersect(union(intersect(sin_cells, dc_cells), good_pics), names(which(patient == p & site == s)))
                if (length(sub_cells) == 0) { next}
                l = min(table(factor(vecsplit(comb[ sub_cells], ":",2), levels = c("APC", "PIC"))))
                resampled_dc = c(resampled_dc, sample(intersect(intersect(sin_cells, dc_cells), sub_cells), l))
                resampled_pics = c(resampled_pics, sample(intersect(good_pics, sub_cells), l))
        }
}

dc_all = rbind(table(paste0(site[resampled_pics], "@PIC"),  parser_dc[resampled_pics]),
        table(paste0(site[ resampled_dc], "@APC"),  sin_names[resampled_dc]))
dc_all = dc_all[c("normal@APC", "normal@PIC", "tumor@APC", "tumor@PIC"), intersect(lin_ord, colnames(dc_all))]
dc_all_n = dc_all / rowSums(dc_all)

dc_pvals = matrix(NA, nrow=ncol(dc_all), ncol=2, dimnames=list(colnames(dc_all), c("normal", "tumor")))
for (pop in colnames(dc_n)) {
        X = cbind(dc_all[,pop], rowSums(dc_all[, setdiff(colnames(dc_all), pop)]))
        dc_pvals[pop, "normal"] = fisher.test(X[grep("normal", rownames(X)),])$p.value
        dc_pvals[pop, "tumor"]  = fisher.test(X[grep("tumor", rownames(X)),])$p.value
}

disp_pops = c("Monocytes_CD31", "MonMac", "TAM_TREM2", "cDC2", "cDC1", "MigDC")
n=length(disp_pops)
pdf(paste0(outdir, "/Fig2d.pdf"), useDingbats=F, height=7, width=5 * n)
par(mfrow = c(1,n), mar=c(5,10,5,1))
for (pop in disp_pops) {
        pop_melt = both_melt[ both_melt$pop == pop,]
        pop_melt$left_coord = ifelse(pop_melt$site == "tumor", 2, 0)
        pop_melt$right_coord = 1 + pop_melt$left_coord
        with(pop_melt, plot(c(left_coord, right_coord), c(left_freq, right_freq), type="n", axes=F, xlab="", ylab="", main = unique(anno_names[pop])))
        with(pop_melt, segments(left_coord, left_freq, right_coord, right_freq, lwd=1))
        with(pop_melt, points(c(left_coord, right_coord), c(left_freq, right_freq),
                pch= 21, cex=3, bg=name2color[pop]))
	pvals = rep(NA, 2); names(pvals) = names(table(pop_melt$site))
	for (s in names(table(pop_melt$site))) { pvals[s] =  with(pop_melt[ pop_melt$site == s,], wilcox.test(left_freq, right_freq, paired=T)$p.value)}
        #pvals = with(pop_melt, tapply(site, site, function(x) wilcox.test(left_freq[site == x], right_freq[ site == x], paired = TRUE, alternative = "two.sided")$p.value))
	axis(3, at = c(0.5,2.5), labels = paste0(c("Normal tissue\n", "TME\n"), "(P=",round(pvals, 4), ")"), cex.axis=1)
        axis(2, las=2)
	axis(1,at=0:3,labels=rep(c("Myeloid", "PIC"),2))
}
dev.off()

pdf(paste0(outdir, "/Fig2f.pdf"), useDingbats=F)
par(lwd=1)
dc_all_n = dc_all / rowSums(dc_all)
dc_all_n = dc_all_n[ rev(rownames(dc_all_n)),]
X = barplot(t(dc_all_n), las=2, col = name2color[ colnames(dc_all_n)], space = c(2,1), axes = F, names.arg = rep("", nrow(dc_all_n)), horiz=T)
names(X) = rownames(dc_all_n)
dc_all_sum = t(apply(dc_all_n,1,cumsum))
dc_all_melt = melt(dc_all_sum)
dc_all_melt$y = X[as.vector(dc_all_melt$Var1)]
left_melt = dc_all_melt[grep("@PIC", dc_all_melt$Var1),]
right_melt = dc_all_melt[grep("@APC", dc_all_melt$Var1),]
both_melt = cbind(left_melt, right_melt[,c(3,4)])
colnames(both_melt) = c("name", "pop", "left_val", "left_y", "right_val", "right_y")
with(both_melt, segments(left_val, left_y + 0.5, right_val, right_y - 0.5, lty=2))
both_melt$desc = vecsplit(as.vector(both_melt$name), ":", 1)
both_melt$mid_y = with(both_melt, rowMeans(cbind(left_y, right_y)))
#with(both_melt, axis(2, at = mid_y, labels = desc, las = 2, cex.axis=2))
axis(1)
dev.off()

###############

cell_ord = sample(good_pics[ order(factor(parser_t, levels = lin_ord))])
disp_genes = rev(read.table( paste0(outdir, "/Fig2b.txt"), stringsAsFactor=F)[[1]])

IM = log(1 + 7 * sweep(ds[disp_genes,cell_ord],2,alpha[cell_ord],"*"))
vct = factor(parser_t, levels = lin_ord[1:11]); names(vct) = good_pics
IM2 = log(1 + 7 * ds[disp_genes, sample(intersect(colnames(ds), t_cells))])
vct2 = factor(sin_names[t_cells], levels = lin_ord[1:11]); names(vct) = good_pics
grad = colorRampPalette(c("white", "orange", "tomato", "mediumorchid4", "midnightblue"))(1000)
pdf(paste0(outdir, "/Fig2b.pdf"), useDingbats=F, height=5, width=10)
par(fig = c(0,0.5,0.1,1), mar = c(0.5,5,0.5,0.5))
image.2(IM2, col = grad, vct = vct2[colnames(IM2)], annotate="rows", zlim = c(0, max(c(IM, IM2)))); box(lwd=1)
par(fig = c(0,0.5,0,0.1), new=T)
image(matrix(seq_len(ncol(IM2))), axes=F, col = name2color[ sin_names[ colnames(IM2)[ order(vct2[ colnames(IM2)])]]]); box(lwd=1)
par(fig = c(0.5,1,0.1,1), new=T) #, mar = rep(0.5,4))
image.2(IM, col = grad, vct = vct[colnames(IM)], annotate="none", zlim = c(0, max(c(IM, IM2)))); box(lwd=1)
par(fig = c(0.5,1,0,0.1), new=T)
image(matrix(seq_len(ncol(IM))), axes=F, col = name2color[ parser_t[ colnames(IM)[ order(vct[ colnames(IM)])]]]); box(lwd=1)
dev.off()

pdf(paste0(supdir, "/gene_colorbar.pdf"), useDingbats=F)
par(mar = rep(0.5,4))
image(matrix(seq_len(1000)), col = grad, axes=F); box()
dev.off()

cell_ord = sample(good_pics[ order(factor(parser_dc, levels = lin_ord))])
disp_genes = rev(read.table( paste0(outdir, "/Fig2e.txt"), stringsAsFactor=F)[[1]])
IM = log(1 + 7 * sweep(ds[disp_genes,cell_ord],2,alpha[cell_ord],"*"))
vct = factor(parser_dc, levels = lin_ord[14:23]); names(vct) = good_pics

IM2 = log(1 + 7 * ds[disp_genes, sample(intersect(colnames(ds), setdiff(dc_cells, names(sin_names)[ sin_names %in% c("NK", "Cyto_NK")])))])
vct2 = factor(sin_names[dc_cells], levels = lin_ord[14:23]); names(vct2) = dc_cells

grad = colorRampPalette(c("white", "orange", "tomato", "mediumorchid4", "midnightblue"))(1000)
pdf(paste0(outdir, "/Fig2e.pdf"), useDingbats=F, height=5, width=10)
par(fig = c(0,0.5,0.1,1), mar = c(0.5,5,0.5,0.5))
image.2(IM2, col = grad, vct = vct2[colnames(IM2)], annotate="rows", zlim = c(0, max(c(IM, IM2)))); box(lwd=1)
par(fig = c(0,0.5,0,0.1), new=T)
image(matrix(seq_len(ncol(IM2))), axes=F, col = name2color[ sin_names[ colnames(IM2)[ order(vct2[ colnames(IM2)])]]]); box(lwd=1)
par(fig = c(0.5,1,0.1,1), new=T) #, mar = rep(0.5,4))
image.2(IM, col = grad, vct = vct[colnames(IM)], annotate="none", zlim = c(0, max(c(IM, IM2)))); box(lwd=1)
par(fig = c(0.5,1,0,0.1), new=T)
image(matrix(seq_len(ncol(IM))), axes=F, col = name2color[ parser_dc[ colnames(IM)[ order(vct[ colnames(IM)])]]]); box(lwd=1)
dev.off()

############

mix_data = as.matrix(read.delim("Source_data/Source_data_mixing.txt", stringsAsFactors = F, row.names=1))
mix_t = c("FITC", "PE-Cy7", "FITC", "PE-Cy7")
mix_dc = c("PE", "PE", "APC-Cy7", "APC-Cy7")
colnames(mix_data) = paste0(mix_t, ":", mix_dc)

t_mat = matrix(F, nrow=2, ncol=4, dimnames = list(c("FITC", "PE-Cy7"), 1:4))
t_mat[ cbind(mix_t, 1:4)] = T

dc_mat = matrix(F, nrow=2, ncol=4, dimnames = list(c("PE", "APC-Cy7"), 1:4))
dc_mat[ cbind(mix_dc, 1:4)] = T

cols = c("white", "#728EC7","#2C3F86")

pdf(paste0(supdir, "/FigS2e.pdf"), useDingbats = F, height=10, width=7)
par(fig = c(0,1,0.4,1), mar = c(1,5,1,1))
barplot(mix_data, beside=T, xaxs="i", las=2, col = c("gray80", "gray20"))
par(fig = c(0,1,0.2,0.4), new=T)
image.2(t_mat[2:1,] * c(1,2), col = cols, annotate="rows"); box()
grid(nx=4, ny=2, col="black", lty=1)
par(fig = c(0,1,0,0.2), new=T)
image.2(dc_mat[2:1,] * c(1,2), col = cols, annotate="rows"); box()
grid(nx=4, ny=2, col="black", lty=1)
dev.off()
