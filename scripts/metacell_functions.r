


vecsplit = function(strvec, del, i) {
	 unlist(lapply(sapply(strvec, strsplit, del), "[[", i))
}

image.2 = function(X, col = colorRampPalette(c("blue", "white", "red"))(1000), balance=F, annotate="both", zlim = NULL,
        hct = NULL, vct = NULL, lwd=1, lty=1, cex = 1, text=F, text_mat = X, las=2) {
        if (is.null(zlim)) {
                if (balance) {
                                zlim = max(abs(X), na.rm=T) * c(-1,1)
                } else {
                        zlim = quantile(X, c(0,1), na.rm=T)
                }
	}
        hcls = NULL; vcls = NULL; hticks = NULL; vticks = NULL
        nrow = nrow(X); ncol = ncol(X)
        if (!is.null(hct)) {
                X = X[order(hct),]; text_mat = text_mat[order(hct),]
                hticks = seq(-1 / (2 * (nrow-1)),1 + 1 / (2 * (nrow-1)), length.out = nrow + 1)
		hcls = cumsum(table(hct)) + 1
        }
        if (!is.null(vct)) {
                X = X[,order(vct)]; text_mat = text_mat[,order(vct)]
                vticks = seq(-1 / (2 * (ncol-1)),1 + 1 / (2 * (ncol-1)), length.out = ncol + 1)
		vcls = cumsum(table(vct)) + 1
        }
        message("zlim: ", zlim[1], "<>", zlim[2])
        image(t(X), axes = F, col = col, zlim = zlim)
        abline(h = hticks[hcls], v = vticks[vcls], lwd=lwd, lty=lty)
        if (annotate %in% c("rows", "both")) {
	        mtext(rownames(X), las = las, side=2, at = (1 - seq_len(nrow(X))) / (1 - nrow(X)), cex = cex)
                if (!is.null(hct)) {
                        mtext(names(hcls), side = 4, las = las, at = rowMeans(cbind(hticks[c(1,hcls[-length(hcls)])], hticks[hcls])), cex = cex)}
	}
	if (annotate %in% c("columns", "both")) {
	        mtext(colnames(X), las = las, side=1, at = (1 - seq_len(ncol(X))) / (1 - ncol(X)), cex = cex)
                if (!is.null(vct)) {
                                   mtext(names(vcls), side = 3, las = las, at = rowMeans(cbind(vticks[c(1,vcls[-length(vcls)])], vticks[vcls])), cex = cex)}
        }
	if (text) {
           hmed = seq(0,1,length.out=nrow); vmed = seq(0,1,length.out=ncol)
           text(rep(vmed, each = nrow), rep(hmed, ncol), text_mat)
	}
}

bi_correlation_plot = function(X, m="spearman", fname=NULL,
	row_cols = NULL, col_cols = NULL, row_mar = 5, col_mar = 5,
	min_height=3000, min_width=3000,
	cor_palette = colorRampPalette(c("blue", "white", "red"))(1000),
	X_palette = colorRampPalette(c("purple", "gray20", "green"))(1000), balance_X = T) {

	row_C = cor(t(X), m=m); 
	diag(row_C) = NA
	row_hc = hclust(as.dist(1-row_C), "ward.D2"); row_ord = row_hc$order

	col_C = cor(X, m="spearman"); diag(col_C) = NA
	col_hc = hclust(as.dist(1-col_C), "ward.D2");
	col_hc = as.hclust(reorder(as.dendrogram(col_hc),
        	row_ord, agglo.FUN=mean))
	col_ord = col_hc$order

	row_n = nrow(row_C); col_n = nrow(col_C)
	row_p = 0.95 * row_n / (row_n + col_n)
	if (!is.null(fname)) {
		png(fname, height= max(min_height, (col_n + row_n) * 12), width= max(min_width, (col_n + row_n) * 12))
	}
	par(mar = c(row_mar,row_mar,1,1),lwd=2, fig = c(0,row_p,0,row_p))
	image.2(row_C[row_ord, row_ord], balance=T, col = cor_palette); box()
	par(mar = c(col_mar,row_mar,1,1), fig = c(row_p,0.95,row_p,0.95), new=T)
	image.2(col_C[col_ord, col_ord], balance=T, col = cor_palette); box()
	par(mar = c(row_mar,row_mar,1,1), fig = c(row_p,0.95,0,row_p), new=T)
	image.2(X[row_ord, col_ord], balance=T, col = X_palette); box()
	par(mar = c(col_mar,1,1,1), fig = c(0.95,1,row_p, 0.95), new=T);
	image(t(seq_along(col_ord)), axes=F, col = col_cols[ colnames(X)[col_ord]]); box()
	par(mar = c(row_mar,1,1,1), fig = c(0.95,1,0,row_p), new=T);
	image(t(seq_along(row_ord)), axes=F, col = row_cols[ rownames(X)[row_ord]]); box()
	if (!is.null(fname)) {
		dev.off()
	}
	list(row_hc = row_hc, col_hc = col_hc)
}


cluster_sep = function(X, ct, dim = "h", names = T) {
	n = length(ct)
	ticks = seq(-1 / (2 * (n-1)),1 + 1 / (2 * (n-1)), length.out = n + 1)
	bcls = cumsum(table(ct)) + 1; names(bcls) = comb_levels
        abline(h = ticks[bcls], lwd = 3)
        mtext(names(bcls), side = 2, las = 2, at = rowMeans(cbind(ticks[c(1,bcls[-length(bcls)])], ticks[bcls])))	
}

sc_pipeline = function(basedir, outdir, batches, index_fn, umi_dir,filt_outliers_on_clusts = T,
	batch_meta_attr = "Amp_batch_ID", amb_epsilon = 0.03, mark_blacklist_terms = c(),
	T_edge_2dproj = 0.05, force_max_deg_2dproj = 8, read_and_clean_only = F, cells = NULL) {

	dir.create(basedir)
	dir.create(outdir)

	scdb = scdb_init(basedir = basedir)

	sc_raw_mat = sc_pipe_build_mat(scdb,
        index_fn = index_fn,
        batch_meta_attr = batch_meta_attr,
        base_dir= umi_dir,
        mars_batches = batches,
        outdir = outdir)

	sc_clean_mat = sc_pipe_mat_to_clean_mat(scdb, sc_raw_mat,
        batch_meta_attr = batch_meta_attr,
        filt_amb_on_clusts = T,
        filt_outliers_on_clusts = filt_outliers_on_clusts,
        remove_outliers_before_clustering = F,
        mark_blacklist_terms = mark_blacklist_terms,
        amb_epsilon = amb_epsilon)
	
	if (!read_and_clean_only) {
		if (!is.null(cells)) {
			sc_clean_mat@cells = cells; sc_clean_mat@ncells = length(cells)
			sc_clean_mat@mat = sc_clean_mat@mat[,cells]
		}
		sc_cl = sc_pipe_mat_to_cmods(scdb, sc_clean_mat,
	        mark_blacklist_terms = mark_blacklist_terms,
	        outdir = outdir,
	        filt_outliers_on_clusts = filt_outliers_on_clusts,
	        mc_fp_metadata_fields = NA,
	        tab_mc_fp_fn = NA)

		sc_2d = sc_pipe_plots(scdb, sc_cl,
	      	outdir = outdir,
	        mark_blacklist_terms = mark_blacklist_terms,
	        T_edge_2dproj = T_edge_2dproj,
	        force_max_deg_2dproj = force_max_deg_2dproj)
	}
	scdb
}

read_large_umis = function(mat_id, bs = 1e4, cells = NULL) { 
	mat = scdb_mat(mat_id)
	if (is.null(cells)) { cells = mat@cells}
	ncells = length(cells)
	umis = NULL
	for (i in seq_len(ncells %/% bs + 1)) {
		from = (i - 1) * bs + 1
		to = min(i * bs, ncells)
		umis = cbind(umis, as.matrix(mat@mat[,cells[from:to]]))
	}
	umis
}

perform_bootstrap = function(scdb, n, outdir, mark_blacklist_terms, boot_ratio = 0.7, T_edge_2dproj = 0.05, force_max_deg_2dproj = 8, 
	n_procs = 16, analysis_only = F) {
	sc_cl = sc_pipe_mat_to_cmods(scdb)
	filename = paste0(scdb@basedir, "/scrdb_data_cocluster.Rda")
	if (!file.exists(filename)) {
	   cc = scc_bootstrap_clusts(sc_cl,
        	scdb = scdb,
		k_cut = get_param("scc_bootstrap_k_cut"),
        	boot_ratio = boot_ratio,
		n_procs = n_procs,
        	min_clust_size = 20,
        	verbose = T,
        	rseed=1)

	   save(cc, file = filename)
	} else { load(filename)}
	if (!analysis_only) {
		bs_cl = scc_coclust_to_cell_modules(scdb, sc_cl, cc, n, 20, method="hclust")
		file.remove(paste0(scdb@basedir, "/scrdb_data_plot.RData"))
		sc_2d = sc_pipe_plots(scdb, bs_cl,
        		outdir = outdir,
	        	mark_blacklist_terms = mark_blacklist_terms,
	        	T_edge_2dproj = T_edge_2dproj,
	        	force_max_deg_2dproj = force_max_deg_2dproj)
	}
}

plot_bootstrapping = function(scdb) {

}

score_on_gene_prograns = function(umis, gene_map, thresh = 0.5) {

	genes = intersect(rownames(umis), names(gene_map))
	umis = umis[genes,]; gene_map = gene_map[genes]
	cell_modules = apply(umis, 2, tapply, gene_map, sum)
	cell_modules = cell_modules[ -nrow(cell_modules),]
	below_med = cell_modules <= apply(cell_modules,1,quantile, thresh)
	cell_bg = apply(matrix(rownames(cell_modules)),1, function(x) 
		{ifelse(sum(umis[,below_med[x,]]) > 0, sum(cell_modules[x, below_med[x,]]) / sum(umis[,below_med[x,]]), 0)})
	names(cell_bg) = rownames(cell_modules)
	cell_uc = colSums(umis)
	zscore = apply(matrix(rownames(cell_modules)),1,
               function(x) { (cell_modules[x,] - cell_bg[x] * cell_uc) / sqrt(cell_bg[x] * (1 - cell_bg[x]) * cell_uc)})

	zscore
}

xy_scatter_genes = function(x,y, bad_genes = c(), disp_genes = c(), cols = c("navyblue", "chocolate2"), text = F, fc = 1, reg = 10, lwd = 1) {
	if (length(disp_genes) == 0) {
		z = log2((y + reg) / (x + reg))
		disp_genes = names(which(abs(z) > fc))
	}
        good_genes = setdiff(names(x), bad_genes)
	lim = log2(c(reg, reg + max(c(x,y))))
        plot(log2(x[good_genes] + reg), log2(y[good_genes] + reg), pch = 20, cex = 2, col = cols[1 + good_genes %in% disp_genes],
        xlim = lim, ylim = lim, axes = F, xlab = "", ylab = "")
        if (text & length(disp_genes) > 0) { text(log2(x[disp_genes] + reg), log2(y[disp_genes] + reg), disp_genes)}
        abline(coef = c(fc,1), lty = 2); abline(coef = c(-fc,1), lty = 2)
        axis(1); axis(2)
}

choose_genes_from_clust = function(mc_id, mat_id, good_clusts = colnames(sc_cl@mc_fp), 
	nms_per_clust = 5, nms_thresh = 5, max_num = Inf, bad_genes = c(), must_haves = c(),
	ord = "none") {

	sc_cl = scdb_mc(mc_id); sc_mat = scdb_mat(mat_id)
	lfp = log2(sc_cl@mc_fp[,good_clusts])
	nms = unique(as.vector(apply(lfp,2,function(x){ head(names(sort(x, decreasing = T)),nms_per_clust)})))
	nms = setdiff(nms, c(bad_genes, names(which(apply(lfp[nms,],1,max) < nms_thresh))))
	nms = union(must_haves, head(names(sort(apply(lfp[nms, ], 1, max),T)), max_num - length(must_haves)))
	if (ord == "hc") {
		nms = nms[ hclust(dist(cor(t(lfp[nms,]))), "ward.D2")$order]
	} else if (ord == "max.col") {
		nms = nms[ order(max.col(lfp[nms,]), rowSums(as.matrix(sc_mat@mat[nms,])))]
	}
	nms
}

sc_to_bulk = function(mc_id, mat_id, comb=NULL, bad_genes = c(), cells = names(comb), min_comb = 0, choose_genes = T, normalize = T,
	g1 = NULL, g2 = NULL) {

	if (is.null(comb)) { cells = union(g1,g2)}	
        sc_cl = scdb_mc(mc_id); sc_mat = scdb_mat(mat_id)
	umis = read_large_umis(mat_id, cells = cells)
	umis_n = sweep(umis,2,colSums(umis),"/") * 1000
	
	if (!normalize) { umis_n = umis}
	if (is.null(comb)) {
		comb = 1 + (cells %in% g1)
		names(comb) = cells
	}
	MAP = as.numeric(factor(comb[cells])); names(MAP) = cells
	if (choose_genes) {
		genes = setdiff(scr_chi_square_diff_genes(umis[,cells], MAP = MAP[cells], fdr = T, pval = 1e-3), bad_genes)
	} else {
		genes = rownames(umis)
	}
	m = t(apply(umis_n[genes, cells],1,tapply,comb[cells],sum))
	sizes = table(comb[cells])
	good_clusts = names(which(sizes >= min_comb))
	m = m[,good_clusts]; sizes = sizes[good_clusts]
	m = sweep(m,2,as.vector(sizes),"/") * min(sizes)
	m
}


create_batch_matrix = function(mc_id, mat_id, comb, fname, batch_iden = "amp_batch_id", cord = colnames(sc_cl@mc_fp), clusts = NULL,
        batch_shades = colorRampPalette(c("white", "navyblue"))(1000), color_batches = F, txt_file = F, norm_by = "col",
        mar = c(5,10,0,20), zlim = c(0,1)) {

        sc_cl = scdb_mc(mc_id); sc_mat = scdb_mat(mat_id)
        cells = intersect(names(sc_cl@mc), names(comb))
        cell_stats = sc_mat@cell_metadata[cells,]
        if (is.null(clusts)) { clusts = sc_cl@mc[cells]}
        comb_levels = names(table(comb))
        if (is.null(levels(comb))) { comb = factor(comb); names(comb) = cells}
        batch_meta = unique(cbind(as.vector(cell_stats[cells,batch_iden]), as.vector(comb)))
        B = batch_meta[,1]
        batch_dist = table(cell_stats[cells,batch_iden], clusts[cells])
        batch_dist = batch_dist[B,]
        bct = factor(batch_meta[,2], levels = comb_levels); names(bct) = batch_meta[,1]
        if (color_batches) {
                cols = sc_cl@colors[ as.numeric(cord)]
        } else {cols = rep("black", length(cord))}
        #write.table(colnames(batch_dist))
	batch_dist = batch_dist[, as.character(cord)]
        if (!is.null(fname)) {
                png(fname, height = max(nrow(batch_dist) * 12, 1000), width = max(ncol(batch_dist) * 20, 1500))
                par(mar = mar)
        }
        if (norm_by == "col") {
                dist_n = sweep(batch_dist,2,colSums(batch_dist),"/")
        } else if (norm_by == "row") {
                dist_n = sweep(batch_dist,1,rowSums(batch_dist),"/")
        } else {
                dist_n = sweep(batch_dist,1,rowSums(batch_dist),"/")
                dist_n = sweep(dist_n,2,colSums(dist_n),"/")
        }
        rownames(dist_n) = paste0(rownames(batch_dist), "(", rowSums(batch_dist), ")")
        colnames(dist_n) = paste0(colnames(batch_dist), "(", colSums(batch_dist), ")")
        image.2(dist_n, col = batch_shades, , zlim = zlim, hct = bct, text = T, text_mat = batch_dist)
        if (!is.null(fname)) {
                dev.off()
        }
        rownames(batch_meta) = batch_meta[,1]
        df = cbind(batch = rownames(batch_dist), source = batch_meta[rownames(batch_dist),2], as.matrix(batch_dist))
        if (txt_file) {
                write.table(df, sep = "\t", quote = F, row.names=F, file = gsub(".png", ".txt", fname))
        }
}

plot_sc_heatmap = function(mc_id, mat_id, nms, clusts = NULL, good_clusts = NULL, cells = NULL, fname = NULL, annotate=F, annotate_what = "both", mar = rep(0,4),
	genes_shades = colorRampPalette(c("white", "orange", "tomato","mediumorchid4", "midnightblue"))(1000), draw_cls = T, normalize = T,
	lwd=1, lty=2) {

	sc_cl = scdb_mc(mc_id); sc_mat = scdb_mat(mat_id)
	if (is.null(good_clusts)) {
		good_clusts = colnames(sc_cl@mc_fp)
	}
	if (is.null(cells)) {
		cells = names(which(sc_cl@mc > 0 & sc_cl@mc %in% good_clusts))
	}
	if (is.null(clusts)) {
		clusts = sc_cl@mc
	}
	#umis = as.matrix(sc_mat@mat[,cells])
	umis = read_large_umis(mat_id, cells=cells)
	umis_n = sweep(umis,2,colSums(umis),"/") * 1000
	if (normalize) {
		foc = log(1 + 7 * umis_n)
	} else {
	        foc = log(1 + 7 * umis)
	}
	cls = cumsum(table(factor(clusts[cells], levels = good_clusts))) / length(cells)	
	if (!is.null(fname)) { 
		png(fname, height = max(2000, 12*length(nms)), width = 3300)
		par(mar = mar)
	}
	cell_ord = cells[order(factor(clusts[cells], levels = good_clusts), sample(length(cells)))]
	IM = foc[nms, cell_ord]
	image(t(IM), col = genes_shades, axes = F)
	zlim = quantile(IM, c(0,1))
	message("zlim: ", zlim[1], " - ", zlim[2])
	if (annotate & (annotate_what == "rows" | annotate_what == "both")) {
		mtext(nms, side = 2, las = 2, at = (1 - seq_along(nms)) / (1 - length(nms)))
	}
        if (annotate & (annotate_what == "columns" | annotate_what == "both")) {
		mtext(names(cls), side = 1, las = 2, at = rowMeans(cbind(c(0,cls[-length(cls)]), cls)))
	}
	if (draw_cls) {
		abline(v = cls, lty = lty, lwd=lwd)
	}
        if (!is.null(fname)) {
		dev.off()
	}
	cell_ord
}

.graph_con_comp = function(amat)
{
        if(nrow(amat) > 500) {
                message("graph_con_comp work on small matrices only")
                retrun(NA)
        }
        diag(amat) = TRUE
        ps = amat/rowSums(amat)
        for(i in 1:20) {
                ps = ps %*% ps
        }
        hc = hclust(dist(ps>0))
        compo = cutree(hc, h=0.01)
        return(compo)
}

reposition_cc = function(sc_2d, sc_cl, good_clusts = colnames(sc_cl@mc_fp), margin = 500, coords = NULL) {
        xlim = quantile(sc_2d@sc_x, c(0,1), na.rm = T) + c(-margin, margin)
        ylim = quantile(sc_2d@sc_y, c(0,1), na.rm = T) + c(-margin, margin)
	clusts = colnames(sc_cl@mc_fp)
	good_cells = names(sc_cl@mc)[ sc_cl@mc %in% good_clusts]
	graph = sc_2d@graph
	graph = graph[ graph$mc1 %in% good_clusts & graph$mc2 %in% good_clusts,]
	graph_mat = matrix(0, length(clusts), length(clusts), dimnames = list(clusts, clusts))
	graph_mat[cbind(as.character(graph$mc1), as.character(graph$mc2))] = 1
	P = .graph_con_comp(graph_mat)
	cols = rainbow(max(P))
	if (is.null(coords)) {
		plot(sc_2d@sc_x[good_cells], sc_2d@sc_y[good_cells], pch = 21, bg = cols[P[sc_cl@mc[good_cells]]], xlim = xlim, ylim = ylim)
        	coords <- unlist(locator(2, type="l"))
	}
	nearest_cc = P[sc_cl@mc[names(which.min( rowSums(cbind(sc_2d@sc_x[good_cells] - coords[1], sc_2d@sc_y[good_cells] - coords[3]) ^ 2)))]]
        new_2d = sc_2d
	new_2d@sc_x = new_2d@sc_x + (P[sc_cl@mc] == nearest_cc) * (coords[2] - coords[1])
        new_2d@sc_y = new_2d@sc_y + (P[sc_cl@mc] == nearest_cc) * (coords[4] - coords[3])
	new_2d@mc_x = tapply(new_2d@sc_x, sc_cl@mc, mean)
	new_2d@mc_y = tapply(new_2d@sc_y, sc_cl@mc, mean)
        list(sc_2d = new_2d, coords = coords)
}

scr_write_models_file = function(mc_id, mat_id, filename = "models.txt") {
	sc_cl = scdb_mc(mc_id); sc_mat = scdb_mat(mat_id)	
	umis = as.matrix(sc_mat@mat[, names(sc_cl@mc)])
        genes = rownames(sc_cl@mc_fp)
        write.table(data.frame(log2(sc_cl@mc_fp), umicount = rowSums(umis[genes,])), col.names = NA, quote = F, sep = "\t", file = filename)
}

proj_sc_on_clusts = function(umis, new_cells, old_cells, clusts, markers, knn_for_cor = 10, K = 10) {

	umis_n = sweep(umis,2,colSums(umis),"/") * 1000
	ranks = rev(1 - (0:knn_for_cor+1)/(knn_for_cor + 2))
	d_old = t(log2(1+7*umis_n))[old_cells ,markers]
	d_new = t(log2(1+7*umis_n))[new_cells ,markers]
	dcor = as.matrix(cor(t(d_new), t(d_old)))

	clusts = clusts[old_cells]
	clusts1 = diag(max(clusts))[,clusts]    # 1 encoding of the clusters, rows= clusters, cols =nodes
	csizes = rowSums(clusts1); csizes[csizes==0] = 1 # add 1 to avoid division by 0
	clusts1 = clusts1/csizes

	m_knn = t(apply(dcor, 1, function(x) ranks[ rank(-x)])) #1 - pmin(rank(-x)/knn_for_cor,1) ))
	m_knn[ is.na(m_knn)] = 0

	t_knn = t(apply(t(dcor), 1, function(x) ranks[ rank(-x)])) #1 - pmin(rank(-x)/knn_for_cor,1) ))
	t_knn[ is.na(t_knn)] = 0

	forward_votes = t(clusts1 %*% t(m_knn)) *1000
	backward_votes = tcrossprod(t(t_knn), clusts1) * 1000
	votes = backward_votes * forward_votes + 1e-10 * forward_votes  # + 1e-10 * forward_votes is for nodes that no one wants
	clusts = apply(votes,1,which.max)
	names(clusts) = rownames(d_new)
	clusts
}

gate_facs = function(facs, xlab, ylab, gate = TRUE, colgate = NULL, roof = 1e4, rect = FALSE, log = "", polygon = NULL, n = 1000) {

  if (is.null(polygon)) {
    disp_facs = facs[gate, ]
    if (roof < dim(disp_facs)[1]) {
      ind = sample(dim(disp_facs)[1], roof)
    } else {
      ind = 1:dim(disp_facs)[1]
    }

    x = disp_facs[ind, xlab]
    y = disp_facs[ind, ylab]

    if (!is.null(colgate)) {
      colgate = colgate[gate, ]
      colgate = colgate[ind, ]
      colors = unbinary(apply(matrix(1:dim(colgate)[1]), MARGIN = 1, FUN = function(x){paste(colgate[x,] * 1, collapse = "")}))
      #smoothScatter(x, y)
      #points(x, y, col = colors, pch=21, cex = 0.4)
      plot(x, y, pch = 21, cex = 0.4, log = log, col = colors + 1)
    } else {
      plot(x, y, pch = 21, cex = 0.4, log = log)
    }

    coords <- locator(n, type="l") # add lines
   C = unlist(coords)
    if (!rect) {
      n = length(coords$x)
      c = rep(0, n*4)
      dim(c) = c(n*2,2)
      c[1:n,] = C
      c[(n+1):(n*2),] = C
    } else {
      n = 4
      c = matrix(0,n*2,2)
      c[1,] = c(C[1], C[3])
      c[2,] = c(C[1], C[4])
      c[3,] = c(C[2], C[4])
      c[4,] = c(C[2], C[3])
      c[5:8,] = c[1:4,]
    }
    lines(c, col = "red")
  } else {
    c = polygon
    n = nrow(c) / 2
  }

  x = facs[, xlab]
  y = facs[, ylab]

  all_pos = apply(matrix(1:n), MARGIN = 1, FUN =
                    function(i){sign( (c[i + 1,1]- c[i,1])*(y-c[i,2]) - (c[i + 1,2]-c[i,2])*(x-c[i,1]) )})
  w = unlist(all_pos)
  within = apply(w < 0, 1, prod) == 1

  return(list(gate = within, polygon = c))

}

scr_chi_square_diff_genes = function(umis, MAP = NULL, g1 = NULL, g2 = NULL, pval, fdr = F) {

  if (is.null(MAP)) {
    MAP = c(rep(1,length(g1)), rep(2, length(g2)))
    names(MAP) = c(g1, g2)
  }
  cells = names(MAP)
  umis = umis[,cells]
  uniform_a = rowSums(umis)/sum(umis)
  exp_count = matrix(uniform_a, ncol = 1) %*% matrix(colSums(umis),1) # exp_counts per cell
  dimnames(exp_count)  = dimnames(umis)
  ex = t(daply(.data= data.frame(cbind(V1 = MAP, t(exp_count)), check.names = F), .(V1), colSums))[-1,]
  obs = t(daply(.data= data.frame(cbind(V1 = MAP, t(umis)), check.names = F), .(V1), colSums))[-1,]

  x2 = rowSums(((obs-ex)^2 )/ex ) # chi^2 with df = ncluster-1

  if (!fdr) {
    sig_genes = x2 > qchisq(1-pval,df= length(unique(MAP)) - 1)
  } else {
    pvals = p.adjust(1 - pchisq(x2, df = length(unique(MAP)) - 1), "fdr")
    sig_genes = pvals < pval
  }
  sig_genes[ is.na(sig_genes)] = F
  return (names(sig_genes)[sig_genes])

}

scr_find_gene_modules = function(umis, nmodules = 30, min_gene_umis = 5, min_cell_umis = 500,
                                 min_var_mean = 1.2, min_cor_within_module = 0.05, rseed = 1) {
  set.seed(rseed)
  u = .downsamp(umis, min_cell_umis)
  rownames(u) = rownames(umis)
  vm = apply(u, 1, var) / rowMeans(u)
  genes = rowSums(u)>min_gene_umis & vm >= min_var_mean
  u = u[genes,]
  u = log2(u+1)
  #   u = u-rowMeans(u)
  corr = cor(t(u))
  hc = hclust(as.dist(1-corr), method = "ward.D2")
  cls = cutree(hc, nmodules)
  #   cor_means = laply(1:nmodules, function(x) {c = corr[cls==x,cls==x]; mean(c[lower.tri(c)])})
  modules = list()
  for(k in 1:nmodules) {
    c = corr[cls==k,cls==k];
    if (mean(c[lower.tri(c)]) > min_cor_within_module) {
      modules[[length(modules)+1]] = names(cls[cls==k])
    }
  }
  return(list(modules = modules, cls = cls, corr = corr))
}


scr_find_outline = function(clusts_2d, reg = 0.7, bw=50, cells = names(clusts_2d@x)) {
	#cells_graph_x = clusts_2d@x[ !is.nan(clusts_2d@x)]
        #cells_graph_y = clusts_2d@y[ !is.nan(clusts_2d@y)]
	xl = c(min(clusts_2d@x[cells], na.rm=T), max(clusts_2d@x, na.rm=T))
        yl = c(min(clusts_2d@y[cells], na.rm=T), max(clusts_2d@y, na.rm=T))

        rngx = xl[2]-xl[1]
        rngy = yl[2]-yl[1]
        pt_reg = merge(seq(xl[1],xl[2],length.out=bw)
                        ,seq(yl[1],yl[2],length.out=bw))
	pt_reg = matrix(unlist(pt_reg),ncol=2)
        back <- bkde2D(rbind(pt_reg, cbind(clusts_2d@x[cells], clusts_2d@y[cells])),
                        b=c(rngx/bw, rngy/bw),
                        gridsize=c(500, 500))

        outline = matrix(0, nrow = 500, ncol = 500)
        Z = back$fhat > quantile(back$fhat, reg)
        outline[ Z[,c(500, 1:499)] != Z] = 1
	outline = melt(outline)

	outline[,1] = back$x1[ outline[,1]]
        outline[,2] = back$x2[ outline[,2]]
        outline = outline[ outline[,3] == 1,1:2]
        neigh =nn2(outline,k = 20)

        ord = c(); i = 1;
        while(length(ord) != nrow(outline)){
                ord = c(ord,i);
                i = neigh$nn.idx[i,which(!(neigh$nn.idx[i,] %in% ord))[1]]
        }
        ord = ord[ !is.na(ord)]
        outline[ord,]
}

plot_gene_2d = function(scl2d,
        gene_nm, w, h, base_dir = NULL, rna_mat = NULL,
        top_color_q = 0.95,
        cex=1.5,
        n_reg_grid=50,
        reg_factor=1,
        low_umi=0,
        mid_umi=1,
        cont_levels = 0,
        bw_bins = 50,
        positive_psize = 0.5,
        negative_psize = 0.5,
        min_rna_to_color = 0,
        outline = NULL,
        return_values = F,
        rna_shades = colorRampPalette(c("white", "white", "lightgray", "darkorange1", "darkgoldenrod1", "darkgoldenrod4")),
        pt_shades = colorRampPalette(c("white", "white", "lightgray", "darkgray", "orange", "burlywood1", "chocolate4")),
        bg_shades = colorRampPalette(c("white", "white", "lightgray", "darkgray", "cyan", "blue1", "blue3")),
        bg_cells = c())
{
        graph_x = scl2d@sc_x
        graph_y = scl2d@sc_y
        if(is.null(rna_mat)) {
                rna = scl2d@scl@scmat@mat[gene_nm,names(graph_x)]
                rna_tot = colSums(as.matrix(scl2d@scl@scmat@mat)[,names(graph_x)])
        } else {
                rna = rna_mat[gene_nm,names(graph_x)]
                rna_tot = colSums(rna_mat[,names(graph_x)])
        }
        med_n = median(rna_tot)
        rna_base = floor(med_n*rna/rna_tot)
        rna = rna_base + sapply(med_n*rna/rna_tot-rna_base, function(p) rbinom(1,1,p))
        rna[is.nan(rna)] = 0

        xl = c(min(graph_x,na.rm=T), max(graph_x, na.rm=T))
        yl = c(min(graph_y,na.rm=T), max(graph_y, na.rm=T))
        rngx = xl[2]-xl[1]
        rngy = yl[2]-yl[1]
        pt_reg = merge(seq(xl[1],xl[2],length.out=n_reg_grid)
                        ,seq(yl[1],yl[2],length.out=n_reg_grid))
        pt_reg = matrix(unlist(pt_reg),ncol=2)
        epsilon=0
        if(!is.null(base_dir)) {
                fnm = gene_nm
                if(nchar(fnm) > 15) {
                        fnm = substr(gene_nm, 1, 15)
                }
                fnm = sub("/", "_", gene_nm)
                fn = sprintf("%s/%s.png", base_dir, fnm)
                png(fn, heigh=h, w=w);
                par(mar=c(0,0,0,0))
        }
        rna_x = unlist(apply(cbind(graph_x,rna),1,function(x) rep(x[1],each=x[2]*reg_factor)))
        rna_y = unlist(apply(cbind(graph_y,rna),1,function(x) rep(x[1],each=x[2]*reg_factor)))
        cells = sapply(strsplit(names(rna_x), "\\."), "[", 1)
        fore <- bkde2D(rbind(pt_reg, cbind(rna_x, rna_y)[!(cells %in% bg_cells),]),
                b=c(rngx/bw_bins, rngy/bw_bins),
                gridsize=c(500, 500))
        if (length(bg_cells) > 0) {
           bg_rna = cbind(rna_x, rna_y)[cells %in% bg_cells,]
        } else {
           bg_rna = cbind(graph_x, graph_y)
        }
        back <- bkde2D(rbind(pt_reg, bg_rna),
                b=c(rngx/bw_bins, rngy/bw_bins),
                gridsize=c(500, 500))
        fore$fhat = fore$fhat * sum(back$fhat)/sum(fore$fhat)
#               smoothScatter(rna_x, rna_y, colramp=rna_col, xlim=xl, ylim=yl)
        lrs = log2(fore$fhat/back$fhat)
        if (length(bg_cells) > 0) {
           zlim = max(abs(quantile(lrs))); zlim = c(-zlim, zlim)
        } else {
                if(median(lrs,na.rm=T)>-1) { #background regul is too dominant
                        lrs = lrs - median(lrs,na.rm=T) - 1
                }
                zlim = c(-1,max(4, max(lrs)))
        }
        message("plot ", gene_nm, " tot ", sum(rna))

#        if (is.null(outline)) {
           xlim = quantile(back$x1, c(0,1))
           ylim = quantile(back$x2, c(0,1))
#        } else {
#           xlim = quantile(outline[,1], c(0,1))
#           ylim = quantile(outline[,2], c(0,1))
#        }
        image(x=back$x1, y=back$x2, z=lrs,zlim=zlim, col=rna_shades(1000),xaxt='n', yaxt='n', xlim = xlim, ylim = ylim)
        if(cont_levels > 0) {
                contour(x=back$x1, y=back$x2, z=lrs, nlevels=cont_levels, lwd=3,add=T, drawlabels=F)
        }

        low_umi_g = ifelse(low_umi == -1, floor(median(rna)), low_umi)
        high_rna_t = as.numeric(quantile(rna, top_color_q))
#               pt_cols = pt_shades[ifelse(rna <= low_umi_g, 1, ifelse(rna <= mid_umi, 2, ifelse(rna <= high_rna_t, 3,4)))]
        pt_val = round(1+999/length(rna) * rank(pmax(rna,min_rna_to_color), ties.method="min"))
        pt_cols = ifelse(names(rna) %in% bg_cells, bg_shades(1000)[pt_val], pt_shades(1000)[pt_val])

#                       cex=ifelse(pt_cols=="white", 0.5, 1),
        points(graph_x,
                graph_y, pch=21,
                bg=pt_cols,
                cex=ifelse(rna>0, positive_psize, negative_psize),
                col=ifelse(pt_cols=="white", "black", "black"),
                lwd=0.5)
        if (is.null(outline)) {
                grid(lwd=3,col="black",lty=1)
        } else {
          points(outline[,1], outline[,2], type = "l", lwd = 4)
        }
        if(!is.null(base_dir)) {
                dev.off()
        }
        if (return_values) {
           return (list(rna_cols = pt_cols, back = back, lrs = lrs))
        }
}

proj_cells_on_clust_graph = function(sc_2d, clusts, blur = 0, omit_clust=-1, nn.idx, main = T, use_top_k = -1)
{
        if(use_top_k == -1) {
                use_top_k = ncol(nn.idx)
        }

        x_cl = sc_2d@mc_x
        y_cl = sc_2d@mc_y

        blurx = blur*(max(x_cl) - min(x_cl))
        blury = blur*(max(y_cl) - min(y_cl))

        omit_x = NA
        omit_y = NA
        omit_clust = as.character(omit_clust)
        clusts = as.character(clusts)
        if(omit_clust != -1) {
                omit_x = x_cl[omit_clust]
                omit_y = y_cl[omit_clust]
                x_cl[omit_clust] = NA
                y_cl[omit_clust] = NA
        }

        px = apply(nn.idx[,1:use_top_k], 1,
                                function(x) ifelse(sum(clusts[x]!=omit_clust)>0,
                                                mean(x_cl[clusts[x]], na.rm=T),
                                                omit_x))
        py = apply(nn.idx[,1:use_top_k], 1,
                                function(x) ifelse(sum(clusts[x]!=omit_clust)>0,
                                                mean(y_cl[clusts[x]], na.rm=T),
                                                omit_y))


        message("Blur x ", blurx, " y ", blury)
	px = px + rnorm(mean=0, sd=blurx, n=length(px))
        py = py + rnorm(mean=0, sd=blury, n=length(py))
        list(x = px, y = py)
}

proj_ds_on_graph = function(id_2d, mc_id, mat_id, markers, umis, K = 50, fn = "new_umis.png", bg_cells = NULL, coords = NULL, knn_for_cor = 100,
	     outline = NULL, bw = 50, cex = 1.5, reg = 10, bg_reg = 1, lwd = 1, 
	     clust_shades = colorRampPalette(c("gray88", "orange", "red"))(101))
{
	
	sc_2d = scdb_mc2d(id_2d); sc_cl = scdb_mc(mc_id); sc_mat = scdb_mat(mat_id)
	cells = names(sc_2d@sc_x)
	old_umis = as.matrix(sc_mat@mat[,cells])
	old_n = sweep(old_umis,2,colSums(old_umis), "/") * 1000
	umis_n = sweep(umis,2,colSums(umis), "/") * 1000
	# create gmod_fp for new ds
	clust_names = names(sc_2d@mc_x)
	if (is.null(coords)) {

		# compute knn
		d_old = t(log2(1+7*old_n[markers, cells]))
		d_new = t(log2(1+7*umis_n[markers,]))
                dcor = as.matrix(cor(t(d_new), t(d_old)))
		m_knn = t(apply(dcor, 1, function(x) 1 - pmin(rank(-x)/knn_for_cor,1) ))
		nn.idx = t(apply(m_knn, 1, function(x) as.numeric(factor(names(tail(sort(x),n=K)), levels = cells))))
		ass_clust = matrix(sc_cl@mc[cells[nn.idx]], ncol = K, dimnames = dimnames(nn.idx))
       		clust = as.numeric(apply(ass_clust,1,function(x){ clust_names[ which.max(table(factor(x, levels = clust_names)))]}))
       		names(clust) = colnames(umis)
	# compute coordinates
		coords = proj_cells_on_clust_graph(sc_2d, sc_cl@mc, 0.02, nn.idx = nn.idx, main = F)	
	}
        return(list( coords = coords, clust = clust))	
}
#	x = sc_2d@sc_x
#      	y = sc_2d@sc_y
#	x_cl = sc_2d@mc_x #nodeRenderInfo(scr_clust_render)$nodeX
#	y_cl = sc_2d@mc_y #nodeRenderInfo(scr_clust_render)$nodeY

#	if (!is.null(bg_cells)) {
#	  wfg = !(colnames(umis) %in% bg_cells)
#	  fg_coords = list(x = coords$x[wfg], y = coords$y[wfg])
#	  bg_coords = list(x = coords$x[!wfg], y = coords$y[!wfg])
#	  fg_clust = clust[ wfg]; bg_clust = clust[ !wfg]
#	} else {
#	  wfg = T
#	  fg_coords = coords
#	  bg_coords = list(x = sc_2d@sc_x, y = sc_2d@sc_y)
#	  fg_clust = clust; bg_clust = sc_cl@mc
#	}

#	xl = c(min(x, na.rm = T), max(x, na.rm = T))
#   	yl = c(min(y, na.rm = T), max(y, na.rm = T))
#        rngx = xl[2]-xl[1]
#	rngy = yl[2]-yl[1]
#        pt_reg = merge(seq(xl[1],xl[2],length.out=bw)
#                        ,seq(yl[1],yl[2],length.out=bw))
#        pt_reg = matrix(unlist(pt_reg),ncol=2)
#	fore <- bkde2D(rbind(pt_reg, 
#			cbind(rep(bg_coords$x, each = bg_reg), rep(bg_coords$y, each = bg_reg)),
#	     	 	cbind(rep(fg_coords$x, each = reg), rep(fg_coords$y, each = reg))),
#                        b=c(rngx/bw, rngy/bw),
#                        gridsize=c(500, 500))
#        back <- bkde2D(rbind(pt_reg, 
#	     		cbind(rep(bg_coords$x, each = bg_reg), rep(bg_coords$y, each = bg_reg))),
#                        b=c(rngx/bw, rngy/bw),
#                        gridsize=c(500, 500))
#        fore$fhat = fore$fhat * sum(back$fhat)/sum(fore$fhat)
#        lrs = log2(fore$fhat/back$fhat)
#	lrs[ is.nan(lrs)] = 0
#        if(median(lrs, na.rm = T)>-1) {
#       #                 lrs = lrs - median(lrs, na.rm = T) - 1
#        }
#        lrs = pmax(pmin(lrs,4,  na.rm = T),-1,  na.rm = T)
        
#	gridx = findInterval(coords$x, seq(min(x, na.rm=T), max(x, na.rm=T), length.out=500), rightmost.closed=T)
#	gridy = findInterval(coords$y, seq(min(y, na.rm=T), max(y, na.rm=T), length.out=500), rightmost.closed=T)
#	lrs_range = (lrs - min(lrs)) / (max(lrs) - min(lrs)) * 100

#	png(fn, height = 2000, width = 2000)
#	plot(x, y, axes = F, xlab = "", ylab = "", type = "n")
#	points(coords$x, coords$y, pch = 20, cex = ifelse(wfg, cex,0), lwd = 2,
#	    col = clust_shades[ 1 + lrs_range[ cbind(gridx, gridy)]])
#	if (!is.null(outline)) {points(outline[,1], outline[,2], type = "l", lwd = 4)}
#	dev.off()
#	return(list( coords = coords, clust = clust))
#}

plot_virtual_facs = function(plate_facs, xlab, ylab, fname, filter = T, gates = list()) {
	png(fname, height=1500, width=1500)
	short_facs = plate_facs[ filter & plate_facs[,xlab] > 0 & plate_facs[,ylab] > 0 & !is.na(plate_facs[,xlab]) & !is.na(plate_facs[,ylab]),]
	k = bkde2D(cbind(log10(short_facs[,xlab]), log10(short_facs[,ylab])), bandwidth = 0.1)
	plot(log10(short_facs[,xlab]), log10(short_facs[,ylab]), pch = 20, cex = 1.5, col = "gray", log = "", xlim = c(1,5), axes = F, xlab = "", ylab = "")
	axis(1, at = 1:5, labels = 10^(1:5))
	axis(2, at = 1:5, labels = 10^(1:5))
	contour(k$x1, k$x2, k$fhat, nlevels = 20, drawlabels=F, lwd = 5, add = T)
	lapply(gates, function(x) lines(log10(x), lwd = 5, col = "red"))
	dev.off()
}

detect_batchy_genes = function(mat_id, gset_id, group_by, marks = NULL, batch_meta_attr = "amp_batch_id", cell_per_batch = 384) {

	sc_mat = scdb_mat(mat_id); sc_gset = scdb_gset(gset_id)
	marks = names(sc_gset@gene_set)
	umis = as.matrix(sc_mat@mat)
	cell_stats = sc_mat@cell_metadata[colnames(umis),]
	m = t(apply(umis[marks,] > 1,1,tapply, as.vector(cell_stats[, batch_meta_attr]), sum))
	m = sweep(m,2, table(as.vector(cell_stats[, batch_meta_attr])), "/") * cell_per_batch
	batch_temp = unique(cbind(as.vector(cell_stats[,c(batch_meta_attr, group_by)]),1))
	batch2area = apply(cbind(as.vector(batch_temp[,group_by]), "") ,1, paste0,collapse="."); names(batch2area) = batch_temp[,batch_meta_attr]
	mm = t(apply(m,1,tapply,batch2area,mean)); m2 = log2((m + 10) / (mm[,batch2area] + 10))
	m2 = log2((m + 10) / (apply(mm,1,mean) + 10))
	m2
}

.downsamp = function (umis, n, replace = F) {
        m = nrow(umis)
        .downsamp_one=function(v,n, replace = F){
                a = tabulate(sample(rep(1:length(v),times=v),replace=replace,size=n),nbins=m)
                return (a)
        }
	ret = apply(umis[, colSums(umis) >= n], 2, .downsamp_one, n)
	rownames(ret) = rownames(umis)
	return(ret)
}

mc_compute_unnorm_fp = function(mc, us, mc_cores = 16)
{
        f_g_cov = rowSums(us) > 10
	cells = intersect(names(mc), colnames(us))
#        mc_cores = get_param("mc_cores")
#        doMC::registerDoMC(mc_cores)
        all_gs = rownames(us[f_g_cov, cells])
        n_g = length(all_gs)
        g_splts = split(all_gs, 1+floor(mc_cores*(1:n_g)/(n_g+1)))
        fnc = function(gs) {
                                        .row_stats_by_factor(us[gs, cells],
                                                                        mc[cells],
                                                                        function(y) {exp(rowMeans(log(1+y)))-1}) }

        clust_geomean = do.call(rbind, mclapply(g_splts, fnc, mc.cores = mc_cores))

        mc_meansize = tapply(colSums(us[,cells]), mc[cells], mean)
        ideal_cell_size = pmin(1000, median(mc_meansize))
        g_fp = t(ideal_cell_size*t(clust_geomean)/as.vector(mc_meansize))
        return(g_fp)
}

mc_compute_fp_amir = function(mc, us, norm_by_mc_size=T, min_total_umi=10)
{
	f_g_cov = rowSums(us) > min_total_umi
	cells = intersect(names(mc), colnames(us))
	if(0) {
		mc_cores = get_param("mc_cores")
		doMC::registerDoMC(mc_cores)
		all_gs = rownames(us[f_g_cov,])
		n_g = length(all_gs)
		g_splts = split(all_gs, 1+floor(mc_cores*(1:n_g)/(n_g+1)))
		fnc = function(gs) {
						.row_stats_by_factor(us[gs,],
										mc@mc,
										function(y) {exp(rowMeans(log(1+y)))-1}) }

		clust_geomean = do.call(rbind, mclapply(g_splts, fnc, mc.cores = mc_cores))
	}
	clust_geomean = t(tgs_matrix_tapply(us[f_g_cov, cells], mc[cells],
									  function(y) {exp(mean(log(1+y)))-1}))
	rownames(clust_geomean) = rownames(us)[f_g_cov]

#	clust_geomean = .row_stats_by_factor(us[f_g_cov,],
#									mc@mc,
#									function(y) {exp(rowMeans(log(1+y)))-1})

	if (norm_by_mc_size) {
		mc_meansize = tapply(colSums(us[,cells]), mc[cells], mean)
		ideal_cell_size = pmin(1000, median(mc_meansize))
		g_fp = t(ideal_cell_size*t(clust_geomean)/as.vector(mc_meansize))
	}
	else {
		g_fp = clust_geomean
	}
	#normalize each gene
	fp_reg = 0.1
	#0.1 is defined here because 0.1*mean_num_of_cells_in_cluster
	#is epxected to be 3-7, which means that we regulairze
	#umicount in the cluster by 3-7.
	g_fp_n = (fp_reg+g_fp)/apply(fp_reg+g_fp, 1, median)

	return(g_fp_n)
}

plot_density_points = function(X, fac, with_box = F, polygon = T, grid = T, cex=1,
		    col_by = "gradient", col_fac = NULL, col_vec = NULL, inflate = 0.5, las = 1, legend = T,
		    pol.col = "white", gradient = colorRampPalette(c("black", "yellow", "orange", "tomato3"))(101),
		    direction = "both", horiz = F, pt.alpha = 0.7, with_points = T, adjust_inflate = T) {

  if(is.null(levels(fac))) {fac = factor(fac)}
  num_fac = as.numeric(fac)
  Y = tapply(X, fac, density)
  Z = matrix(unlist(lapply(Y, function(y) as.numeric(cut(X, y$x)))), ncol = length(Y))
  colnames(Z) = levels(fac)
  Z = Z[ cbind(seq_along(X), fac)]
  y = matrix(unlist(lapply(Y, function(x) x$y)), ncol = length(Y))
  colnames(y) = levels(fac)
  if (adjust_inflate) {inflate = inflate / max(sapply(Y, function(x) max(x$y)))}
  min_val = ifelse(direction == "right", 0, -1)
  max_val = ifelse(direction == "left", 0, 1)
  lim = c(1 + min_val * inflate, max(num_fac) + max_val * inflate) 
  if (!horiz) {plot(num_fac + runif(length(X), min = min_val, max = max_val) * y[cbind(Z,fac)] * inflate, X, type = "n", axes = F, xlab = "", ylab = "")} #, xlim = lim)}
	else {plot(X, num_fac + runif(length(X), min = min_val, max = max_val) * y[cbind(Z,fac)] * inflate, type = "n", axes = F, xlab = "", ylab = "")}  #, ylim = lim)}
  grid()
  if (polygon) {
    if (length(pol.col) == 1) { pol.col = rep(pol.col, length(unique(num_fac)))}
    if (direction == "both") {
    	for (i in rev(seq_along(levels(fac)))) {
      		len = length(Y[[i]]$y)
      		idx = c(seq_len(len), rev(seq_len(len)))
      		if (!horiz) {polygon(i + rep(c(1,-1), each = len) * Y[[i]]$y[idx] * inflate, Y[[i]]$x[idx], col = pol.col[i])}
		else {polygon(Y[[i]]$x[idx], i + rep(c(1,-1), each = len) * Y[[i]]$y[idx] * inflate, col = pol.col[i])}
        }
    } else {
    	for (i in rev(seq_along(levels(fac)))) {
      		sgn = ifelse(direction == "left", -1, 1)
      		if (!horiz) {polygon(i + Y[[i]]$y * inflate * sgn, Y[[i]]$x, col = pol.col[i])}
		else {polygon(Y[[i]]$x, i + Y[[i]]$y * inflate * sgn, col = pol.col[i])}
    	}
    }
  }
  cols = "white"
  if (col_by == "gradient") {
     vals = round(t((t(y) - apply(y,2,min)) / (apply(y,2,max)) - apply(y,2,min)) * 100) + 1
     cols = gradient[vals[cbind(Z,fac)]]
  } else if (col_by == "factor") {
     cols = rainbow(length(levels(col_fac)))
     cols = cols[ as.numeric(col_fac)]
  } else if(col_by == "vector") {
     cols = col_vec
  }
  if (with_points) {
	if (!horiz) {
		points(num_fac + runif(length(X), min = min_val, max = max_val) * y[cbind(Z,fac)] * inflate, X, pch = 21, cex=cex, bg = adjustcolor(cols, alpha.f = pt.alpha))
	} else {
		points(X, num_fac + runif(length(X), min = min_val, max = max_val) * y[cbind(Z,fac)] * inflate, pch = 21, cex=cex, bg = adjustcolor(cols, alpha.f = pt.alpha))
  	}
  }
  axis((!horiz) + 1)
  axis(horiz + 1, at = seq_along(unique(fac)), labels = names(table(fac)), las = las)
  if (with_box) {boxplot(X ~ fac, axes = F, boxwex = rep(0.1, length(table(fac))), add = T, horizontal = horiz, outline = F)}
  if (col_by == "factor") {legend("bottomleft", levels(col_fac), pch = 21, pt.bg = rainbow(length(levels(col_fac))), cex = 2)}
  
}

.row_stats_by_factor = function (data, fact, rowFunction = rowMeans) {
        u = as.character(sort(unique(fact)))
        fact[is.na(fact)] = F
        n=length(u)
        centers = matrix(NA,dim(data)[1], n, dimnames = list(rownames(data), u))
        for (i in u) {
                if(sum(fact==i, na.rm=T)>1) {
                        centers[,i] = rowFunction(data[,fact==i,drop=F])
                } else {
                        centers[,i] = data[,fact==i]
                }
        } # much faster than tapply
        return(centers)
}

mc_compute_fp_abs = function(mc, us, norm_to_med = T)
{
        f_g_cov = rowSums(us) > 10

        mc_cores = 16
        doMC::registerDoMC(mc_cores)
        all_gs = rownames(us[f_g_cov,])
        n_g = length(all_gs)
        g_splts = split(all_gs, 1+floor(mc_cores*(1:n_g)/(n_g+1)))
        fnc = function(gs) {
                                        .row_stats_by_factor(us[gs,],
                                                                        mc,
                                                                        function(y) {exp(rowMeans(log(1+y)))-1}) }
        clust_geomean = do.call(rbind, lapply(g_splts, fnc))
        mc_meansize = tapply(colSums(us), mc, mean)
        ideal_cell_size = pmin(1000, median(mc_meansize))
        g_fp = clust_geomean #t(ideal_cell_size*t(clust_geomean)/as.vector(mc_meansize))
        #normalize each gene
        fp_reg = 0.1
        #0.1 is defined here because 0.1*mean_num_of_cells_in_cluster
        #is epxected to be 3-7, which means that we regulairze
        #umicount in the cluster by 3-7.
	if (norm_to_med) {
	        g_fp = (fp_reg+g_fp)/apply(fp_reg+g_fp, 1, median)
	}
        return(g_fp)
}

plot_two_genes_fp = function(mc_id, ga, gb, log = T, col_mc = c(), show_mc=NULL) {
	sc_cl = scdb_mc(mc_id)
	fp = sc_cl@mc_fp
	if (log) { fp = log2(fp)}
	if (is.null(show_mc)) {show_mc = colnames(fp)}
	show_mc = as.numeric(show_mc)
	if (length(ga) == 1) {a = fp[ga, show_mc]; } else { a = colSums(fp[ga, show_mc])}
	if (length(gb) == 1) {b = fp[gb, show_mc]; } else { b = colSums(fp[gb, show_mc])}
	if (length(col_mc) > 0) {
		cols = ifelse(names(a) %in% col_mc, "yellow", "tomato4")
	} else {cols = sc_cl@colors[ show_mc]}
	plot(a,b, xlab = ga, ylab = gb, pch = 21, cex = 3, bg = cols)
	text(a,b, names(a))
	return(data.frame(a = a, b = b))
}

plot_two_mc_fp = function(mc_id, ca, cb) {
	sc_cl = scdb_mc(mc_id)
	fp = log2(sc_cl@mc_fp)
	a = fp[,ca]; b = fp[,cb]
	plot(a,b, xlab = ca, ylab = cb, type="n")
	grid()
	text(a,b, names(a))
	#return(data.frame(a = a, b = b))
}

mcell_mc_export_tab = function (mc_id, gstat_id = NULL, mat_id = NULL, T_gene_tot = 50,
	  T_fold = 1, metadata_fields = c("batch_set_id"), tab_clust_fp_fn)
{
  mc = scdb_mc(mc_id)
  gstat = scdb_gstat(gstat_id)
  scmat = scdb_mat(mat_id)
  if (is.null(mc)) {
    stop("MC-ERR non existing mc_id ", mc_id, " when trying to export fp table")
  }
  if (is.null(gstat)) {
    stop("MC-ERR non existing gstat id ", gstat_id, " when trying to export fp table")
  }
  if (is.null(scmat)) {
    stop("MC-ERR non existing mat id ", mat_id, " when trying to export fp table")
  }
  if (nrow(mc@color_key) > 0) {
    col2group = as.character(mc@color_key$group)
    names(col2group) = as.character(mc@color_key$color)
    groups = col2group[mc@colors]
  }
  else {
    groups = rep(NA, max(mc@mc))
  }
  mat = as.matrix(scmat@mat)
  out_df = rbind(tapply(colSums(mat[, names(mc@mc)]),
    mc@mc, mean), table(mc@mc), groups, seq_along(groups))
  rownames(out_df) = c("mean_umis", "n_cells", "group", "mc_id")
  if (!is.null(metadata_fields)) {
    for (s in metadata_fields) {
      out_df = rbind(table(scmat@cell_metadata[names(mc@mc),
        s], mc@mc), out_df)
    }
  }
  fp_max = apply(mc@mc_fp, 1, max)
  fp_tot = gstat[intersect(rownames(mc@mc_fp), rownames(gstat)),
    "tot"]
 f = fp_max > T_fold & fp_tot > T_gene_tot
  out_df = rbind(out_df, round(log2(mc@mc_fp[f, ]), 2))
  write.table(out_df, tab_clust_fp_fn, sep = "\t", quote = F, col.names = NA)
}


###################

mcell_mat_rpt_cor_anchors = function (mat_id, gene_anchors, gene_anti = c(), cor_thresh,
    tab_fn, sz_cor_thresh = NA, downsample_n = NA, cells = NULL)
{
    mat = scdb_mat(mat_id)
    if (is.null(mat)) {
        stop("missing mat ", mat_id)
    }
    if (!is.null(cells)) {
        #umis = as.matrix(mat@mat[, cells])
	umis = read_large_umis(mat_id, cells=cells)
    }
    else {
        #umis = as.matrix(mat@mat)
	umis = read_large_umis(mat_id)
    }
    if (is.na(downsample_n)) {
        downsample_n = quantile(colSums(umis), 0.05)
	message(downsample_n)
    }
    #mat_ds = scm_downsamp(umis, downsample_n)
    mat_ds = .downsamp(umis, downsample_n)
    mat_ds = mat_ds[rowSums(mat_ds) > 10, ]
    csize = colSums(umis[, colnames(mat_ds)])
    gcors = data.frame(sz_cor = apply(mat_ds, 1, cor, csize))
    gene_anchors = intersect(gene_anchors, rownames(gcors))
    gene_anti = intersect(gene_anti, rownames(gcors))
    message(length(gene_anchors), " genes")
    for (g in gene_anchors) {
        gcors[, g] = apply(mat_ds, 1, cor, mat_ds[g, ])
    }
    for (g in gene_anti) {
        gcors[, g] = apply(mat_ds, 1, cor, mat_ds[g, ])
    }
    N1 = length(gene_anchors) + 1
    N2 = N1 + length(gene_anti)
    if (length(gene_anchors) == 1) {
        gcors$max = gcors[, 2]
    }
    else {
        gcors$max = apply(gcors[, 2:N1], 1, max)
    }
    if (length(gene_anti) == 0) {
        gcors$neg_max = 0
    }
    else if (length(gene_anti) == 1) {
        gcors$neg_max = gcors[, N2]
    }
    else {
        gcors$neg_max = apply(gcors[, (N1 + 1):N2], 1, max)
    }
    f = !is.na(gcors$max) & gcors$max > cor_thresh & gcors$max >
        gcors$neg_max
    if (!is.na(sz_cor_thresh)) {
        f = f | (gcors$sz_cor > sz_cor_thresh)
    }
    write.table(gcors[f, ], tab_fn, sep = "\t", quote = F)
}

stirlings_lmulti = function(vec) { 
	vec = vec[vec > 0]
	n = sum(vec)
	log_stirl = 0.5 * (log(n) - (length(vec) - 1) * log(2 * pi) - sum(log(vec))) + 
		n * log(n) - sum(vec * log(vec))
	log_stirl
}


explore_gene_set = function(mc_id, mat_id, cells, ct, outdir, prefix, comb,
	downsamp_n = 500, pos_gene = NULL, neg_gene = NULL, batch_attr = "amp_batch_id",
	gradient = colorRampPalette(c("blue", "white", "red"))(1000)) {

        dir.create(outdir)
	sc_cl = scdb_mc(mc_id); sc_mat = scdb_mat(mat_id)
	umis = as.matrix(sc_mat@mat[,cells])
	nms = names(ct)
	mat_ds = as.matrix(scm_downsamp(umis, downsamp_n))

	C = cor(t(log2(1 + mat_ds[nms, ])), m = "spearman")
	diag(C) = NA
	ord = order(ct, -rowSums(umis[nms, ]))
	cls = cumsum(table(ct)) / length(ct)

	modules_n = apply(mat_ds[nms,], 2, tapply, as.vector(ct), sum)
	modules_c = t(apply(modules_n,1,tapply, sc_cl@mc[colnames(modules_n)], mean))
	modules_c = log2((1 + modules_c) / (1 + rowMedians(modules_c)))
	modules_c[ is.infinite(modules_c)] = NA

	chc = hclust(dist(cor(modules_c))); 
	if (!is.null(pos_gene)) {
		chc = as.hclust(reorder(as.dendrogram(chc),
                      modules_c[ct[pos_gene],] - modules_c[ct[neg_gene],],
                      agglo.FUN=mean))
	}
	cord = chc$order; ccls = NA
	modules_c = modules_c[,cord]
	modules_c = modules_c[ order(max.col(modules_c)),]

	message("Plotting modules on clusters")
	png(paste0(outdir, "/", prefix, "_modules_on_clusts.png"), height=2000, width=2000)
	par (mar = c(5,20,0,20), fig = c(0,1,0.8,1))	
	image.2(modules_c, col = gradient, balance = T)
	abline(v = ccls)
	par(mar = c(5,20,0,20), fig = c(0,1,0,0.8), new=T)
	create_batch_matrix(mc_id, mat_id, comb, NULL, batch_attr,
		cord = colnames(modules_c), color_batches=T)
	abline(v = ccls)
	dev.off()

#	clusts = colnames(modules_c)
	df = data.frame(ct = ct, umicount = rowSums(umis[nms, cells]))
	df = df[ order(factor(df$ct, levels = rownames(modules_c)), -df$umicount),]
	write.table(df, sep = "\t", quote = F, col.names = NA, file = paste0(outdir, "/", prefix, "_modules.txt"))

	message("Plotting Truth matrix")
	cls = cumsum(table(factor(ct, levels = rownames(modules_c)))) / length(ct)
	IM = matrix(NA, nrow = length(ct), ncol = ncol(modules_c), dimnames = list(names(ct), colnames(modules_c)))
	shared_genes = intersect(rownames(sc_cl@mc_fp), names(ct))
	IM[shared_genes,] = log2(sc_cl@mc_fp[shared_genes, colnames(modules_c)])
	IM = IM[rownames(df), colnames(modules_c)]
	height = 10 * nrow(IM)
	png(paste0(outdir, "/", prefix, "_trurh.png"), height=height, width=2000)
	image.2(IM, col = gradient, balance = T)
	abline(v=ccls, h = cls)
	mtext(names(cls), las = 2, side = 4, at = rowMeans(cbind(c(0,cls[-length(cls)]),cls)))
	dev.off()

	message("Plotting CT correlation")
	png(paste0(outdir, "/", prefix, "_cor.png"), height=height, width=height)
	par(mar=c(10,10,3,3))
	image.2(C[rownames(df),rownames(df)], col = gradient, balance = T)
	abline(h=cls,v=cls)
	mtext(names(cls), las = 2, side = 3, at = rowMeans(cbind(c(0,cls[-length(cls)]),cls)))
	mtext(names(cls), las = 2, side = 4, at = rowMeans(cbind(c(0,cls[-length(cls)]),cls)))
	dev.off()
	
	colnames(modules_c)
}


mc_colorize_sup_hierarchy = function (mc_id, supmc, supmc_key, gene_key = NULL)
{
    mc = scdb_mc(mc_id)
    if (is.null(mc)) {
        stop("MC-ERR metacell object is not avaialble in scdb, id = ",
            mc_id)
    }
    lfp = log2(mc@mc_fp)
    mc@colors = rep("white", ncol(mc@mc_fp))
    if (!file.exists(supmc_key)) {
        stop("Sup mc key file ", supmc_key, " does not exist")
    }
    key = read.delim(supmc_key, h = T, sep = "\t", stringsAsFactors = F)
    if (class(key)[1] != "data.frame" | length(intersect(c("supid",
        "color", "name"), colnames(key))) != 3) {
        stop("MC-ERR sup id color key must be a data frame with fields supid, color, name")
    }
    for (i in 1:nrow(key)) {
        mcs = supmc[[key$supid[i]]]$mcs
        mc@colors[mcs] = key$color[i]
    }
    color_key = data.frame(gene = rep("", times = nrow(key)),
        group = as.character(key$name), color = as.character(key$color))
    if (!is.null(gene_key)) {
        if (!file.exists(gene_key)) {
            stop("Gene color key file ", gene_key, " does not exist")
        }
        gkey = read.delim(gene_key, h = T, sep = "\t", stringsAsFactors = F)
        if (class(gkey)[1] != "data.frame" | length(intersect(c("name",
            "gene", "color", "T_fold"), colnames(gkey))) != 4) {
            stop("MC-ERR sup id color key must be a data frame with fields gene, name, color, T_fold")
        }
        for (i in 1:nrow(gkey)) {
            gene = gkey$gene[i]
            if (gene %in% rownames(lfp)) {
                T_fold = gkey$T_fold[i]
                mcs = which(lfp[gene, ] > T_fold)
                if (length(mcs) > 0) {
                  mc@colors[mcs] = gkey$color[i]
                  color_key = rbind(as.matrix(color_key), matrix(c(gene = gene,
                    group = gkey$name[i], color = gkey$color[i]),
                    nrow = 1))
                }
            }
        }
    }
    colnames(color_key) = c("gene", "group", "color")
    mc@color_key = as.data.frame(color_key)
    scdb_add_mc(mc_id, mc)
}


import_metacell_structure = function(id, folder=id, all_id="all", bad_genes = c(), url=NULL, mc2d=T, graph=T) {

	if (is.null(url)) {
		sin_tab = read.delim(paste0(folder, "/mc.txt"), stringsAsFactor=F)
		col_tab = read.delim(paste0(folder, "/colors.txt"), stringsAsFactor=F)
		color_key = read.delim(paste0(folder, "/color_key.txt"), stringsAsFactor=F)
                gset_tab = read.delim(paste0(folder, "/gset.txt"), stringsAsFactor=F)
                if (graph) {cgraph = read.delim(paste0(folder, "/cgraph.txt"), stringsAsFactor=T)}
		if (mc2d) {
			mc2d_mc = read.delim(paste0(folder, "/mc2d_mc.txt"), stringsAsFactor=F, row.names=1)
			mc2d_sc = read.delim(paste0(folder, "/mc2d_sc.txt"), stringsAsFactor=F, row.names=1)
			mc2d_graph = read.delim(paste0(folder, "/mc2d_graph.txt"), stringsAsFactor=F)
		}
        } else {
                sin_tab = read.delim(text = getURL(paste0(url, "/mc.txt")), stringsAsFactor=F)
                col_tab = read.delim(text = getURL(paste0(url, "/colors.txt")), stringsAsFactor=F)
                color_key = read.delim(text = getURL(paste0(url, "/color_key.txt")), stringsAsFactor=F)
                gset_tab = read.delim(text = getURL(paste0(url, "/gset.txt")), stringsAsFactor=F)
                if (graph) {cgraph = read.delim(text = getURL(paste0(url, "/cgraph.txt")), stringsAsFactor=T)}
		if (mc2d) {
	                mc2d_mc = read.delim(text = getURL(paste0(url, "/mc2d_mc.txt")), stringsAsFactor=F)
	                mc2d_sc = read.delim(text = getURL(paste0(url, "/mc2d_sc.txt")), stringsAsFactor=F)
        	        mc2d_graph = read.delim(text = getURL(paste0(url, "/mc2d_graph.txt")), stringsAsFactor=F)
		}
        }
	all_mat = scdb_mat(all_id)
        sin_mc = sin_tab[,2]; names(sin_mc) = sin_tab[,1]

        cells = intersect(names(sin_mc), all_mat@cells)
	sin_mc = sin_mc[ cells]
        mcell_mat_ignore_genes(new_mat_id=id, mat_id=all_id, bad_genes, reverse=F)
	mcell_mat_ignore_cells(id, id, cells, reverse=T)

        sin_mat = scdb_mat(id)
        sin_cl = tgMCCov(sin_mc, setdiff(sin_mat@cells, names(sin_mc)), sin_mat)
        colors = rep(NA, ncol(sin_cl@mc_fp))
        col_vec = col_tab[,2]; names(col_vec) = col_tab[,1]
        colors[as.numeric(names(col_vec))] = col_vec
        sin_cl@colors = colors
        sin_cl@color_key = rbind(sin_cl@color_key, color_key)
        scdb_add_mc(id, sin_cl)
	if (mc2d) {
	        mc_x = mc2d_mc$x; names(mc_x) = rownames(mc2d_mc)
	        mc_y = mc2d_mc$y; names(mc_y) = rownames(mc2d_mc)
        	sc_x = mc2d_sc[ cells, "x"]; names(sc_x) = cells
	        sc_y = mc2d_sc[ cells, "y"]; names(sc_y) = cells
	        sin_2d = tgMC2D(mc_id=id, mc_x=mc_x, mc_y=mc_y, sc_x=sc_x, sc_y=sc_y, graph=mc2d_graph)
	        scdb_add_mc2d(id, sin_2d)
	}
	if (graph) {
	        sin_cgraph = tgCellGraph(cgraph, sin_mat@cells)
        	scdb_add_cgraph(id, sin_cgraph)
	}

        gset = gset_tab[,2]; names(gset) = gset_tab[,1]
        sc_gset = tgGeneSets(gset)
        scdb_add_gset(id, sc_gset)
}

