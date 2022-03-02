id_s = paste0(id, "_singlets")
id_d = paste0(id, "_PIC")

sc_mat = scdb_mat(id)
sin_2d = scdb_mc2d(id_s); sin_cl = scdb_mc(paste0(id_s, "_f")); sin_mat = scdb_mat(id_s)
db_mat = scdb_mat(id_d)

cells = union(names(sin_cl@mc), db_mat@cells)
umis = as.matrix(sc_mat@mat[,cells])
cell_stats = sc_mat@cell_metadata[cells,]
fp = sin_cl@mc_fp
lfp = log2(sin_cl@mc_fp)

dir.create("revision_figs")
outdir = paste0("revision_figs/figure1")
supdir = paste0("revision_figs/figureS1")
dir.create(outdir)
dir.create(supdir)

sin_stats = sin_mat@cell_metadata[names(sin_cl@mc),]
sin_umis = as.matrix(sin_mat@mat[, names(sin_cl@mc)])
sin_n = sweep(sin_umis,2,colSums(sin_umis),"/") * 1000

color_scheme = sin_cl@color_key
color2name = as.vector(color_scheme$group); names(color2name) = color_scheme$color
name2color = as.vector(color_scheme$color); names(name2color) = color_scheme$group
sin_names = color2name[ sin_cl@colors[ sin_cl@mc]]; names(sin_names) = names(sin_cl@mc)

