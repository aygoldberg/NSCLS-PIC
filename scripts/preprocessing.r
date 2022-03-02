set.seed(1111)
library(metacell)
library(scales)
library(ape)
library(ggrepel)
library(tglkmeans)
library(Hmisc)
library(dendextend)
require("KernSmooth")
require("reshape2")
require("RANN")
require("plyr")
require("plotrix")
require("gplots")
require("parallel")
library("compositions")
require(glmnet)

source("scripts/metacell_functions.r")
source("scripts/pic_parser.r")

dir.create("saved_work")
dir.create("figures")
scdb_init("saved_work", force_reinit=T)

##################################
# Import the entire count matrix
message ("Preprocessing")

metadata = read.delim("annotations/TableS1.txt", stringsAsFactor=F)
write.table(metadata[ metadata$organism == "Human",], quote=F, sep = "\t", row.names=F, file = "annotations/metadata_human.txt")
write.table(metadata[ metadata$organism == "Mouse",], quote=F, sep = "\t", row.names=F, file = "annotations/metadata_mouse.txt")

mcell_import_multi_mars("all_human", "annotations/metadata_human.txt", "output/umi.tab/", force = T)
mcell_import_multi_mars("all_mouse", "annotations/metadata_mouse.txt", "output/umi.tab/", force = T)

human_mat = scdb_mat("all_human")
human_stats = human_mat@cell_metadata
singlets = rownames(human_stats)[ human_stats$PIC == "Singlets"]
doublets = rownames(human_stats)[ human_stats$PIC == "PIC"]

mcell_mat_ignore_cells("human_singlets", "all_human", singlets, reverse=T)
mcell_mat_ignore_cells("human_PIC", "all_human", doublets, reverse=T)

mouse_mat = scdb_mat("all_mouse")
mouse_stats = mouse_mat@cell_metadata
singlets = rownames(mouse_stats)[ mouse_stats$PIC == "Singlets"]
doublets = rownames(mouse_stats)[ mouse_stats$PIC == "PIC"]

mcell_mat_ignore_cells("mouse_singlets", "all_mouse", singlets, reverse=T)
mcell_mat_ignore_cells("mouse_PIC", "all_mouse", doublets, reverse=T)

mat = scdb_mat("all_human")
nms = c(rownames(mat@mat), rownames(mat@ignore_gmat))
bad_genes = unique(c(grep("ERCC", nms, v=T), grep("^MT-", nms, v=T)))

import_metacell_structure("human_singlets", "import/human/", "human_singlets", bad_genes)
mcell_import_multi_mars("melanoma_all", "annotations/metadata_melanoma.txt", "output/published_data/melanoma", force = T)
import_metacell_structure("melanoma", "import/melanoma/", "melanoma_all", bad_genes, mc2d=F, graph=F)

#########
# mouse

mat = scdb_mat("all_mouse")
nms = c(rownames(mat@mat), rownames(mat@ignore_gmat))
bad_genes = unique(c(grep("ERCC", nms, v=T), grep("^mt-", nms, v=T)))
import_metacell_structure("in_vitro", "import/in_vitro/", "mouse_singlets", bad_genes, mc2d=F, graph=F)
import_metacell_structure("ot2", "import/ot2/", "mouse_singlets", bad_genes, mc2d=F, graph=F)
import_metacell_structure("ot2_ln", "import/ot2_ln/", "mouse_singlets", bad_genes, mc2d=F, graph=F)
import_metacell_structure("virus", "import/virus/", "mouse_singlets", bad_genes, mc2d=F, graph=F)
import_metacell_structure("ot2_apd-1", "import/ot2_apd-1/", "mouse_singlets", bad_genes, mc2d=F, graph=F)

mcell_import_multi_mars("pic-seq_all", "annotations/metadata_pic-seq.txt", "output/published_data/pic-seq", force = T)

#########################################################
#							#
# Import the breast dataset				#
#							#
#########################################################

indir = "output/published_data/breast/"
cohort2 = readRDS(paste0(indir, "/1867-counts_cells_cohort2.rds"))
cohort1 = readRDS(paste0(indir, "/1864-counts_tcell_cohort1.rds"))
cohort1_stats = read.delim(paste0(indir, "/1872-BIOKEY_metaData_cohort1_web.csv"), sep = ",", stringsAsFactor=F, row.names=1)
cohort2_stats = read.delim(paste0(indir, "/1871-BIOKEY_metaData_cohort2_web.csv"), sep = ",", stringsAsFactor=F, row.names=1)
cohort_stats = rbind(cohort1_stats, cohort2_stats)
cells = intersect(rownames(cohort_stats)[ cohort_stats$cellType == "T_cell"], union(colnames(cohort1), colnames(cohort2)))
genes = union(rownames(cohort1), rownames(cohort2))
data = matrix(0, nrow = length(genes), ncol = length(cells), dimnames = list(genes, cells))

data[ rownames(cohort1), intersect(cells, colnames(cohort1))] = as.matrix(cohort1[, intersect(cells, colnames(cohort1))])
data[ rownames(cohort2), intersect(cells, colnames(cohort2))] = as.matrix(cohort2[, intersect(cells, colnames(cohort2))])
cell_metadata = cohort_stats[ colnames(data),]

scdb_add_mat("breast_all", tgScMat(mat=data, cell_metadata=cell_metadata))
mat = scdb_mat("breast_all")
nms = c(rownames(mat@mat), rownames(mat@ignore_gmat))
bad_genes = unique(c(grep("ERCC", nms, v=T), grep("^MT-", nms, v=T)))
import_metacell_structure("breast", "import/breast/", "breast_all", bad_genes, mc2d=F, graph=F)

##################


