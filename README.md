This repository contains the PIC-seq all the needed code and metadata files to analyze and generate figures for Cohen M. & Giladi A. et al. Nature Cancer 2022

In order to run the scripts, download processed data from GSE160903 to the folder output/umi.tab

Unzip and change file names by running these shell commands (in output/umi.tab):
gzip -d *
ls -1 | awk -F'_' '{print $0,$2}' | xargs -n 2 mv

Download published data into their respective folders:
GSE123139	output/published_data/melanoma/
EGAD00001006608 	output/published_data/breast/
GSE135382	output/published_data/pic-seq/

To start analysis, run from the root directory: Rscript scripts/run.r

Please send questions to Amir Giladi: aygoldberg@gmail.com

Required R packages:

Package	version
glmnet	2.0-16
foreach	1.4.4
Matrix	1.2-18
compositions	1.40-2
bayesm	3.1-0.1
energy	1.7-5
robustbase	0.93-3
tensorA	0.36.1
gplots	3.0.1.1
plotrix	3.7-4
plyr	1.8.4
RANN	2.6.1
reshape2	1.4.3
KernSmooth	2.23-15
dendextend	1.9.0
Hmisc	4.2-0
Formula	1.2-3
survival	3.2-3
lattice	0.20-38
tglkmeans	0.2.0
ggrepel	0.8.1
ggplot2	3.3.2.9000
ape	5.2
scales	1.0.0
metacell	0.3.41
tgstat	2.3.5
misha	4.0.10