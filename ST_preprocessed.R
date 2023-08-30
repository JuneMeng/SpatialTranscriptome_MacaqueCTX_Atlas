suppressMessages(library(Matrix))
suppressMessages(library(Seurat))
suppressMessages(library(data.table))
setwd('/home/share/share_baiyq/spa_matrix_1118') #

args = commandArgs(T)
samplelist_ori <- unique(unlist(lapply(list.files('/home/share/ST_macaque/total_gene_2D/macaque3-MQ121-20221206/', pattern = '.txt'), function(item) strsplit(item, split = '_')[[1]][3]))) # macaque 3 new layer
samplelist <- unlist(lapply(list.files('.', pattern = 'cell-meta-matrix2[.]T*'), function(item) strsplit(item, split = '[.]')[[1]][2]))
samplelist <- setdiff(samplelist_ori, samplelist)
sampletag = samplelist[as.integer(args[1])]
print(sampletag)
### spatial data
# input raw data
filename <- paste0('/home/share/ST_macaque/total_gene_2D/macaque3-MQ121-20221206/total_gene_', sampletag, '_macaque_f001_2D_macaque3-MQ121-20221206.txt') # macaque 3 new layer
expr <- fread(filename, header = TRUE, sep = "\t", data.table = F)
colnames(expr)[colnames(expr) == 'cell_label'] <- 'cell_id'
expr <- expr[expr$cell_id != 0,]
# interpret gene_area
region_map <- fread('/home/share/ST_macaque/total_gene_2D/macaque3-MQ121-20221206/regions-macaque3-MQ121-20221206.csv', header = TRUE, sep = ',', data.table = F) # macaque 3 new layer
region_map <- region_map[region_map$chip == sampletag,]
region_map <- region_map[!(duplicated(region_map$global_id)),]
expr2 <- merge(expr, region_map[,c('origin_name','global_id')], by.x = 'gene_area', by.y = 'global_id')
expr <- expr2
expr2 <- NULL
expr$gene_area <- NULL
colnames(expr)[colnames(expr) == 'origin_name'] <- 'gene_area'
# prefilter
expr <- expr[!(expr$gene_area %in% c('','background')),]
idx <- grep("-l[1-6]", expr$gene_area)	# filter cortex cells
expr <- expr[idx,]
# area
mat_area <- expr[,c('cell_id', 'gene_area')]
mat_area <- mat_area[!(duplicated(mat_area$cell_id)),]
print(dim(mat_area))
# locus
list_xmin <- unlist(tapply(expr[,'x'], expr[,'cell_id'], min))
list_xmax <- unlist(tapply(expr[,'x'], expr[,'cell_id'], max))
list_ymin <- unlist(tapply(expr[,'y'], expr[,'cell_id'], min))
list_ymax <- unlist(tapply(expr[,'y'], expr[,'cell_id'], max))
#list_rx <- unlist(tapply(expr[,'rx'], expr[,'cell_id'], min)) # if no rx/ry should comment
#list_ry <- unlist(tapply(expr[,'ry'], expr[,'cell_id'], min)) # if no rx/ry should comment
mat_loc <- data.frame(xmin = list_xmin, xmax = list_xmax, ymin = list_ymin, ymax = list_ymax) # if rx/ry exists should comment
rownames(mat_loc) <- names(list_xmin)
print(dim(mat_loc))
# volume
mat_loc <- apply(mat_loc, c(1,2), as.numeric)
volume <- unlist(apply(mat_loc, 1, function(entry) (entry[2] - entry[1]) * (entry[4] - entry[3])))
mat_volume <- data.frame(cell_id = rownames(mat_loc), volume = volume)
print(dim(mat_volume))
# create seurat object
binData = expr[,c('gene', 'umi_count', 'cell_id')]
binData[,1] <- factor(binData[,1])
binData[,3] <- factor(binData[,3])
dimnames <- list()
dimnames[[1]] <- levels(binData[,1])
dimnames[[2]] <- levels(binData[,3])
smatrix <- sparseMatrix(as.integer(binData[,1]), as.integer(binData[,3]), x = binData[,2], dimnames = dimnames)
spa_matrix <- CreateSeuratObject(count = smatrix,
                                 project = "Spatial",
                                 min.cells = 3,
                                 min.features = 0)
print(dim(spa_matrix@assays$RNA@counts))
print(nrow(spa_matrix@meta.data))

# percent.mt
MTG = c("ND6","COX3","COX1","ND5","ND4","ND2","ND4L","ATP8","CYTB","COX2","ND3","ATP6","ND1")
spa_matrix <- PercentageFeatureSet(spa_matrix, features = MTG, col.name = "percent.mt")
# meta.data
spa_matrix@meta.data$cell_id <- rownames(spa_matrix@meta.data)
spa_matrix@meta.data <- merge(spa_matrix@meta.data, mat_volume, by = 'cell_id')
mat_loc <- as.data.frame(mat_loc)
print(head(mat_loc))
mat_loc$cell_id <- rownames(mat_loc)
spa_matrix@meta.data <- merge(spa_matrix@meta.data, mat_loc, by = 'cell_id')
spa_matrix@meta.data <- merge(spa_matrix@meta.data, mat_area, by = 'cell_id')
rownames(spa_matrix@meta.data) <- spa_matrix@meta.data$cell_id
spa_matrix@meta.data <- spa_matrix@meta.data[colnames(spa_matrix@assays$RNA@counts),]
# QC (volume 100-3median, nfeature 100, ncount/nfeature Q1, percent.mt 15)
print(table(spa_matrix@meta.data$nFeature_RNA > 100))
print(table(spa_matrix@meta.data$percent.mt < 15))
print(table(spa_matrix@meta.data$volume > (100 / 0.5**2)))
print(table(spa_matrix@meta.data$volume < 3 * median(spa_matrix@meta.data$volume)))
print(table(spa_matrix@meta.data$nCount_RNA/spa_matrix@meta.data$nFeature_RNA > 1.5))
spa_matrix <- subset(spa_matrix, subset = (nFeature_RNA > 100) & (percent.mt < 15))
print(nrow(spa_matrix@meta.data))
print(paste0('filter_ratio:', sampletag, ',',as.character(nrow(mat_loc)), ',', as.character(nrow(spa_matrix@meta.data)), ',', as.character(nrow(spa_matrix@meta.data)/nrow(mat_loc))))

# sct & regressout
spa_matrix <- SCTransform(spa_matrix, vars.to.regress = "percent.mt", verbose = FALSE)
# saverds
saveRDS(spa_matrix, paste0('spa_matrix.', sampletag, '.rds'))
