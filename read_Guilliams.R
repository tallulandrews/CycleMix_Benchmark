Guiliams_meta <- read.delim("/home/tandrew6/scratch/CycleMix_Manuscript/Guilliams/annot_humanMyeloid.csv", sep=",", header=TRUE)

require(Seurat)
Guilliams_counts <- Read10X("/home/tandrew6/scratch/CycleMix_Manuscript/Guilliams/rawData_human/countTable_human/")


