# 1、找出136套数据所有细胞类型的markers（1vn）
suppressWarnings(library(Seurat))

rawdata_dir = "/sibcb1/bioinformatics/hongyuyang/dataset/Tres/2.tisch_count_data2"
dataset_markers_dir = "/sibcb1/bioinformatics/hongyuyang/dataset/Tres/2.tisch_data/6.Marker/DatasetMarker2"
rawdata_list = list.files(rawdata_dir, pattern = "\\.csv$")
for (i in 1:length(rawdata_list)) {
    rawdata_filename = rawdata_list[i]
    print(paste0(rawdata_filename, " process start."))
    markers_filename = paste0(dataset_markers_dir, "/", rawdata_filename)
    if (file.exists(markers_filename)) {
        print(paste0(markers_filename, " exists."))
        next
    }

    rawdata_path = paste0(rawdata_dir, "/", rawdata_filename)
    rawdata_df = read.delim(file=rawdata_path, header = TRUE, sep = ",", row.names = 1, check.names=FALSE)
    seurat_data = CreateSeuratObject(counts = rawdata_df, names.field = 1, names.delim = '\\.')
    # result = try({
    #     immune_seurat = subset(seurat, orig.ident %in% c("CD8T", "Mono/Macro", "Neutrophils", "NK", "B"))
    # }, silent = TRUE)
    # if (inherits(result, "try-error")) {
    #     print(paste0(rawdata_filename, " has no immune cell."))
    #     next
    # }

    seurat_data = NormalizeData(seurat_data)
    # all cell type: positive markers
    dataset_markers <- FindAllMarkers(seurat_data, only.pos=T)
    if (length(dataset_markers) > 0) {
        write.table(dataset_markers[c(6,7,2,3,4,1,5)], file=markers_filename, row.names=F, col.names=T, sep=',', quote=F)
    }
    print(paste0(rawdata_filename, " process end."))
}