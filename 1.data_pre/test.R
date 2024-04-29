library(Matrix)
gem_tag = "exprmat"
load('/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/1.neutrophil_data/Gao2024/NEU_exprmat.rda')
gem_sparse = get(gem_tag)
gem = as.matrix(gem_sparse)

sample_list <- list()
sample_flag <- list()
for (column in colnames(gem)) {
    column_list <- unlist(strsplit(column, "_"))
    if (!grepl("OA902", column)) {
        sample_flag <- c(sample_flag, TRUE)
        if (length(column_list) == 2) {
            if (grepl("OA089", column)) {
                sample_list <- c(sample_list, paste(column_list[1], unlist(strsplit(column, "-"))[length(unlist(strsplit(column, "-")))]))
            } else {
                sample_list <- c(sample_list, column_list[1])
            }
        } else if (length(column_list) == 3 || length(column_list) == 4 || length(column_list) == 5 || length(column_list) == 6) {
            sample_list <- c(sample_list, paste(column_list[1:(length(column_list) - 1)], collapse = "_"))
        } else if (length(column_list) == 8 || length(column_list) == 9) {
            sample_list <- c(sample_list, paste(c(column_list[1:4], column_list[length(column_list) - 1], column_list[length(column_list)]), collapse = "_"))
        } else if (length(column_list) == 10) {
            sample_list <- c(sample_list, column_list[1])
        }
    } else {
        sample_flag <- c(sample_flag, FALSE)
    }
}
gem_filtered <- gem[, sample_flag]
colnames(gem_filtered) <- sapply(paste("Neutrophils", sample_list, colnames(gem_filtered), sep = ".0"), identity)

column_sums <- colSums(gem_filtered)
non_zero_columns <- which(column_sums != 0)
for (i in non_zero_columns) {
  gem_filtered[, i] <- gem_filtered[, i] / gem_filtered(gem[, i]) * 10000
}
gem_filtered <- log(gem_filtered + 1, base = 2)
row_means <- rowMeans(gem_filtered)
gem_filtered <- sweep(gem_filtered, 1, row_means, "-")

gem_filename = '/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/1.neutrophil_data/Gao2024/new_gem.csv'
write.csv(gem_filtered, gem_filename, row.names = TRUE)