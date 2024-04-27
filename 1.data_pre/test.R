rawdata_path = "/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/4.Macrophage_analysis/TCGA/raw_RData/"
output_path = "/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/4.Macrophage_analysis/TCGA/gem/"
annatation_path = "/sibcb2/bioinformatics/KnowledgeBase/Firehose_Methylation/RnBeads_450K_hg19_Probes_GR.RData"

rawdata_list = list.files(rawdata_path)
for (rawdata_shortname in rawdata_list) {
  rawdata_tag = strsplit(rawdata_shortname, "\\.")[[1]][1]
  rawdata_filename = paste0(rawdata_path, rawdata_shortname)
  rawdata = readRDS(rawdata_filename)

  rawdata = get(rawdata_tag)
  gem_filename = paste0(output_path, rawdata_tag, '.csv')
  write.csv(rawdata, gem_filename, row.names=TRUE)

  print(rawdata_shortname)
}

rawdata_tag = "ACC"
load('/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/4.Macrophage_analysis/TCGA/raw_RData/ACC.RData')
rawdata = get(rawdata_tag)
gem_filename = paste0(output_path, rawdata_tag, '.csv')
write.csv(rawdata, gem_filename, row.names=TRUE)


library(Matrix)
gem_tag = "exprmat"
load('/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2-1.neutrophil_data/1.raw_data/NEU_exprmat.rda')
gem_sparse = get(gem_tag)
gem_matrix = as.matrix(gem_sparse)
gem_filename = '/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2-1.neutrophil_data/2.gem_data/gem.csv'
write.csv(gem_matrix, gem_filename, row.names = TRUE)