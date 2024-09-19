library("optparse")
library(Seurat)

option_list = list(make_option("--tag", type = "character", default = NULL, help = "Output tag"))
args = parse_args(OptionParser(option_list=option_list))
tag = args$tag

# 写入CSV文件
output_path = file.path("/sibcb1/bioinformatics/hongyuyang/dataset/Tres/2.tisch_count_data", paste0(tag, ".csv"))
if(file.exists(output_path)){
  print("文件已存在")
  q()
}

rds_path = file.path("/sibcb1/bioinformatics/shijiantao/rasc/RDS", paste0(tag, ".RDS"))
rds_Seurat = readRDS(rds_path)
count_data <- as.data.frame(rds_Seurat@assays$RNA@counts)

write.csv(count_data, output_path, row.names=TRUE)
print(paste(tag, "process end!"))