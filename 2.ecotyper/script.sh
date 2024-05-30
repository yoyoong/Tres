options("repos" = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

install.packages("ComplexHeatmap")
install.packages("NMF")
install.packages("RColorBrewer")
install.packages("cluster")
install.packages("circlize")
install.packages("cowplot")
install.packages("data.table")
install.packages("doParallel")
install.packages("grid")
install.packages("reshape2")
install.packages("viridis")
install.packages("config")
install.packages("argparse")
install.packages("colorspace")
install.packages("plyr")
install.packages("Biobase")



gem_path='/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/5.Analysis_data/Braun2020/Braun2020.Express.tsv'
/sibcb/program/install/r-4.1/bin/Rscript EcoTyper_recovery_bulk.R -d Carcinoma -m ${gem_path} -o RecoveryOutput
