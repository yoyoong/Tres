dataset_tag=NSCLC_GSE176021_aPD1
input_file=/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch2_data/1.gem_data/${dataset_tag}.csv
output_file_directory=/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch2_data/1.gem_data/${dataset_tag}

rm ${output_file_directory}
mkdir ${output_file_directory}

# 计算总列数
total_columns=$(head -n 1 "$input_file" | tr ',' '\n' | wc -l)

# 计算需要创建的文件数
split_size=10000
num_files=$(( (total_columns + ${split_size} - 1) / ${split_size} ))

# 循环拆分文件
for (( i=0; i<num_files; i++ )); do
    # 计算起始和结束列索引
    start_column=$(( i * ${split_size} + 2 ))
    end_column=$(( (i + 1) * ${split_size} + 1))

    # 构造输出文件名
    output_file=${output_file_directory}/${dataset_tag}_$((i + 1)).csv

    # 提取对应列并写入输出文件
    cut -d ',' -f 1,$start_column-$end_column "$input_file" > "$output_file"

    echo "${i} Processing end."
done