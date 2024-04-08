dataset_tag=NSCLC_GSE176021_aPD1
input_file=/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/1.gem_data/${dataset_tag}.csv
output_file_directory=/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/1.gem_data/${dataset_tag}

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



# 输入CSV文件名
input_file=/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/1.gem_data/NSCLC_GSE127471.csv
output_directory=/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/1.gem_data/NSCLC_GSE127471
# 读取列名
header=$(head -n 1 "$input_file")
# 以逗号为分隔符将列名分割成数组
IFS=',' read -r -a header_list <<< "$header"

# 遍历列名
i=2
for column in "${header_list[@]}"; do
    # 提取第一个.之前的字段作为分组依据
    group_field=$(echo ${column} | cut -d '.' -f 1)

    # 判断是否已经创建了该分组文件，若未创建则创建新文件并写入第一列
    group_output_file="${output_directory}/${group_field}.csv"
    if [ ! -e ${group_output_file} ]; then
        cut -d ',' -f 1 ${input_file} >> ${group_output_file}
    fi

    # 写入当前列数据到对应分组文件
    column_data=$(awk -F ',' -v col="$i" '{print ","$col}' ${input_file} | tail -n +2)
    awk -F ',' -v col="$i" '{print $col}' ${input_file} | tail -n +2 >> ${group_output_file}
    awk -F ',' -v col="$i" '{print $col}' ${input_file} | tail -n +2 | paste -d ',' ${group_output_file} - > ${group_output_file}

    pipe=$(mktemp -u)
    mkfifo ${pipe}
    echo ${column_data} > ${pipe} &

    paste -d ',' ${pipe} ${group_output_file} > ${group_output_file}
    i=$((i + 1))
done