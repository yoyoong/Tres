expression_dir=/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch2_data/1.gem_data
expression_list=$(ls ${expression_dir})
# compute Prolifertion and Signaling
for expression_filename in ${expression_list[*]}
do
    expression_path=${expression_dir}/${expression_filename}
    expression_tag=$(echo "$expression_filename" | cut -d '.' -f1)
    echo "Processing file: $expression_tag"

    prolifertion_output_directory=/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch2_data/2-1.Prolifertion
    if [ ! -d ${prolifertion_output_directory} ]; then
      mkdir ${prolifertion_output_directory}
    fi
    log_path=/sibcb2/bioinformatics2/hongyuyang/code/Tres/log/2.tisch2_data/qc/Prolifertion/${expression_tag}.log
    rm ${log_path}
    echo "python3 /sibcb2/bioinformatics2/hongyuyang/code/Tres/2.mTres/mtres_response.py -E $expression_path -D $prolifertion_output_directory -O $expression_tag" | \
    qsub -q b1.q@fnode014.sibcb.ac.cn -N ${expression_tag}.mtres_response -V -cwd -o ${log_path} -j y

    signaling_output_directory=/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch2_data/2-2.Signaling
    if [ ! -d ${signaling_output_directory} ]; then
      mkdir ${signaling_output_directory}
    fi
    log_path=/sibcb2/bioinformatics2/hongyuyang/code/Tres/log/2.tisch2_data/qc/Signaling/${expression_tag}.log
    rm ${log_path}
    echo "python3 /sibcb2/bioinformatics2/hongyuyang/code/Tres/2.mTres/mtres_signaling.py -E $expression_path -D $signaling_output_directory -O $expression_tag" | \
    qsub -q b1.q@fnode015.sibcb.ac.cn -N ${expression_tag}.mtres_signaling -V -cwd -o ${log_path} -j y

    sleep 10s
done

expression_list=("NSCLC_GSE176021_aPD1")
# compute Prolifertion and Signaling
for expression_filename in ${expression_list[*]}
do
    expression_path=${expression_dir}/${expression_filename}
    expression_tag=$(echo "$expression_filename" | cut -d '.' -f1)
    echo "Processing file: $expression_tag"

    prolifertion_output_directory=/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch2_data/2-1.Prolifertion
    if [ ! -d ${prolifertion_output_directory} ]; then
      mkdir ${prolifertion_output_directory}
    fi
    log_path=/sibcb2/bioinformatics2/hongyuyang/code/Tres/log/2.tisch2_data/qc/Prolifertion/${expression_tag}.log
    rm ${log_path}
    echo "python3 /sibcb2/bioinformatics2/hongyuyang/code/Tres/2.mTres/mtres_response.py -E $expression_path -D $prolifertion_output_directory -O $expression_tag" | \
    qsub -q b1.q@fnode005.sibcb.ac.cn -N ${expression_tag}.mtres_response -V -cwd -o ${log_path} -j y
done

# compute qc
celltype_list=("CD8T")
for celltype in ${celltype_list[*]}
do
    echo "Processing celltype: ${celltype}"

    output_file_directory=/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch2_data/3.qc_result
    if [ ! -d ${output_file_directory} ]; then
      mkdir ${output_file_directory}
    fi

    response_path=/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch2_data/2-1.Prolifertion
    signaling_path=/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch2_data/2-2.Signaling

    log_path=/sibcb2/bioinformatics2/hongyuyang/code/Tres/log/2.tisch2_data/qc/${celltype}.mtres_qc.log
    rm ${log_path}
    echo "python3 /sibcb2/bioinformatics2/hongyuyang/code/Tres/2.mTres/mtres_qc.py -CT ${celltype} -R ${response_path} -S ${signaling_path}" | \
    qsub -q b1.q@fnode008.sibcb.ac.cn  -N ${celltype}.mtres_qc -V -cwd -o ${log_path} -j y
done
