expression_list=$(find /sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch2_data/1.gem_data -name "*.csv")
# expression_path=/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch2_data/1.gem_data/NSCLC_GSE176021_aPD1
# expression_tag=NSCLC_GSE176021_aPD1

# compute signaling
for expression_path in ${expression_list[*]}
do
    expression_filename=$(basename "$expression_path")
    # expression_tag=$(echo "$expression_filename" | cut -d '.' -f1-3)
    expression_tag=$(echo "$expression_filename" | cut -d '.' -f1)
    echo "Processing file: $expression_tag"

    signaling_output_directory=/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch2_data/2-2.Signaling
    if [ ! -d ${signaling_output_directory} ]; then
      mkdir ${signaling_output_directory}
    fi
    rm /sibcb2/bioinformatics2/hongyuyang/code/Tres/log/2.tisch2_data/qc/${expression_tag}.mtres_signaling.log
    echo "python3 /sibcb2/bioinformatics2/hongyuyang/code/Tres/mtres_signaling.py -E $expression_path -D $signaling_output_directory -O $expression_tag" | \
    qsub -N ${expression_tag}.mtres_signaling -V -cwd -o /sibcb2/bioinformatics2/hongyuyang/code/Tres/log/2.tisch2_data/qc/${expression_tag}.mtres_signaling.log -j y
done

# compute response
expression_list=$(find /sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch2_data/1.gem_data -name "*.csv")
for expression_path in ${expression_list[*]}
do
    expression_filename=$(basename "$expression_path")
    # expression_tag=$(echo "$expression_filename" | cut -d '.' -f1-3)
    expression_tag=$(echo "$expression_filename" | cut -d '.' -f1)
    echo "Processing file: $expression_tag"

    prolifertion_output_directory=/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch2_data/2-1.Prolifertion
    if [ ! -d ${prolifertion_output_directory} ]; then
      mkdir ${prolifertion_output_directory}
    fi
    rm /sibcb2/bioinformatics2/hongyuyang/code/Tres/log/2.tisch2_data/qc/response/${expression_tag}.mtres_response.log
    echo "python3 /sibcb2/bioinformatics2/hongyuyang/code/Tres/mtres_response.py -E $expression_path -D $prolifertion_output_directory -O $expression_tag" | \
    qsub -N ${expression_tag}.mtres_response -V -cwd -o /sibcb2/bioinformatics2/hongyuyang/code/Tres/log/2.tisch2_data/qc/response/${expression_tag}.mtres_response.log -j y
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

    rm /sibcb2/bioinformatics2/hongyuyang/code/Tres/log/2.tisch2_data/qc/${celltype}.mtres_qc.log
    echo "python3 /sibcb2/bioinformatics2/hongyuyang/code/Tres/mtres_qc.py -CT ${celltype} -R ${response_path} -S ${signaling_path}" | \
    qsub -q g5.q -N ${celltype}.mtres_qc -V -cwd -o /sibcb2/bioinformatics2/hongyuyang/code/Tres/log/2.tisch2_data/qc/${celltype}.mtres_qc.log -j y
done
