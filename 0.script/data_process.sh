# interaction_cytokine_summary
cytokine_list=("IL6" "FGF2" "CXCL12" "IL1B" "TRAIL" "TGFB1" "TGFB2" "BDNF" "IL17A" "IL12" "IL2" "IL15")
for cytokine in ${cytokine_list[*]}
do
    echo "Processing cytokine: ${cytokine}"

    interaction_path=/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch2_data/4.Interaction/dataset_interaction
    output_file_directory=/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch2_data/4.Interaction/cytokine_summary
    if [ ! -d ${output_file_directory} ]; then
      mkdir ${output_file_directory}
    fi
    output_tag=${cytokine}.summary

    log_directory=/sibcb2/bioinformatics2/hongyuyang/code/Tres/log/2.tisch2_data/data_process
    if [ ! -d ${log_directory} ]; then
      mkdir ${log_directory}
    fi
    log_filename=${log_directory}/${output_tag}.log
    rm ${log_filename}

    echo "python3 /sibcb2/bioinformatics2/hongyuyang/code/Tres/data_process/interaction_cytokine_summary.py -I ${interaction_path} -D ${output_file_directory} -O ${output_tag} -C ${cytokine}" | \
      qsub -q g1.q -N ${output_tag} -V -cwd -o ${log_filename} -j y
done
