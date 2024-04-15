# mtres_interaction
celltype=Macrophage
expression_dir=/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/1.gem_data
expression_list=$(ls ${expression_dir})
# expression_list=("NSCLC_GSE176021_aPD1")
for expression_filename in ${expression_list[*]}
do
    expression_path=${expression_dir}/${expression_filename}
    # output_tag=$(echo "$expression_filename" | cut -d '.' -f1-3)
    output_tag=$(echo "$expression_filename" | cut -d '.' -f1)
    echo "Processing file: ${output_tag}"

    if [ "${celltype}" == "CD8T" ]; then
      outdir=/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/5-1.CD8T_Interaction/dataset_interaction
      response_data=/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/3-1.Prolifertion/${output_tag}.csv
      signaling_data=/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/2.Signaling/${output_tag}.csv
      log_path=/sibcb2/bioinformatics2/hongyuyang/code/Tres/log/2.tisch_data/qc/Prolifertion/${output_tag}.log
    elif [ "${celltype}" == "Macrophage" ]; then
      outdir=/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/5-2.Macrophage_Interaction/dataset_interaction
      response_data=/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/3-2.Polarization/${output_tag}.csv
      signaling_data=/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/2.Signaling/${output_tag}.csv
      log_path=/sibcb2/bioinformatics2/hongyuyang/code/Tres/log/2.tisch_data/5-2.Macrophage_Interaction/${output_tag}.log
    fi

    if [ -f ${response_data} ] && [ -f ${signaling_data} ]; then
      rm ${log_path}
      echo "python3 /sibcb2/bioinformatics2/hongyuyang/code/Tres/2.mTres/mtres_interaction.py -E ${expression_path} -R ${response_data} -S ${signaling_data} -D ${outdir} -O ${output_tag}" | \
        qsub -q g5.q  -N ${output_tag} -V -cwd -o ${log_path} -j y
    fi
done


# interaction_cytokine_summary
cytokine_list=("IL6" "FGF2" "CXCL12" "IL1B" "TRAIL" "TGFB1" "TGFB2" "BDNF" "IL17A" "IL12" "IL2" "IL15")
for cytokine in ${cytokine_list[*]}
do
    echo "Processing cytokine: ${cytokine}"

    interaction_path=/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/4.Interaction/dataset_interaction
    output_file_directory=/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/4.Interaction/cytokine_summary
    if [ ! -d ${output_file_directory} ]; then
      mkdir ${output_file_directory}
    fi
    output_tag=${cytokine}.summary

    log_directory=/sibcb2/bioinformatics2/hongyuyang/code/Tres/log/2.tisch_data/data_process
    if [ ! -d ${log_directory} ]; then
      mkdir ${log_directory}
    fi
    log_filename=${log_directory}/${output_tag}.log
    rm ${log_filename}

    echo "python3 /sibcb2/bioinformatics2/hongyuyang/code/Tres/data_process/interaction_cytokine_summary.py -I ${interaction_path} -D ${output_file_directory} -O ${output_tag} -C ${cytokine}" | \
      qsub -q g1.q -N ${output_tag} -V -cwd -o ${log_filename} -j y
done
