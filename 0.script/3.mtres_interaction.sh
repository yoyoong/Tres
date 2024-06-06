# mtres_interaction
celltype=CD8T
expression_dir=/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/1.new_gem_data
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
      response_path=/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/3-1.Proliferation/${output_tag}.csv
      signaling_path=/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/2.Signaling/${output_tag}.csv
      log_path=/sibcb2/bioinformatics2/hongyuyang/code/Tres/log/2.tisch_data/5-1.CD8T_Interaction/${output_tag}.log
    elif [ "${celltype}" == "Macrophage" ]; then
      outdir=/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/5-2.Macrophage_Interaction/dataset_interaction
      response_path=/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/3-2.Polarization/${output_tag}.csv
      signaling_path=/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/2.Signaling/${output_tag}.csv
      log_path=/sibcb2/bioinformatics2/hongyuyang/code/Tres/log/2.tisch_data/5-2.Macrophage_Interaction/${output_tag}.log
    elif [ "${celltype}" == "Neutrophils" ]; then
      outdir=/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/5-3.Neutrophils_Interaction
      response_path=/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/3-3.Neutrophils_response/${output_tag}.csv
      signaling_path=/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/2.Signaling/${output_tag}.csv
      log_path=/sibcb2/bioinformatics2/hongyuyang/code/Tres/log/2.tisch_data/5-3.Neutrophils_Interaction/${output_tag}.log
    elif [ "${celltype}" == "NK" ]; then
      outdir=/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/5-4.NK_Interaction
      response_path=/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/3-4.NK_response/${output_tag}.csv
      signaling_path=/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/2.Signaling/${output_tag}.csv
      log_path=/sibcb2/bioinformatics2/hongyuyang/code/Tres/log/2.tisch_data/5-4.NK_Interaction/${output_tag}.log
    fi


    if [ -f ${response_path} ] && [ -f ${signaling_path} ]; then
      # if [ ! -f ${log_path} ]; then
        rm ${log_path}
        echo "python3 /sibcb2/bioinformatics2/hongyuyang/code/Tres/2.mTres/mtres_interaction.py -E ${expression_path} -R ${response_path} -S ${signaling_path} -D ${outdir} -O ${output_tag} -CT ${celltype}" | \
          qsub -q g5.q -N ${output_tag} -V -cwd -o ${log_path} -j y
        sleep 30s
      # fi
    fi
done


# interaction_cytokine_summary
celltype_list=("CD8T" "Macrophage" "Neutrophils" "NK")
for celltype in ${celltype_list[*]}
do
    echo "Processing celltype: ${celltype}"

    log_directory=/sibcb2/bioinformatics2/hongyuyang/code/Tres/log/2.tisch_data/interaction_cytokine_summary
    if [ ! -d ${log_directory} ]; then
      mkdir ${log_directory}
    fi
    log_filename=${log_directory}/${celltype}.log
    rm ${log_filename}

    echo "python3 /sibcb2/bioinformatics2/hongyuyang/code/Tres/2.mTres/interaction_cytokine_summary.py -CT ${celltype}" | \
      qsub -q g5.q -N ${celltype} -V -cwd -o ${log_filename} -j y
done
