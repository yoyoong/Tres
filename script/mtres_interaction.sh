expression_dir=/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch2_data/1.gem_data
expression_list=$(ls ${expression_dir})
expression_list=("BRCA_GSE161529" "ESCA_GSE160269" "Glioma_GSE148842" "HNSC_GSE139324" "MM_GSE161801" "NPC_GSE162025" "NSCLC_GSE131907" "NSCLC_GSE176021_aPD1" "PRAD_GSE172301" "UVM_GSE139829")

for expression_filename in ${expression_list[*]}
do
    expression_path=${expression_dir}/${expression_filename}
    expression_tag=$(echo "$expression_filename" | cut -d '.' -f1)
    echo "Processing file: $expression_tag"

    output_file_directory=/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch2_data/4.Interaction
    if [ ! -d ${output_file_directory} ]; then
      mkdir ${output_file_directory}
    fi

    response_data=/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch2_data/2-1.Prolifertion/${expression_tag}.csv
    signaling_data=/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch2_data/2-2.Signaling/${expression_tag}.csv

    if [ -f ${response_data} ] && [ -f ${signaling_data} ]; then
      log_directory=/sibcb2/bioinformatics2/hongyuyang/code/Tres/log/2.tisch2_data/mtres_interaction
      if [ ! -d ${log_directory} ]; then
        mkdir ${log_directory}
      fi
      log_filename=${log_directory}/${expression_tag}.log
      rm ${log_filename}
      echo "python3 /sibcb2/bioinformatics2/hongyuyang/code/Tres/mtres_interaction.py -E ${expression_path} -R ${response_data} -S ${signaling_data} -D ${output_file_directory} -O ${expression_tag}" | \
        qsub -q b1.q -N ${expression_tag} -V -cwd -o ${log_filename} -j y
    fi
done
