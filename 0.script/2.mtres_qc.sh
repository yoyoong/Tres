# compute Signaling
expression_dir=/sibcb1/bioinformatics/hongyuyang/dataset/Tres/2.tisch_data/1.new_gem_data
expression_list=$(ls ${expression_dir})

expression_list=(SARC_GSE119352_mouse_aPD1aCTLA4)
for expression_filename in ${expression_list[*]}
do
    expression_path=${expression_dir}/${expression_filename}
    expression_tag=$(echo "$expression_filename" | cut -d '.' -f1)
    echo "Processing file: $expression_tag"

    signaling_output_directory=/sibcb1/bioinformatics/hongyuyang/dataset/Tres/2.tisch_data/2.Signaling
    log_path=/sibcb1/bioinformatics/hongyuyang/code/Tres/log/2.tisch_data/2.Signaling/${expression_tag}.log
    rm ${log_path}
    echo "python3 /sibcb1/bioinformatics/hongyuyang/code/Tres/2.mTres/mtres_signaling.py -E $expression_path -D $signaling_output_directory -O $expression_tag" | \
    qsub -q g5.q -N ${expression_tag}.mtres_signaling -V -cwd -o ${log_path} -j y

    sleep 10s
done

# compute Response
celltype=B
expression_dir=/sibcb1/bioinformatics/hongyuyang/dataset/Tres/2.tisch_data/1.new_gem_data
expression_list=$(ls ${expression_dir})
for expression_filename in ${expression_list[*]}
do
    expression_path=${expression_dir}/${expression_filename}
    output_tag=$(echo "${expression_filename}" | cut -d '.' -f1)
    echo "Processing file: ${output_tag}"

    if [ "${celltype}" == "CD8T" ]; then
      response_outdir=/sibcb1/bioinformatics/hongyuyang/dataset/Tres/2.tisch_data/2-1.Proliferation
      genesets_GMT_file=/sibcb1/bioinformatics/hongyuyang/dataset/Tres/0.model_file/SMaRT_geneset.txt
      log_path=/sibcb1/bioinformatics/hongyuyang/code/Tres/log/2.tisch_data/qc/Proliferation/${output_tag}.log
    elif [ "${celltype}" == "Macrophage" ]; then
      response_outdir=/sibcb1/bioinformatics/hongyuyang/dataset/Tres/2.tisch_data/3-2.Polarization
      genesets_GMT_file=/sibcb1/bioinformatics/hongyuyang/dataset/Tres/0.model_file/SMaRT_geneset.txt
      log_path=/sibcb1/bioinformatics/hongyuyang/code/Tres/log/2.tisch_data/3-2.Polarization/${output_tag}.log
    elif [ "${celltype}" == "Neutrophils" ]; then
      response_outdir=/sibcb1/bioinformatics/hongyuyang/dataset/Tres/2.tisch_data/3-3.Neutrophils_response
      genesets_GMT_file=/sibcb1/bioinformatics/hongyuyang/dataset/Tres/0.model_file/Tres_kegg.Neutrophils.txt
      log_path=/sibcb1/bioinformatics/hongyuyang/code/Tres/log/2.tisch_data/3-3.Neutrophils_response/${output_tag}.log
    elif [ "${celltype}" == "NK" ]; then
      response_outdir=/sibcb1/bioinformatics/hongyuyang/dataset/Tres/2.tisch_data/3-4.NK_response
      genesets_GMT_file=/sibcb1/bioinformatics/hongyuyang/dataset/Tres/0.model_file/Tres_kegg.NK.txt
      log_path=/sibcb1/bioinformatics/hongyuyang/code/Tres/log/2.tisch_data/3-4.NK_response/${output_tag}.log
    elif [ "${celltype}" == "NK_act" ]; then
      response_outdir=/sibcb1/bioinformatics/hongyuyang/dataset/Tres/2.tisch_data/3-4-0.NK_act_response
      genesets_GMT_file=/sibcb1/bioinformatics/hongyuyang/dataset/Tres/0.model_file/Tres_kegg.NK_act.txt
      log_path=/sibcb1/bioinformatics/hongyuyang/code/Tres/log/2.tisch_data/3-4-0.NK_act_response/${output_tag}.log
    elif [ "${celltype}" == "B" ]; then
      response_outdir=/sibcb1/bioinformatics/hongyuyang/dataset/Tres/2.tisch_data/3-5.B_response
      genesets_GMT_file=/sibcb1/bioinformatics/hongyuyang/dataset/Tres/0.model_file/Tres_kegg.B.txt
      log_path=/sibcb1/bioinformatics/hongyuyang/code/Tres/log/2.tisch_data/3-5.B_response/${output_tag}.log
    fi
    signature_name_file=/sibcb1/bioinformatics/hongyuyang/dataset/Tres/0.model_file/signature_name_file.txt

    rm ${log_path}
    echo "python3 /sibcb1/bioinformatics/hongyuyang/code/Tres/2.mTres/mtres_response.py -E ${expression_path} -G ${genesets_GMT_file} -S ${signature_name_file} -CT ${celltype} -D ${response_outdir} -O ${output_tag}" | \
    qsub -q b1.q -N ${output_tag}.response -V -cwd -o ${log_path} -j y
    sleep 3s
done