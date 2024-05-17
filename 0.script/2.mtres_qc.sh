# compute Signaling
expression_dir=/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/1.gem_data
expression_list=$(ls ${expression_dir})
for expression_filename in ${expression_list[*]}
do
    expression_path=${expression_dir}/${expression_filename}
    expression_tag=$(echo "$expression_filename" | cut -d '.' -f1)
    echo "Processing file: $expression_tag"

    signaling_output_directory=/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/2-2.Signaling
    log_path=/sibcb2/bioinformatics2/hongyuyang/code/Tres/log/2.tisch_data/qc/Signaling/${expression_tag}.log
    rm ${log_path}
    echo "python3 /sibcb2/bioinformatics2/hongyuyang/code/Tres/2.mTres/mtres_signaling.py -E $expression_path -D $signaling_output_directory -O $expression_tag" | \
    qsub -q b1.q@fnode015.sibcb.ac.cn -N ${expression_tag}.mtres_signaling -V -cwd -o ${log_path} -j y

    sleep 10s
done

# compute Response
celltype=NK
expression_dir=/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/1.new_gem_data
expression_list=$(ls ${expression_dir})
for expression_filename in ${expression_list[*]}
do
    expression_path=${expression_dir}/${expression_filename}
    output_tag=$(echo "${expression_filename}" | cut -d '.' -f1)
    echo "Processing file: ${output_tag}"

    if [ "${celltype}" == "CD8T" ]; then
      response_outdir=/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/2-1.Prolifertion
      genesets_GMT_file=/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/0.model_file/SMaRT_geneset.txt
      log_path=/sibcb2/bioinformatics2/hongyuyang/code/Tres/log/2.tisch_data/qc/Prolifertion/${output_tag}.log
    elif [ "${celltype}" == "Macrophage" ]; then
      response_outdir=/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/3-2.Polarization
      genesets_GMT_file=/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/0.model_file/SMaRT_geneset.txt
      log_path=/sibcb2/bioinformatics2/hongyuyang/code/Tres/log/2.tisch_data/3-2.Polarization/${output_tag}.log
    elif [ "${celltype}" == "Neutrophils" ]; then
      response_outdir=/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/3-3.Neutrophils_response
      genesets_GMT_file=/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/0.model_file/Tres_kegg.Neutrophils.txt
      log_path=/sibcb2/bioinformatics2/hongyuyang/code/Tres/log/2.tisch_data/3-3.Neutrophils_response/${output_tag}.log
    elif [ "${celltype}" == "NK" ]; then
      response_outdir=/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/3-4.NK_response
      genesets_GMT_file=/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/0.model_file/Tres_kegg.NK.txt
      log_path=/sibcb2/bioinformatics2/hongyuyang/code/Tres/log/2.tisch_data/3-4.NK_response/${output_tag}.log
    fi
    signature_name_file=/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/0.model_file/signature_name_file.txt

    rm ${log_path}
    echo "python3 /sibcb2/bioinformatics2/hongyuyang/code/Tres/2.mTres/mtres_response.py -E ${expression_path} -G ${genesets_GMT_file} -S ${signature_name_file} -CT ${celltype} -D ${response_outdir} -O ${output_tag}" | \
    qsub -q g5.q -N ${output_tag}.response -V -cwd -o ${log_path} -j y
    # sleep 10s
done

# compute qc
celltype_list=("CD8T")
for celltype in ${celltype_list[*]}
do
    echo "Processing celltype: ${celltype}"

    output_file_directory=/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/4.qc_result
    if [ ! -d ${output_file_directory} ]; then
      mkdir ${output_file_directory}
    fi

    response_path=/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/2-1.Prolifertion
    signaling_path=/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/2-2.Signaling

    log_path=/sibcb2/bioinformatics2/hongyuyang/code/Tres/log/2.tisch_data/qc/${celltype}.mtres_qc.log
    rm ${log_path}
    echo "python3 /sibcb2/bioinformatics2/hongyuyang/code/Tres/2.mTres/mtres_qc.py -CT ${celltype} -R ${response_path} -S ${signaling_path}" | \
    qsub -q b1.q@fnode008.sibcb.ac.cn  -N ${celltype}.mtres_qc -V -cwd -o ${log_path} -j y
done
