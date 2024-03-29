h5_list=$(find /sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch2_data/0.raw_data -name "*.h5")

h5_path=/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch2_data/0.raw_data/BRCA_GSE161529_expression.h5
for h5_path in ${h5_list[@]}
do
  h5_filename=$(basename "${h5_path}")
  h5_tag=$(echo "${h5_filename}" | sed 's/_expression.h5//g')
  echo "Processing file: ${h5_tag}"

  celltype_mapping_rules_file=/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/0.model_file/tisch2_celltype_mapping_rule.txt
  output_file_directory=/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch2_data/1.gem_data

  log_path=/sibcb2/bioinformatics2/hongyuyang/code/Tres/log/2.tisch2_data/convert_mtx/${h5_tag}.log
  rm ${log_path}
  echo "python3 /sibcb2/bioinformatics2/hongyuyang/code/Tres/1.data_pre/convert_mtx.py -I ${h5_path} -CTR ${celltype_mapping_rules_file} -D ${output_file_directory} -O ${h5_tag}" | \
    qsub -q b1.q@fnode004.sibcb.ac.cn -N ${h5_tag} -V -cwd -o ${log_path} -j y
  sleep 1m
done

