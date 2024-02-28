h5_list=$(find /sibcb2/bioinformatics2/hongyuyang/dataset/Tres/0.raw_data/tisch2 -name "*.h5")

for h5_path in ${h5_list}
do
  h5_filename=$(basename "${h5_path}")
  h5_tag=$(echo "${h5_filename}" | sed 's/_expression.h5//g')
  echo "Processing file: ${h5_tag}"

  celltype_mapping_rules_file=/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/0.model_data/tisch2_celltype_mapping_rule.txt
  output_file_directory=/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.gem_data/tisch2

  echo "python3 /sibcb2/bioinformatics2/hongyuyang/code/Tres/pre_process/convert_mtx.py -I ${h5_path} -R ${celltype_mapping_rules_file} -D ${output_file_directory} -O ${h5_tag}" | \
    qsub -N ${h5_tag} -V -cwd -o /sibcb2/bioinformatics2/hongyuyang/code/Tres/log/convert_mtx/${h5_tag}.log -j y
done
