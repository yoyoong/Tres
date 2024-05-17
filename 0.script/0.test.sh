log_path=/sibcb2/bioinformatics2/hongyuyang/code/Tres/log/test/mtres_response.log
rm $log_path
echo "python3 /sibcb2/bioinformatics2/hongyuyang/code/Tres/2.mTres/mtres_response.py" | qsub -q b1.q@fnode006.sibcb.ac.cn -N mtres_response -V -cwd -o ${log_path} -j y

log_path=/sibcb2/bioinformatics2/hongyuyang/code/Tres/log/test/mtres_signaling.log
rm $log_path
echo "python3 /sibcb2/bioinformatics2/hongyuyang/code/Tres/2.mTres/mtres_signaling.py" | qsub -q g5.q@fnode003.sibcb.ac.cn -N mtres_response -V -cwd -o ${log_path} -j y

for i in {1..5}
do
  expression_path=/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.neutrophil_data/Gao2024.Neutrophils.csv

  outdir=/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.neutrophil_data
  response_path=/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.neutrophil_data/Gao2024.response.csv
  signaling_path=/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.neutrophil_data/Gao2024.signaling.csv
  log_path=/sibcb2/bioinformatics2/hongyuyang/code/Tres/log/test/mtres_interaction.log
  output_tag=Gao2024.interaction

  if [ -f ${response_path} ] && [ -f ${signaling_path} ]; then
    if [ ! -f ${log_path} ]; then
      rm ${log_path}
      echo "python3 /sibcb2/bioinformatics2/hongyuyang/code/Tres/2.mTres/mtres_interaction.py -E ${expression_path} -R ${response_path} -S ${signaling_path} -D ${outdir} -O ${output_tag}" | \
        qsub -q b1.q@fnode006.sibcb.ac.cn -N mtres_interaction -V -cwd -o ${log_path} -j y
    fi
  fi

  echo "Processing: ${i}"
  sleep 1h
done

log_path=/sibcb2/bioinformatics2/hongyuyang/code/Tres/log/test/interaction_cytokine_summary.log
rm ${log_path}
echo "python3 /sibcb2/bioinformatics2/hongyuyang/code/Tres/2.mTres/interaction_cytokine_summary.py" | qsub -q b1.q@fnode007.sibcb.ac.cn -N cytokine_summary -V -cwd -o ${log_path} -j y