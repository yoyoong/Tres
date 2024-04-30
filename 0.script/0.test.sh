log_path=/sibcb2/bioinformatics2/hongyuyang/code/Tres/log/test/pre_neutrophil2.log
rm log_path
echo "python3 /sibcb2/bioinformatics2/hongyuyang/code/Tres/1.data_pre/pre_neutrophil.py" | qsub -q g5.q -N pre_neutrophil -V -cwd -o ${log_path} -j y

log_path=/sibcb2/bioinformatics2/hongyuyang/code/Tres/log/test/test.R.log
rm log_path
echo "/sibcb/program/install/r-4.1/bin/Rscript /sibcb2/bioinformatics2/hongyuyang/code/Tres/1.data_pre/test.R" | qsub -q g5.q@fnode013.sibcb.ac.cn  -N test -V -cwd -o ${log_path} -j y

log_path=/sibcb2/bioinformatics2/hongyuyang/code/Tres/log/test/mtres_response.log
rm log_path
echo "python3 /sibcb2/bioinformatics2/hongyuyang/code/Tres/2.mTres/mtres_response.py" | qsub -q g5.q@fnode010.sibcb.ac.cn -N mtres_response -V -cwd -o ${log_path} -j y

log_path=/sibcb2/bioinformatics2/hongyuyang/code/Tres/log/test/mtres_signaling.log
rm log_path
echo "python3 /sibcb2/bioinformatics2/hongyuyang/code/Tres/2.mTres/mtres_signaling.py" | qsub -q g5.q@fnode003.sibcb.ac.cn -N mtres_response -V -cwd -o ${log_path} -j y

for i in {1..5}
do
  expression_path=/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/1.neutrophil_data/Gao2024/Gao2024.Neutrophils.csv
  outdir=/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/1.neutrophil_data/Gao2024
  response_path=/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/1.neutrophil_data/Gao2024/Gao2024.signaling.csv
  signaling_path=/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/1.neutrophil_data/Gao2024/Gao2024.signaling.csv
  log_path=/sibcb2/bioinformatics2/hongyuyang/code/Tres/log/test/mtres_interaction.log
  output_tag=Gao2024.interaction

  if [ -f ${response_path} ] && [ -f ${signaling_path} ]; then
    if [ ! -f ${log_path} ]; then
      rm ${log_path}
      echo "python3 /sibcb2/bioinformatics2/hongyuyang/code/Tres/2.mTres/mtres_interaction.py -E ${expression_path} -R ${response_path} -S ${signaling_path} -D ${outdir} -O ${output_tag}" | \
        qsub -q g5.q@fnode010.sibcb.ac.cn -N mtres_interaction -V -cwd -o ${log_path} -j y
    fi
  fi

  echo "Processing: ${i}"
  sleep 1h
done