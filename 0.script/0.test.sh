log_path=/sibcb2/bioinformatics2/hongyuyang/code/Tres/log/test/pre_neutrophil2.log
rm log_path
echo "python3 /sibcb2/bioinformatics2/hongyuyang/code/Tres/1.data_pre/pre_neutrophil.py" | qsub -q g5.q -N pre_neutrophil -V -cwd -o ${log_path} -j y

log_path=/sibcb2/bioinformatics2/hongyuyang/code/Tres/log/test/test.R.log
rm log_path
echo "/sibcb/program/install/r-4.1/bin/Rscript /sibcb2/bioinformatics2/hongyuyang/code/Tres/1.data_pre/test.R" | qsub -q g5.q@fnode013.sibcb.ac.cn  -N test -V -cwd -o ${log_path} -j y