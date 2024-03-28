dataset_list=("Zhang2021" "SadeFeldman2018" "Yost2019" "Fraietta2018")
# get_bulk_profile
for dataset in ${dataset_list[*]}
do
    log_directory=/sibcb2/bioinformatics2/hongyuyang/code/Tres/log/3.clinical_validation
    if [ ! -d ${log_directory} ]; then
      mkdir ${log_directory}
    fi
    log_filename=${log_directory}/get_bulk_profile.${dataset}.log
    rm ${log_filename}

    echo "python3 /sibcb2/bioinformatics2/hongyuyang/code/Tres/3.clinical_validation/get_bulk_profile.py -D ${dataset}" | \
      qsub -q b1.q -N ${dataset} -V -cwd -o ${log_filename}  -j y
done

# correlation_calculate
for dataset in ${dataset_list[*]}
do
    log_directory=/sibcb2/bioinformatics2/hongyuyang/code/Tres/log/3.clinical_validation
    if [ ! -d ${log_directory} ]; then
      mkdir ${log_directory}
    fi
    log_filename=${log_directory}/correlation_calculate.${dataset}.log
    rm ${log_filename}

    echo "python3 /sibcb2/bioinformatics2/hongyuyang/code/Tres/3.clinical_validation/correlation_calculate.py -D ${dataset}" | \
      qsub -q g5.q@fnode010.sibcb.ac.cn -N ${dataset} -V -cwd -o ${log_filename}  -j y
done

