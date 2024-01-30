echo "python3 /sibcb2/bioinformatics2/hongyuyang/code/Tres/csn/get_csndm.py" | \
qsub -q gpu.q -N get_csndm -V -cwd -o /sibcb2/bioinformatics2/hongyuyang/code/Tres/log/test/get_csndm.log -j y