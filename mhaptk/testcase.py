from mhaptk import main
import argparse, sys

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    args = parser.parse_args()
    args.metrics = " MM PDR CHALM MHL MCR MBS Entropy"
    args.cpgPath = "/sibcb1/bioinformatics/hongyuyang/code/Tres/mhaptk/testcase1.mhap.gz"
    args.cpgPath = "/sibcb1/bioinformatics/hongyuyang/code/Tres/mhaptk/hg38_CpG.gz"
    args.bedPath = "/sibcb1/bioinformatics/hongyuyang/code/Tres/mhaptk/MHB.bed"
    args.outputFile = "/sibcb1/bioinformatics/hongyuyang/code/Tres/mhaptk/testcase1.stat.tsv"
    main(args)