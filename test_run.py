
import argparse
import time
import resource
import os, sys, pandas, numpy, pathlib
import CytoSig
import pandas as pd
import gzip
import pickle

input_file = '/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/0.Tres_data/sc_cohorts/Nasopharyngeal.GSE162025.10x.pickle.gz'
with open(input_file, 'rb') as f:
    loaded_data = pickle.load(f)
    print("aaaa")
expression = pd.read_pickle(input_file)

